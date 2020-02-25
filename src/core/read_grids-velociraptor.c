#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <gsl/gsl_sort_int.h>
#include <hdf5_hl.h>
#include <math.h>

#define MIN(i, j) ((i) < (j) ? (i) : (j))

// WARNING:  This is a hack introduced to fix the reading of VELOCIraptor
// velocity grid files which have no `Num_files` attribute.  It relies on the
// density files always having been read first to set this static variable
// which is then used for reading the velocity files.
static int vr_num_files_hack_ = 0;

int read_grid__velociraptor(
    const enum grid_prop property,
    const int snapshot,
    float* slab)
{
    // N.B. We assume in this function that the slab has the fftw3 inplace complex dft padding.

    run_params_t* params = &(run_globals.params);
    int mpi_size = run_globals.mpi_size;
    int mpi_rank = run_globals.mpi_rank;

    // Have we read this slab before?
    if ((params->FlagInteractive || params->FlagMCMC) && !load_cached_slab(slab, snapshot, property))
        return 0;

    if((property == X_VELOCITY) || (property == Y_VELOCITY) || (property == Z_VELOCITY)) {

        if(run_globals.params.TsVelocityComponent < 1 || run_globals.params.TsVelocityComponent > 3) {
            mlog("Not a valid velocity direction: 1 - x, 2 - y, 3 - z", MLOG_MESG);
            ABORT(EXIT_FAILURE);
        }

        if(run_globals.params.Flag_ConstructLightcone && run_globals.params.TsVelocityComponent!=3) {
            mlog("Light-cone is generated along the z-direction, therefore the velocity component should be in the z-direction (i.e 3).", MLOG_MESG);
            ABORT(EXIT_FAILURE);
        }

    }

    // read in the number of x values and offsets from every grid file
    // nx : number of x-dim values
    // ix_start : first x index
    // n_cell : number of values in each dim
    int* file_nx = NULL;
    int* file_ix_start = NULL;
    int file_n_cell[3] = { 0, 0, 0 };
    double box_size;
    int n_files = 999;
    const char fname_base[STRLEN] = { "%s/grids/snapshot_%03d.%s.%d" };

    {
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);

        for (int ii = 0; ii < n_files; ii++) {

            char fname[STRLEN];
            switch (property){
                case X_VELOCITY:
                case Y_VELOCITY:
                case Z_VELOCITY:
                    sprintf(fname, fname_base, params->SimulationDir, snapshot, "vel", ii);
                    break;
                case DENSITY:
                    sprintf(fname, fname_base, params->SimulationDir, snapshot, "den", ii);
                    break;
                default:
                    mlog_error("Unrecognised grid property in read_grid__velociraptor!");
                    break;
            }

            hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);

            if (ii == 0) {

                // save current error stack
                herr_t (*old_func)(long, void*);
                void *old_client_data;
                hid_t error_stack = 0;
                H5Eget_auto(error_stack, &old_func, &old_client_data);

                // turn off error handling
                H5Eset_auto(error_stack, NULL, NULL);

                herr_t status = H5LTget_attribute_int(file_id, "/", "Num_files", &n_files);

                // restore error handling
                H5Eset_auto(error_stack, old_func, old_client_data);

                // WARNING: This is a hack!  See the `vr_num_files_hack_`
                // declaration above for details.
                if (status >= 0) {
                    vr_num_files_hack_ = n_files;
                } else {
                    assert(vr_num_files_hack_ > 0);
                    n_files = vr_num_files_hack_;
                }

                file_nx = calloc(n_files, sizeof(int));
                file_ix_start = calloc(n_files, sizeof(int));

                status = H5LTget_attribute_double(file_id, "/", "BoxSize", &box_size);
                assert(status >= 0);

                status = H5LTget_attribute_int(file_id, "/", "Ngrid_X", file_n_cell);
                assert(status >= 0);
                status = H5LTget_attribute_int(file_id, "/", "Ngrid_Y", file_n_cell + 1);
                assert(status >= 0);
                status = H5LTget_attribute_int(file_id, "/", "Ngrid_Z", file_n_cell + 2);
                assert(status >= 0);
            }

            herr_t status = H5LTget_attribute_int(file_id, "/", "Local_x_start", file_ix_start + ii);
            assert(status >= 0);
            status = H5LTget_attribute_int(file_id, "/", "Local_nx", file_nx + ii);
            assert(status >= 0);

            H5Fclose(file_id);
        }

        H5Pclose(plist_id);
    }

    assert((file_n_cell[0] == file_n_cell[1]) && (file_n_cell[1] == file_n_cell[2])
        && "Input grids are not cubic!");

    mlog("Reading VELOCIraptor grid for snapshot %d", MLOG_OPEN | MLOG_TIMERSTART, snapshot);
    mlog("n_cell = [%d, %d, %d]", MLOG_MESG, file_n_cell[0], file_n_cell[1], file_n_cell[2]);
    mlog("box_size = %.2f cMpc/h", MLOG_MESG, box_size * params->Hubble_h);

    double resample_factor = calc_resample_factor(file_n_cell);

    // Malloc the slab for this rank we need given the dimensionality of the input file grid.
    // nI : total number of complex values in slab on this rank
    // nR : total number of real values in slab on this rank
    ptrdiff_t rank_nx[mpi_size];
    ptrdiff_t rank_ix_start[mpi_size];
    ptrdiff_t rank_nI[mpi_size];
    rank_nI[mpi_rank] = fftwf_mpi_local_size_3d(file_n_cell[0], file_n_cell[1], file_n_cell[2] / 2 + 1, run_globals.mpi_comm, &rank_nx[mpi_rank], &rank_ix_start[mpi_rank]);

    {
        // Note: I'm using bytes here as I'm not sure what the equivaalent MPI dataype for a ptrdiff_t.
        int recvcounts[mpi_size];
        int displs[mpi_size];
        for (int ii = 0; ii < mpi_size; ii++) {
            recvcounts[ii] = sizeof(ptrdiff_t);
            displs[ii] = ii * sizeof(ptrdiff_t);
        }
        MPI_Allgatherv(&rank_nx[mpi_rank], 1, MPI_BYTE, rank_nx, recvcounts, displs, MPI_BYTE, run_globals.mpi_comm);
        MPI_Allgatherv(&rank_ix_start[mpi_rank], 1, MPI_BYTE, rank_ix_start, recvcounts, displs, MPI_BYTE, run_globals.mpi_comm);
        MPI_Allgatherv(&rank_nI[mpi_rank], 1, MPI_BYTE, rank_nI, recvcounts, displs, MPI_BYTE, run_globals.mpi_comm);
    }

    fftwf_complex* rank_slab = fftwf_alloc_complex((size_t)rank_nI[mpi_rank]);

    // Initialise (just in case!)
    for (int ii = 0; ii < rank_nI[mpi_rank]; ii++)
        rank_slab[ii] = 0 + 0 * I;

    MPI_Group run_group;
    MPI_Comm_group(MPI_COMM_WORLD, &run_group);

    // We are currently assuming the grids to be float, but the VELOCIraptor
    // grids are doubles.  For the moment, let's just read the doubles into a
    // buffer and change them to float appropriately.
    double* local_buffer = calloc((size_t)(rank_nx[mpi_rank] * file_n_cell[1] * file_n_cell[2]), sizeof(double));

    // loop through each file and work out what cores are needed
    int n_required_ranks[n_files];
    bool rank_used[n_files];
    int required_ranks[n_files * mpi_size];

#define rr_index(ii, jj) ((ii) * (mpi_size) + (jj))

    for (int ii = 0; ii < n_files; ii++) {

        n_required_ranks[ii] = 0;
        rank_used[ii] = false;

        for (int jj = 0; jj < mpi_size; jj++)
            required_ranks[rr_index(ii, jj)] = -1;

        for (int jj = 0; jj < mpi_size; jj++) {
            if ((rank_ix_start[jj] < (file_ix_start[ii] + file_nx[ii])) && (file_ix_start[ii] < (rank_ix_start[jj] + rank_nx[jj]))) {
                required_ranks[rr_index(ii, n_required_ranks[ii]++)] = jj;
                if (jj == mpi_rank)
                    rank_used[ii] = true;
            }
        }
    }

    // sort the files by n_required_ranks
    size_t sort_ind[n_files];
    gsl_sort_int_index(sort_ind, n_required_ranks, 1, (const size_t)n_files);

    char dset_name[32];
    switch (property){
        case X_VELOCITY:
            sprintf(dset_name, "Vx");
            break;
        case Y_VELOCITY:
            sprintf(dset_name, "Vy");
            break;
        case Z_VELOCITY:
            sprintf(dset_name, "Vz");
            break;
        case DENSITY:
            sprintf(dset_name, "Density");
            break;
        default:
            mlog_error("Unrecognised grid property in read_grid__velociraptor!");
            break;
    }

    for (int jj = 0; jj < n_files; jj++) {
        int ii = (int)sort_ind[jj];

        // read in the data
        // create an mpi communicator with the required ranks
        if (rank_used[ii]) {
            MPI_Group file_group;
            MPI_Group_incl(run_group, n_required_ranks[ii], required_ranks + rr_index(ii, 0),
                &file_group);

            MPI_Comm file_comm;
            MPI_Comm_create_group(MPI_COMM_WORLD, file_group, ii, &file_comm);

            // There must be a tidier work out these indices...
            int file_start = 0;
            int rank_start = 0;
            int ix_diff = (int)(rank_ix_start[mpi_rank] - file_ix_start[ii]);
            if (ix_diff >= 0) {
                file_start = ix_diff;
            } else {
                rank_start = -ix_diff;
            }
            int nx = (int)MIN(file_nx[ii] - file_start, rank_nx[mpi_rank] - rank_start);

            // select a hyperslab in the filespace
            hid_t fspace_id = H5Screate_simple(1, (hsize_t[1]) { (hsize_t)(file_nx[ii] * file_n_cell[1] * file_n_cell[2]) },
                NULL);
            H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET,
                (hsize_t[1]) { (hsize_t)(file_start * file_n_cell[1] * file_n_cell[2]) }, NULL,
                (hsize_t[1]) { (hsize_t)(nx * file_n_cell[1] * file_n_cell[2]) }, NULL);

            // create the memspace
            hid_t memspace_id = H5Screate_simple(1,
                (hsize_t[1]) { (hsize_t)(rank_nx[mpi_rank] * file_n_cell[1] * file_n_cell[1]) }, NULL);
            H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET,
                (hsize_t[1]) { (hsize_t)(rank_start * file_n_cell[1] * file_n_cell[2]) }, NULL,
                (hsize_t[1]) { (hsize_t)(nx * file_n_cell[1] * file_n_cell[2]) }, NULL);

            hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, file_comm, MPI_INFO_NULL);

            char fname[STRLEN];
            switch (property){
                case X_VELOCITY:
                case Y_VELOCITY:
                case Z_VELOCITY:
                    sprintf(fname, fname_base, params->SimulationDir, snapshot, "vel", ii);
                    break;
                case DENSITY:
                    sprintf(fname, fname_base, params->SimulationDir, snapshot, "den", ii);
                    break;
                default:
                    mlog_error("Unrecognised grid property in read_grid__velociraptor!");
                    break;
            }

            hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
            H5Pclose(plist_id);

            hid_t dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);

            plist_id = H5Pcreate(H5P_DATASET_XFER);

            // TODO(performance): Is a collective read better here?
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

            H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, plist_id, local_buffer);

            H5Pclose(plist_id);

            H5Dclose(dset_id);
            H5Fclose(file_id);
            H5Sclose(memspace_id);
            H5Sclose(fspace_id);

            MPI_Comm_free(&file_comm);
            MPI_Group_free(&file_group);
        }
    }

    MPI_Group_free(&run_group);

    // move the doubles into the float array, with inplace fftw padding
    for (int ii = 0; ii < rank_nx[mpi_rank]; ++ii)
        for (int jj = 0; jj < file_n_cell[1]; ++jj)
            for (int kk = 0; kk < file_n_cell[2]; ++kk)
                ((float*)rank_slab)[grid_index(ii, jj, kk, file_n_cell[1], INDEX_PADDED)] = (float)(local_buffer[grid_index(ii, jj, kk, file_n_cell[1], INDEX_REAL)]);

    free(local_buffer);
    free(file_ix_start);
    free(file_nx);

    // smooth the grid if needed
    smooth_grid(resample_factor, file_n_cell, rank_slab, rank_nI[mpi_rank], rank_ix_start[mpi_rank], rank_nx[mpi_rank]);
    subsample_grid(resample_factor, file_n_cell, (int)rank_ix_start[mpi_rank], (int)rank_nx[mpi_rank], (float*)rank_slab, slab);

    fftwf_free(rank_slab);

    if (property == DENSITY) {
        // TODO: Discuss this with Pascal and check carefully
        float mean_inv = pow(box_size, 3) * run_globals.params.Hubble_h / ((double)run_globals.params.NPart * run_globals.params.PartMass);

    // At this point grid holds the summed densities in each LR cell
    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    int ReionGridDim = run_globals.params.ReionGridDim;
    int slab_nix = run_globals.reion_grids.slab_nix[mpi_rank];
    for (int ii = 0; ii < slab_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++) {
                float* val = &(slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]);
                // the fmax check here tries to account for negative densities introduced by fftw rounding / aliasing effects
                *val = fmaxf(*val * mean_inv - 1.0, -1.0);
            }
    }

    // Do we need to cache this slab?
    if (params->FlagInteractive || params->FlagMCMC)
        cache_slab(slab, snapshot, property);

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

    return 0;
}
