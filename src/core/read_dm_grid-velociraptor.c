#include <hdf5_hl.h>
#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <math.h>

#define MIN(i, j) ((i) < (j) ? (i) : (j))

int read_dm_grid__velociraptor(
    int snapshot,
    float* slab)
{
    // N.B. We assume in this function that the slab has the fftw3 inplace complex dft padding.

    run_params_t* params = &(run_globals.params);
    int mpi_size = run_globals.mpi_size;
    int mpi_rank = run_globals.mpi_rank;

    // Have we read this slab before?
    if ((params->FlagInteractive || params->FlagMCMC) && !load_cached_deltax_slab(slab, snapshot))
        return 0;

    // read in the number of x values and offsets from every grid file
    // nx : number of x-dim values
    // ix_start : first x index
    // n_cell : number of values in each dim
    int *file_nx = NULL;
    int *file_ix_start = NULL;
    int file_n_cell[3] = { 0, 0, 0 };
    double box_size;
    int n_files = 999;
    const char fname_base[STRLEN] = {"%s/grids/snapshot_%03d.den.%d"};

    {
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);

        for (int ii = 0; ii < n_files; ii++) {
            char fname[STRLEN];
            sprintf(fname, fname_base, params->SimulationDir, snapshot, ii);
            hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);

            if (ii == 0) {
                herr_t status = H5LTget_attribute_int(file_id, "/", "Num_files", &n_files);
                assert(status >= 0);

                file_nx = calloc(n_files, sizeof(int));
                file_ix_start = calloc(n_files, sizeof(int));

                status = H5LTget_attribute_double(file_id, "/", "BoxSize", &box_size);
                assert(status >= 0);

                status = H5LTget_attribute_int(file_id, "/", "Ngrid_X", &file_n_cell[0]);
                assert(status >= 0);
                status = H5LTget_attribute_int(file_id, "/", "Ngrid_Y", &file_n_cell[1]);
                assert(status >= 0);
                status = H5LTget_attribute_int(file_id, "/", "Ngrid_Z", &file_n_cell[2]);
                assert(status >= 0);
            }

            herr_t status = H5LTget_attribute_int(file_id, "/", "Local_x_start", &file_ix_start[ii]);
            assert(status >= 0);
            status = H5LTget_attribute_int(file_id, "/", "Local_nx", &file_nx[ii]);
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
        for(int ii=0; ii < mpi_size; ii++) {
            recvcounts[ii] = sizeof(ptrdiff_t);
            displs[ii] = ii*sizeof(ptrdiff_t);
        }
        MPI_Allgatherv(rank_nx + mpi_rank, 1, MPI_BYTE, rank_nx, recvcounts, displs, MPI_BYTE, run_globals.mpi_comm);
        MPI_Allgatherv(rank_ix_start + mpi_rank, 1, MPI_BYTE, rank_ix_start, recvcounts, displs, MPI_BYTE, run_globals.mpi_comm);
        MPI_Allgatherv(rank_nI + mpi_rank, 1, MPI_BYTE, rank_nI, recvcounts, displs, MPI_BYTE, run_globals.mpi_comm);
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
    double* local_buffer = calloc(rank_nx[mpi_rank], sizeof(double));

    // read in the data
    // loop through each file and work out what cores are needed
    for (int ii = 0; ii < n_files; ii++) {
        int n_required_ranks = 0;
        bool rank_used = false;

        int required_ranks[mpi_size];
        for(int jj=0; jj < mpi_size; jj++)
            required_ranks[jj] = -1;

        for (int jj = 0; jj < mpi_size; jj++) {
            if ((rank_ix_start[jj] < (file_ix_start[ii] + file_nx[ii])) &&
                    (file_ix_start[ii] < rank_ix_start[jj] + rank_nx[jj])) {
                required_ranks[n_required_ranks++] = jj;
                if (jj == mpi_rank)
                    rank_used = true;
            }
        }
        if (mpi_rank == 0) {
            printf("file %d -> required_ranks = [ ", ii);
            for (int jj = 0; jj < n_required_ranks; jj++) {
                printf("%d ", required_ranks[jj]);
            }
            printf("]\n------------\n");
        }

        // create an mpi communicator with the required ranks
        if (rank_used) {
            MPI_Group file_group;
            MPI_Group_incl(run_group, n_required_ranks, required_ranks,
                    &file_group);

            MPI_Comm file_comm;
            MPI_Comm_create_group(MPI_COMM_WORLD, file_group, ii, &file_comm);

            // TODO(tidy): there must be a tidier work out these indices
            int file_start = 0;
            int rank_start = 0;
            int ix_diff = rank_ix_start[mpi_rank] - file_ix_start[ii];
            if (ix_diff >= 0) {
                file_start = ix_diff;
            } else {
                rank_start = -ix_diff;
            }
            int nx =
                MIN(file_nx[ii] - file_start, rank_nx[mpi_rank] - rank_start);

            printf("rank %d: file_start = %d, rank_start = %d, nx = %d\n",
                   mpi_rank, file_start, rank_start, nx);

            // select a hyperslab in the filespace
            hsize_t file_dims[3] = {(hsize_t)file_nx[ii], 16, 16};
            hid_t fspace_id = H5Screate_simple(3, file_dims, NULL);
            hsize_t start[3] = {(hsize_t)file_start, 0, 0};
            hsize_t count[3] = {(hsize_t)nx, 16, 16};
            H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count,
                                NULL);

            // create the memspace
            hsize_t mem_dims[3] = {(hsize_t)rank_nx[mpi_rank], 16, 16};
            hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);
            start[0] = (hsize_t)rank_start;
            H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, start, NULL, count,
                                NULL);

            hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, file_comm, MPI_INFO_NULL);
            char fname[128];
            sprintf(fname, fname_base, params->SimulationDir, snapshot, ii);

            hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
            H5Pclose(plist_id);

            hid_t dset_id = H5Dopen(file_id, "Density", H5P_DEFAULT);

            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5D_MPIO_COLLECTIVE);
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
    for (int ii = (file_nx[mpi_rank] - 1); ii >= 0; ii--)
        for (int jj = file_n_cell[1] - 1; jj >= 0; jj--)
            for (int kk = file_n_cell[2] - 1; kk >= 0; kk--)
                ((float*)rank_slab)[grid_index(ii, jj, kk, file_n_cell[1], INDEX_PADDED)] = (float)(local_buffer[grid_index(ii, jj, kk, file_n_cell[1], INDEX_REAL)]);

    free(local_buffer);
    free(file_ix_start);
    free(file_nx);

    // smooth the grid if needed
    smooth_grid(resample_factor, file_n_cell, rank_slab, rank_nI[mpi_rank], rank_ix_start[mpi_rank], rank_nx[mpi_rank]);

    // Copy the read and smoothed slab into the padded fft slab (already allocated externally)
    // TODO(simon): Cont here...
    ptrdiff_t slab_nix = run_globals.reion_grids.slab_nix[mpi_rank];
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[mpi_rank];
    for (int ii = 0; ii < slab_n_complex * 2; ii++)
        slab[ii] = 0.0;

    int ReionGridDim = run_globals.params.ReionGridDim;
    int n_every = file_n_cell[1] / ReionGridDim;
    for (int ii = 0; ii < slab_nix; ii++) {
        int i_hr = n_every * ii;
        assert((i_hr > -1) && (i_hr < rank_nx[mpi_rank]));
        for (int jj = 0; jj < ReionGridDim; jj++) {
            int j_hr = n_every * jj;
            assert((j_hr > -1) && (j_hr < file_n_cell[1]));
            for (int kk = 0; kk < ReionGridDim; kk++) {
                int k_hr = n_every * kk;
                assert((k_hr > -1) && (k_hr < file_n_cell[2]));

                slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] = ((float*)rank_slab)[grid_index(i_hr, j_hr, k_hr, file_n_cell[1], INDEX_PADDED)];
            }
        }
    }

    fftwf_free(rank_slab);

    // N.B. Hubble factor below to account for incorrect units in input DM grids!
    // TODO(pascal): Discuss this with Pascal and check carefully
    double mean = (double)run_globals.params.NPart * run_globals.params.PartMass / pow(box_size, 3) / run_globals.params.Hubble_h / run_globals.params.Hubble_h;

    // At this point grid holds the summed densities in each LR cell
    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    for (int ii = 0; ii < slab_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++) {
                float* val = &(slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]);
                // the fmax check here tries to account for negative densities introduced by fftw rounding / aliasing effects
                *val = fmaxf((float)(((double)*val / mean) - 1.), -1.0 + REL_TOL);
            }


    // Do we need to cache this slab?
    if (params->FlagInteractive || params->FlagMCMC)
        cache_deltax_slab(slab, snapshot);

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

    return 0;
}
