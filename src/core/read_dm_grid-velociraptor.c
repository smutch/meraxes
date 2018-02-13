#include "meraxes.h"
#include "hdf5_hl.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <math.h>

int read_dm_grid__velociraptor(
    int snapshot,
    float* slab)
{
    // N.B. We assume in this function that the slab has the fftw3 inplace complex dft padding.

    run_params_t* params = &(run_globals.params);

    // Have we read this slab before?
    if (params->FlagInteractive && !load_cached_deltax_slab(slab, snapshot))
        return 0;

    // generate the input filename
    char fname[STRLEN];
    sprintf(fname, "%s/grids/snapshot_%03d.den.0", params->SimulationDir, snapshot);

    // open the file (in parallel)
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
    hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    herr_t status = -1;

    int n_cell[3] = {0, 0, 0};
    status = H5LTget_attribute_int(file_id, "/", "Ngrid_X", &n_cell[0]);
    assert(status >= 0);
    status = H5LTget_attribute_int(file_id, "/", "Ngrid_Y", &n_cell[1]);
    assert(status >= 0);
    status = H5LTget_attribute_int(file_id, "/", "Ngrid_Z", &n_cell[2]);
    assert(status >= 0);

    int local_x_start = -1;
    int local_nx;
    status = H5LTget_attribute_int(file_id, "/", "Local_x_start", &local_x_start);
    assert(status >= 0);
    assert(local_x_start == 0 && "We currently require Local_x_start = 0");
    status = H5LTget_attribute_int(file_id, "/", "Local_nx", &local_nx);
    assert(status >= 0);
    assert(local_nx == n_cell[0] && "We currently require Local_nx == Ngrid_X");

    double box_size;
    status = H5LTget_attribute_double(file_id, "/", "BoxSize", &box_size);
    assert(status >= 0);


    assert((n_cell[0] == n_cell[1]) && (n_cell[1] == n_cell[2])
            && "Input grids are not cubic!");

    double resample_factor = calc_resample_factor(n_cell);

    // Malloc the slab
    ptrdiff_t slab_nix = run_globals.reion_grids.slab_nix[run_globals.mpi_rank];
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

    ptrdiff_t slab_nix_file, slab_ix_start_file;
    ptrdiff_t slab_n_complex_file = fftwf_mpi_local_size_3d(n_cell[0], n_cell[1], n_cell[2] / 2 + 1, run_globals.mpi_comm, &slab_nix_file, &slab_ix_start_file);
    fftwf_complex* slab_file = fftwf_alloc_complex(slab_n_complex_file);
    ptrdiff_t slab_ni_file = slab_nix_file * n_cell[1] * n_cell[2];

    // Initialise (just in case!)
    for (int ii = 0; ii < slab_n_complex_file; ii++)
        slab_file[ii] = 0 + 0 * I;
    // N.B. factor of two for fftw padding
    for (int ii = 0; ii < slab_n_complex * 2; ii++)
        slab[ii] = 0.0;

    // read in the data
    // open the dataset 
    hid_t dset_id = H5Dopen(file_id, "Density", H5P_DEFAULT);

    // select a hyperslab in the filespace
    hsize_t file_dims[3] = { n_cell[0], n_cell[1], n_cell[2] };
    hid_t fspace_id = H5Screate_simple(3, file_dims, NULL);
    hsize_t start[3] = { slab_ix_start_file, 0, 0 };
    hsize_t count[3] = { slab_nix_file, n_cell[1], n_cell[2] };
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

    // create the memspace
    hsize_t mem_dims[3] = { slab_nix_file, n_cell[1], n_cell[2] };
    hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

    // We are currently assuming the grids to be float, but the VELOCIraptor
    // grids are doubles.  For the moment, let's just read the doubles into a
    // buffer and change them to float appropriately.
    double* real_buffer = malloc(sizeof(double) * slab_ni_file);
    for(int ii=0; ii < (int)slab_ni_file; ii++)
        real_buffer[ii] = 0.0;

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, plist_id, real_buffer); 
    H5Pclose(plist_id);

    H5Sclose(memspace_id);
    H5Sclose(fspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);

    // move the doubles into the float array, with inplace fftw padding
    for (int ii = slab_nix_file - 1; ii >= 0; ii--)
        for (int jj = n_cell[1] - 1; jj >= 0; jj--)
            for (int kk = n_cell[2] - 1; kk >= 0; kk--)
                ((float*)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_PADDED)] = (float)(real_buffer[grid_index(ii, jj, kk, n_cell[0], INDEX_REAL)]);


    free(real_buffer);

    // smooth the grid if needed
    smooth_grid(resample_factor, n_cell, slab_file, slab_n_complex_file, slab_ix_start_file, slab_nix_file);

    // Copy the read and smoothed slab into the padded fft slab (already allocated externally)
    int ReionGridDim = run_globals.params.ReionGridDim;
    int n_every = n_cell[0] / ReionGridDim;
    for (int ii = 0; ii < slab_nix; ii++) {
        int i_hr = n_every * ii;
        assert((i_hr > -1) && (i_hr < slab_nix_file));
        for (int jj = 0; jj < ReionGridDim; jj++) {
            int j_hr = n_every * jj;
            assert((j_hr > -1) && (j_hr < n_cell[0]));
            for (int kk = 0; kk < ReionGridDim; kk++) {
                int k_hr = n_every * kk;
                assert((k_hr > -1) && (k_hr < n_cell[0]));

                slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] = ((float*)slab_file)[grid_index(i_hr, j_hr, k_hr, n_cell[0], INDEX_PADDED)];
            }
        }
    }

    // N.B. Hubble factor below to account for incorrect units in input DM grids!
    // TODO: Discuss this with Pascal
    double mean = (double)run_globals.params.NPart * run_globals.params.PartMass / pow(box_size, 3) / run_globals.params.Hubble_h;

    // At this point grid holds the summed densities in each LR cell
    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    for (int ii = 0; ii < slab_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++) {
                float* val = &(slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]);
                *val = (float)(((double)*val / mean) - 1.);
            }

    fftwf_free(slab_file);

    // Do we need to cache this slab?
    if (params->FlagInteractive)
        cache_deltax_slab(slab, snapshot);

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

    return 0;
}
