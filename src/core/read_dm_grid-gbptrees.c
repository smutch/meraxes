#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <math.h>

/*
   ==============================================================================
   MAJOR CODE REVISION by Paul Geil (October 2014)
   ==============================================================================

   - A significant numerical accuracy bug was resolved by performing calculations
    using an array of doubles and then casting it as an array of floats (as
    required by 21cmfast)

   - The TIAMAT velocity grids have not been analysed for anomalies.
 */

static inline void read_identifier(FILE* fin, bool skip_flag)
{
    char identifier[32];

    fread(identifier, sizeof(identifier), 1, fin);

    if (skip_flag)
        mlog_error("Skipping grid: %s...", identifier);
    else
        mlog("Reading grid: %s...", MLOG_MESG, identifier);
}

int read_dm_grid__gbptrees(
    int snapshot,
    float* slab)
{
    // N.B. We assume in this function that the slab has the fftw3 inplace complex dft padding.

    run_params_t* params = &(run_globals.params);

    // Have we read this slab before?
    if ((params->FlagInteractive || params->FlagMCMC) && !load_cached_deltax_slab(slab, snapshot))
        return 0;

    char fname[512];
    int n_cell[3];
    double box_size[3];
    int n_grids;
    int ma_scheme;
    long start_foffset;
    int ReionGridDim = run_globals.params.ReionGridDim;

    // Construct the input filename
    sprintf(fname, "%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, snapshot);

    // Read the header
    if (run_globals.mpi_rank == 0) {
        FILE* fd;
        if ((fd = fopen(fname, "rb")) == NULL) {
            fprintf(stderr, "Failed to open file: %s\n", fname);
            return EXIT_FAILURE;
        }

        fread(n_cell, sizeof(int), 3, fd);
        fread(box_size, sizeof(double), 3, fd);
        fread(&n_grids, sizeof(int), 1, fd);
        fread(&ma_scheme, sizeof(int), 1, fd);

        mlog("Reading grid for snapshot %d", MLOG_OPEN | MLOG_TIMERSTART, snapshot);
        mlog("n_cell = [%d, %d, %d]", MLOG_MESG, n_cell[0], n_cell[1], n_cell[2]);
        mlog("box_size = [%.2f, %.2f, %.2f] cMpc/h", MLOG_MESG, box_size[0], box_size[1], box_size[2]);
        mlog("ma_scheme = %d", MLOG_MESG, ma_scheme);

        if (n_grids != 4) {
            mlog_error("n_grids != 4 as expected...");
            fclose(fd);
            ABORT(EXIT_FAILURE);
        }

        assert((n_cell[0] == n_cell[1]) && (n_cell[1] == n_cell[2])
            && "Input grids are not cubic!");

        // Note that we are expecting the first grid to be the density grid
        read_identifier(fd, false);

        start_foffset = ftell(fd);
        fclose(fd);
    }

    // share the needed information with all ranks
    MPI_Bcast(n_cell, 3, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(box_size, 3, MPI_DOUBLE, 0, run_globals.mpi_comm);
    MPI_Bcast(&start_foffset, 1, MPI_LONG, 0, run_globals.mpi_comm);

    // Check if the grid in the file is higher resolution than we require
    double resample_factor = calc_resample_factor(n_cell);

    // Malloc the slab
    ptrdiff_t slab_nix = run_globals.reion_grids.slab_nix[run_globals.mpi_rank];
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

    ptrdiff_t slab_nix_file, slab_ix_start_file;
    ptrdiff_t slab_n_complex_file = fftwf_mpi_local_size_3d(n_cell[0], n_cell[0], n_cell[0] / 2 + 1, run_globals.mpi_comm, &slab_nix_file, &slab_ix_start_file);
    fftwf_complex* slab_file = fftwf_alloc_complex((size_t)slab_n_complex_file);
    ptrdiff_t slab_ni_file = slab_nix_file * n_cell[0] * n_cell[0];

    // Initialise (just in case!)
    for (int ii = 0; ii < slab_n_complex_file; ii++)
        slab_file[ii] = 0 + 0 * I;
    // N.B. factor of two for fftw padding
    for (int ii = 0; ii < slab_n_complex * 2; ii++)
        slab[ii] = 0.0;

    // Read in the slab for this rank
    MPI_File fin = NULL;
    MPI_Status status;
    MPI_Offset slab_offset = start_foffset / sizeof(float) + (slab_ix_start_file * n_cell[1] * n_cell[2]);

    MPI_File_open(run_globals.mpi_comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
    MPI_File_set_view(fin, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

    ptrdiff_t chunk_size = slab_ni_file;
    int n_reads = 1;
    while (chunk_size > INT_MAX / 16) {
        chunk_size /= 2;
        n_reads *= 2;
    }

    MPI_Offset offset = slab_offset;
    for (int ii = 0; ii < n_reads; ii++, offset += chunk_size) {
        MPI_File_read_at(fin, offset, &(((float*)slab_file)[chunk_size * ii]), (int)chunk_size, MPI_FLOAT, &status);

        int count_check;
        MPI_Get_count(&status, MPI_FLOAT, &count_check);
        if (count_check != chunk_size) {
            mlog_error("Failed to read correct number of elements on rank %d.", run_globals.mpi_rank);
            mlog_error("Expected %d but read %d.", (int)chunk_size, count_check);
            ABORT(EXIT_FAILURE);
        }
    }
    MPI_File_close(&fin);

    // reorder the read slab for inplace fftw padding
    for (int ii = (int)(slab_nix_file - 1); ii >= 0; ii--)
        for (int jj = n_cell[0] - 1; jj >= 0; jj--)
            for (int kk = n_cell[0] - 1; kk >= 0; kk--)
                ((float*)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_PADDED)] = ((float*)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_REAL)];

    // smooth the grid if needed
    smooth_grid(resample_factor, n_cell, slab_file, slab_n_complex_file, slab_ix_start_file, slab_nix_file, snapshot);

    // Glossary
    // ix_hi: the x-axis index on the high res grid read from the input file
    // ix_lo: the x-axis index on the low res, subsampled, grid which will be used for the reionisation calculation
    // hi_rank: the rank which holds the current ix_hi value
    // lo_rank: the rank which will hold the current ix_lo value
    // slice: a single x-index cut through an array. i.e. one x-value and all of the y and z values

    // TODO: This needs tidied...  There are name collisions and some variables should be moved from above down here...
    // TODO: Much of this should be applicable to VR trees. I should pull this out in to a separate function.
    // gather all of the slab_ix_start_file values
    int* ix_hi_start_allranks = calloc(run_globals.mpi_size, sizeof(int));  // TODO: Could this be on the stack instead?
    int* nix_hi_allranks = calloc(run_globals.mpi_size, sizeof(int));
    int ix_hi_start = (int)slab_ix_start_file;
    int nix_hi = (int)slab_nix_file;

    MPI_Allgather(&ix_hi_start, 1, MPI_INT, ix_hi_start_allranks, 1, MPI_INT, run_globals.mpi_comm);
    MPI_Allgather(&nix_hi, 1, MPI_INT, nix_hi_allranks, 1, MPI_INT, run_globals.mpi_comm);

    int n_every = n_cell[0] / ReionGridDim;
    int slice_size = ReionGridDim * ReionGridDim;
    float* slice = calloc(slice_size, sizeof(float));
    int nix_lo = (int)run_globals.reion_grids.slab_nix[run_globals.mpi_rank];

    int ix_lo_start = (int)run_globals.reion_grids.slab_ix_start[0];
    int ix_lo_end = ix_lo_start + (int)run_globals.reion_grids.slab_nix[0] - 1;
    ix_hi_start = ix_hi_start_allranks[0];
    int ix_hi_end = ix_hi_start + nix_hi_allranks[0] - 1;

    for (int ix_hi = 0, ix_lo = 0, hi_rank = 0, lo_rank = 0; ix_hi < n_cell[0]; ix_hi += n_every, ++ix_lo) {

        if (ix_lo > ix_lo_end) lo_rank++;
        if (ix_hi > ix_hi_end) hi_rank++;

        ix_lo_start = (int)run_globals.reion_grids.slab_ix_start[lo_rank];
        ix_lo_end = ix_lo_start + (int)run_globals.reion_grids.slab_nix[lo_rank] - 1;
        ix_hi_start = ix_hi_start_allranks[hi_rank];
        ix_hi_end = ix_hi_start + nix_hi_allranks[hi_rank] - 1;

        int ix_lo_slab = ix_lo - ix_lo_start;
        int ix_hi_slab = ix_hi - ix_hi_start;

        if ((hi_rank == run_globals.mpi_rank) && (lo_rank != run_globals.mpi_rank)) {
            // This rank has the high res x-value but needs to send it somewhere else...

            // pack the slice
            for (int iy_lo = 0; iy_lo < ReionGridDim; iy_lo++) {
                int iy_hi = n_every * iy_lo;
                assert((iy_hi > -1) && (iy_hi < n_cell[0]));

                for (int iz_lo = 0; iz_lo < ReionGridDim; iz_lo++) {
                    int iz_hi = n_every * iz_lo;
                    assert((iz_hi > -1) && (iz_hi < n_cell[0]));

                    slice[grid_index(0, iy_lo, iz_lo, ReionGridDim, INDEX_REAL)] = ((float*)slab_file)[grid_index(ix_hi_slab, iy_hi, iz_hi, n_cell[0], INDEX_PADDED)];
                }

            }

            // send it to the lo rank
            int ident = 400000 + 1000*run_globals.mpi_rank + lo_rank;
            MPI_Send(slice, slice_size, MPI_FLOAT, lo_rank, ident, run_globals.mpi_comm);

        } else if ((hi_rank != run_globals.mpi_rank) && (lo_rank == run_globals.mpi_rank)) {
            // This rank needs the high res x-value but needs to get it from somewhere else...

            // receive the packed slice
            int ident = 400000 + 1000*hi_rank + run_globals.mpi_rank;
            MPI_Recv(slice, slice_size, MPI_FLOAT, hi_rank, ident, run_globals.mpi_comm, MPI_STATUS_IGNORE);

            assert((ix_lo_slab >= 0) && (ix_lo_slab < nix_lo));

            // store it
            for (int iy_lo = 0; iy_lo < ReionGridDim; iy_lo++) {
                int iy_hi = n_every * iy_lo;
                assert((iy_hi > -1) && (iy_hi < n_cell[0]));

                for (int iz_lo = 0; iz_lo < ReionGridDim; iz_lo++) {
                    int iz_hi = n_every * iz_lo;
                    assert((iz_hi > -1) && (iz_hi < n_cell[0]));

                    slab[grid_index(ix_lo_slab, iy_lo, iz_lo, ReionGridDim, INDEX_PADDED)] = slice[grid_index(0, iy_lo, iz_lo, ReionGridDim, INDEX_REAL)];
                }

            }

        } else if ((hi_rank == run_globals.mpi_rank) && (lo_rank == run_globals.mpi_rank)){
            // This rank needs the high res x-value and already has it...

            assert((ix_lo_slab >= 0) && (ix_lo_slab < nix_lo));

            for (int iy_lo = 0; iy_lo < ReionGridDim; iy_lo++) {
                int iy_hi = n_every * iy_lo;
                assert((iy_hi > -1) && (iy_hi < n_cell[0]));

                for (int iz_lo = 0; iz_lo < ReionGridDim; iz_lo++) {
                    int iz_hi = n_every * iz_lo;
                    assert((iz_hi > -1) && (iz_hi < n_cell[0]));

                    slab[grid_index(ix_lo_slab, iy_lo, iz_lo, ReionGridDim, INDEX_PADDED)] = ((float*)slab_file)[grid_index(ix_hi_slab, iy_hi, iz_hi, n_cell[0], INDEX_PADDED)];
                }

            }

        }

        // If we didn't match any of the above cases, we neither have nor want theis slice.
    }

    free(slice);
    free(nix_hi_allranks);
    free(ix_hi_start_allranks);

    // N.B. Hubble factor below to account for incorrect units in input DM grids!
    float mean_inv = pow(box_size[0], 3) * run_globals.params.Hubble_h / run_globals.params.NPart / run_globals.params.PartMass;

    // At this point grid holds the summed densities in each LR cell
    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    for (int ii = 0; ii < slab_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++) {
                float* val = &(slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]);
                // the fmax check here tries to account for negative densities introduced by fftw rounding / aliasing effects
                *val = fmaxf(*val * mean_inv - 1.0, -1.0);
            }

    fftwf_free(slab_file);

    // Do we need to cache this slab?
    if (params->FlagInteractive || params->FlagMCMC)
        cache_deltax_slab(slab, snapshot);

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

    return 0;
}
