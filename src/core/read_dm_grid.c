#include "meraxes.h"
#include <fftw3-mpi.h>

double calc_resample_factor(int n_cell[3])
{
    // Check if the grid in the file is higher resolution than we require
    int ReionGridDim = run_globals.params.ReionGridDim;
    if (n_cell[0] != ReionGridDim) {
        double resample_factor = (double)ReionGridDim / (double)n_cell[0];
        if (resample_factor > 1.0001) {
            mlog_error("Grid has a resolution less than that required! Aborting!");
            ABORT(EXIT_FAILURE);
        }
        mlog("Using resample factor = %.3f", MLOG_MESG, resample_factor);
        return resample_factor;
    } else
        return 1.0;
}

void smooth_grid(double resample_factor, int n_cell[3], fftwf_complex* slab, ptrdiff_t slab_n_complex, ptrdiff_t slab_ix_start, ptrdiff_t slab_nix)
{
    if (resample_factor < 1.0) {
        mlog("Smoothing hi-res grid...", MLOG_OPEN | MLOG_TIMERSTART);
        fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(n_cell[0], n_cell[1], n_cell[2], (float*)slab, slab, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
        // real space to k-space.
        // Note: we will leave off factor of VOLUME, in anticipation of the inverse
        // FFT below
        double total_n_cells = n_cell[0] * n_cell[1] * n_cell[2];
        for (int ii = 0; ii < (int)slab_n_complex; ii++)
            slab[ii] /= total_n_cells;

        filter(slab, (int)slab_ix_start, (int)slab_nix, n_cell[0], (float)(run_globals.params.BoxSize / (double)run_globals.params.ReionGridDim / 2.0), 0); // NOTE: Real space top-hat hard-coded for this

        plan = fftwf_mpi_plan_dft_c2r_3d(n_cell[0], n_cell[1], n_cell[2], slab, (float*)slab, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
    }
}

int load_cached_deltax_slab(float* slab, int snapshot)
{
    if (run_globals.SnapshotDeltax[snapshot] != NULL) {
        ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];
        memcpy(slab, run_globals.SnapshotDeltax[snapshot], sizeof(float) * slab_n_complex * 2);
        mlog("Loaded deltax slab from cache.", MLOG_MESG);
        return 0;
    } else
        return 1;
}

int cache_deltax_slab(float* slab, int snapshot)
{
    if (run_globals.SnapshotDeltax[snapshot] == NULL) {
        float** cache = &run_globals.SnapshotDeltax[snapshot];
        ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];
        ptrdiff_t mem_size = sizeof(float) * slab_n_complex * 2;

        *cache = fftwf_alloc_real((size_t)mem_size);
        memcpy(*cache, slab, mem_size);
        return 0;
    } else
        return 1;
}

void free_grids_cache()
{
    if (run_globals.params.Flag_PatchyReion) {
        float** snapshot_deltax = run_globals.SnapshotDeltax;

        if (run_globals.params.FlagInteractive)
            for (int ii = 0; ii < run_globals.NStoreSnapshots; ii++)
                fftwf_free(snapshot_deltax[ii]);

        free(snapshot_deltax);
    }
}
