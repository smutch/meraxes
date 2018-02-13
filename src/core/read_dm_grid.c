#include "meraxes.h"
#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

int load_cached_deltax_slab(float* slab, int snapshot)
{
    if (run_globals.SnapshotDeltax[snapshot] != NULL) {
        ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];
        memcpy(slab, run_globals.SnapshotDeltax[snapshot], sizeof(float) * slab_n_complex * 2);
        return 0;
    }
    else
        return 1;
}

int cache_deltax_slab(float* slab, int snapshot)
{
    if (run_globals.SnapshotDeltax[snapshot] == NULL) {
        float** cache = &run_globals.SnapshotDeltax[snapshot];
        ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];
        ptrdiff_t mem_size = sizeof(float) * slab_n_complex * 2;

        *cache = fftwf_alloc_real(mem_size);
        memcpy(*cache, slab, mem_size);
        return 0;
    }
    else
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
