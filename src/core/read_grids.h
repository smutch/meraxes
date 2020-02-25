#ifndef READ_GRIDS_H
#define READ_GRIDS_H

#include <fftw3-mpi.h>

// N.B. Don't change these values!
enum grid_prop {
    DENSITY = 0,
    X_VELOCITY = 1,
    Y_VELOCITY = 2,
    Z_VELOCITY = 3
};

#ifdef __cplusplus
extern "C" {
#endif

double calc_resample_factor(int n_cell[3]);
void smooth_grid(double resample_factor, int n_cell[3], fftwf_complex* slab, ptrdiff_t slab_n_complex, ptrdiff_t slab_ix_start, ptrdiff_t slab_nix);
void subsample_grid(double resample_factor, int n_cell[3], int ix_hi_start, int nix_hi, float* slab_file, float* slab);
void read_grid(const enum grid_prop property, const int snapshot, float *slab);
int load_cached_slab(float* slab, int snapshot, const enum grid_prop property);
int cache_slab(float* slab, int snapshot, const enum grid_prop property);
int read_grid__gbptrees(const enum grid_prop property, const int snapshot, float* slab);
int read_grid__velociraptor(const enum grid_prop property, const int snapshot, float* slab);
void free_grids_cache(void);

#ifdef __cplusplus
}
#endif

#endif
