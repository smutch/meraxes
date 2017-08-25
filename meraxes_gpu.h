#ifndef _MERAXES_GPU_H
#define _MERAXES_GPU_H
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "utils.h"
#include "meraxes.h"

#ifdef __NVCC__
typedef float2 Complex;
#else
typedef fftwf_complex Complex;
#endif

#ifdef __NVCC__
__device__ void  inline grid_index2indices(const int idx,const int dim,const int local_ix_start,const int mode,int *ix,int *iy,int *iz);
__device__ int   inline grid_index_gpu(int i, int j, int k, int dim, int type);
__device__ float inline k_mag_gpu(const int n_x,const int n_y,const int n_z,const int dim,const float box_size);
__global__ void  set_array_gpu(float *array,int n_real,float val);
__global__ void  sanity_check_aliasing(Complex *grid,int grid_dim,int local_ix_start,int n_real,float val);
__global__ void  complex_vector_times_scalar(Complex *vector,double scalar,int n_complex);
__global__ void  filter_gpu(Complex *box,int slab_nx,int grid_dim,int local_ix_start,int n_complex,float R,double box_size_in,int filter_type);
__global__ void  find_HII_bubbles_gpu_main_loop(
                   float    redshift,
                   int      n_real,
                   int      flag_last_filter_step,
                   int      flag_ReionUVBFlag,
                   int      ReionGridDim,
                   int      local_ix_start,
                   float    R,
                   float    M,
                   float    ReionEfficiency,
                   float    inv_pixel_volume,
                   float    J_21_aux_constant,
                   double   ReionGammaHaloBias,
                   float   *xH,
                   float   *J_21,
                   float   *r_bubble,
                   float   *J_21_at_ionization,
                   float   *z_at_ionization,
                   Complex *deltax_filtered_device,
                   Complex *stars_filtered_device,
                   Complex *sfr_filtered_device);
#ifdef __cplusplus
extern "C" {
#endif
void  _find_HII_bubbles_gpu(double redshift,const bool flag_write_validation_output);
#ifdef __cplusplus
}
#endif
#endif // ifdef __NVCC__

#ifdef __cplusplus
extern "C" {
#endif
void find_HII_bubbles_driver(
    int        snapshot,
    const char *reference_directory,
    const bool flag_write_validation_data,
    timer_info *timer);
#ifdef __cplusplus
}
#endif
#define _MERAXES_GPU_H
#endif
