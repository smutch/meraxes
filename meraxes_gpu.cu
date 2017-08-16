#include "meraxes.h"
#include "meraxes_gpu.h"
#include "utils.h"
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

__device__ void inline grid_index2indices(const int idx,const int dim,const int mode,int *ix,int *iy,int *iz){
    int n;
    int remainder=idx;
    switch(mode){
        case INDEX_PADDED:
            n=(2*(dim/2+1));
            break;
        case INDEX_REAL:
            n=dim;
            break;
        case INDEX_COMPLEX_HERM:
            n=(dim/2+1);
            break;
        default:
            // THROW ERROR HERE
            break;
    }
    (*iz)=remainder%n;
    remainder=(remainder-(*iz))/n;
    (*iy)=remainder%dim;
    remainder=(remainder-(*iy))/dim;
    (*ix)=remainder;
}

__device__ int grid_index_gpu(int i, int j, int k, int dim, int type)
{
  int ind;

  switch(type)
  {
    case INDEX_PADDED:
      ind = k + (2 * (dim / 2 + 1)) * (j + dim * i);
      break;
    case INDEX_REAL:
      ind = k + dim * (j + dim * i);
      break;
    case INDEX_COMPLEX_HERM:
      ind = k + (dim / 2 + 1) * (j + dim * i);
      break;
  }

  return ind;
}

__global__
void set_array_gpu(float *array,int n_real,float val){
    int i_real = blockIdx.x*blockDim.x + threadIdx.x;
    if (i_real < n_real) array[i_real]=val;
}

__global__
void complex_vector_times_scalar(Complex *vector,double scalar,int n_complex){
    int i_complex = blockIdx.x*blockDim.x + threadIdx.x;
    if (i_complex < n_complex){
        vector[i_complex].x*=scalar;
        vector[i_complex].y*=scalar;
    }
}

__global__
void sanity_check_aliasing(Complex *grid,int grid_dim,int n_real,float val){
    int i_real = blockIdx.x*blockDim.x + threadIdx.x;
    if (i_real < n_real){
        int i_x,i_y,i_z;
        grid_index2indices(i_real,grid_dim,INDEX_REAL,&i_x,&i_y,&i_z);
        const int i_padded = grid_index_gpu(i_x,i_y,i_z,grid_dim,INDEX_PADDED);
        ((float *)grid)[i_padded] = fmaxf(((float *)grid)[i_padded], val);
    }
}

__global__
void filter_gpu(Complex *box,int local_ix_start,int slab_nx,int grid_dim,int n_complex,float R,double box_size_in,int filter_type){
    int i_complex = blockIdx.x*blockDim.x + threadIdx.x;
    if (i_complex < n_complex){

        // Map i_complex to grid indices
        int n_x,n_y,n_z;
        grid_index2indices(i_complex,grid_dim,INDEX_COMPLEX_HERM,&n_x,&n_y,&n_z);

        // Compute k_mag for given grid indices
        float box_size = box_size_in;
        float delta_k  = 2.0 * M_PI / box_size;
        float k_x      = n_x*delta_k;
        float k_y      = n_y*delta_k;
        float k_z      = n_z*delta_k;
        float k_mag    = sqrtf(k_x * k_x + k_y * k_y + k_z * k_z);

        // Calculate filter
        const int   i_complex_herm = grid_index_gpu(n_x,n_y,n_z,grid_dim,INDEX_COMPLEX_HERM);
        float       kR             = k_mag * R;
        float       scalar         = 0.;
        int         support        = false;
        switch(filter_type)
        {
          case 0:   // Real space top-hat
            scalar  = (3.0 * (sinf(kR) / powf(kR, 3) - cosf(kR) / powf(kR, 2)));
            support = (kR>1e-4);
            break;

          case 1:                  // k-space top hat
            kR     *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
            scalar  = 0.f;
            support = (kR>1);
            break;

          case 2:            // Gaussian
            kR     *= 0.643; // Equates integrated volume to the real space top-hat
            scalar  = powf(M_E, -kR * kR / 2.0);
            support = true;
            break;
          // Implement this check before the kernel!!!!!!!!!!
          //default:
          //  if (i==0)
          //  {
          //    mlog_error("ReionFilterType.c: Warning, ReionFilterType type %d is undefined!", filter_type);
          //    ABORT(EXIT_FAILURE);
          //  }
          //  break;
        }

        // Apply filter
        if(support){
            box[i_complex_herm].x*=scalar;
            box[i_complex_herm].y*=scalar;
        }
    }
}

__global__
void find_HII_bubbles_gpu_main_loop(
        float    redshift,
        int      n_real,
        int      flag_last_filter_step,
        int      flag_ReionUVBFlag,
        int      ReionGridDim,
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
        Complex *sfr_filtered_device){
    int i_real = blockIdx.x*blockDim.x + threadIdx.x;
    if(i_real < n_real){
        int ix,iy,iz;
        grid_index2indices(i_real,ReionGridDim,INDEX_REAL,&ix,&iy,&iz); // TODO: need to pass local_i_start to this eventually
        const int i_padded = grid_index_gpu(ix,iy,iz, ReionGridDim, INDEX_PADDED);

        double density_over_mean = 1.0 + (double)((float *)deltax_filtered_device)[i_padded];

        double f_coll_stars      =  (double)((float *)stars_filtered_device)[i_padded] / (M * density_over_mean)
                             * (4.0 / 3.0) * M_PI * (R*R*R)  * inv_pixel_volume;

        double sfr_density       = (double)((float *)sfr_filtered_device)[i_padded] * inv_pixel_volume; // In internal units

        float J_21_aux;
        if (flag_ReionUVBFlag)
          J_21_aux = (float)(sfr_density * J_21_aux_constant);

        // Check if ionised!
        if (f_coll_stars > 1.0 / ReionEfficiency)   // IONISED!!!!
        {
          // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
          if (xH[i_real] > REL_TOL)
            if(flag_ReionUVBFlag)
              J_21[i_real] = J_21_aux;

          // Mark as ionised
          xH[i_real]       = 0;

          // Record radius
          r_bubble[i_real] = (float)R;
        }
        // Check if this is the last filtering step.
        // If so, assign partial ionisations to those cells which aren't fully ionised
        else if (flag_last_filter_step && (xH[i_real] > REL_TOL))
          xH[i_real] = (float)(1.0 - f_coll_stars * ReionEfficiency);

        // Check if new ionisation
        float *z_in = z_at_ionization;
        if ( (xH[i_real] < REL_TOL) && (z_in[i_real] < 0) )   // New ionisation!
        {
          z_in[i_real] = (float)redshift;
          if (flag_ReionUVBFlag)
            J_21_at_ionization[i_real] = J_21_aux * (float)ReionGammaHaloBias;
        }
    }
}

