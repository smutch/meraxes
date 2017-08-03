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

__device__ void inline grid_index2indices(const int idx,const int dim,const int mode,int *ix,int *iy,int *iz);
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

__device__ int grid_index_gpu(int i, int j, int k, int dim, int type);
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
        grid_index2indices(i_real,grid_dim,INDEX_COMPLEX_HERM,&i_x,&i_y,&i_z);
        const int i_padded = grid_index_gpu(i_x,i_y,i_z,grid_dim,INDEX_PADDED);
        ((float *)grid)[i_padded] = fmaxf(((float *)grid)[i_padded], val);
    }
}

__global__
void filter_gpu(Complex *box,int local_ix_start,int slab_nx,int grid_dim,int n_complex,float R,double box_size_in,int filter_type);
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
        Complex *sfr_filtered_device);
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

// Presently, this is just a copy of what's in Meraxes
void _find_HII_bubbles_gpu(
    // input
    double redshift,
    MPI_Comm mpi_comm,
    int mpi_rank,
    double box_size,
    int ReionGridDim,
    int local_nix,
    int flag_ReionUVBFlag,
    double ReionEfficiency,
    double ReionNionPhotPerBary,
    double UnitLength_in_cm,
    double UnitMass_in_g,
    double UnitTime_in_s,
    double ReionRBubbleMax,
    double ReionRBubbleMin,
    double ReionDeltaRFactor,
    double ReionGammaHaloBias,
    double ReionAlphaUV,
    double ReionEscapeFrac,

    bool validation_output,

    // preallocated 1D grids (local_nix * ReionGridDim * ReionGridDim)
    float *J_21,  // real
    float *r_bubble, // real

    // input grids
    float *deltax,  // real & padded
    float *stars,  // real & padded
    float *sfr,  // real & padded

    // preallocated
    Complex *deltax_filtered_device_in,  // complex
    Complex *stars_filtered_device_in,  // complex
    Complex *sfr_filtered_device_in,  // complex

    // length = mpi.size
    ptrdiff_t *slabs_n_complex,
    ptrdiff_t *slabs_ix_start,

    // output - preallocated real grids (local_nix * ReionGridDim * ReionGridDim)
    float *xH, // real
    float *z_at_ionization,
    float *J_21_at_ionization,

    // output - single values
    double *volume_weighted_global_xH,
    double *mass_weighted_global_xH
    )
{
  const double pixel_volume         = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
  const double total_n_cells        = pow((double)ReionGridDim, 3);
  const double inv_total_n_cells    = 1.f/total_n_cells;
  const int    slab_n_real          = local_nix * ReionGridDim * ReionGridDim;
  const int    slab_n_complex       = (int)(slabs_n_complex[mpi_rank]);
  double       cell_length_factor   = L_FACTOR;

  

  if (validation_output)
  {
    // prepare output file
    char fname[STRLEN];
    sprintf(fname, "validation_input-core%03d-z%.2f.h5", mpi_rank, redshift);
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // write all of the input values
    H5LTset_attribute_double(file_id, "/", "redshift", &redshift, 1);
    H5LTset_attribute_int(file_id, "/", "mpi_rank", &mpi_rank, 1);
    H5LTset_attribute_double(file_id, "/", "box_size", &box_size, 1);
    H5LTset_attribute_int(file_id, "/", "ReionGridDim", &ReionGridDim, 1);
    H5LTset_attribute_int(file_id, "/", "local_nix", &local_nix, 1);
    H5LTset_attribute_int(file_id, "/", "flag_ReionUVBFlag", &flag_ReionUVBFlag, 1);
    H5LTset_attribute_double(file_id, "/", "ReionEfficiency", &ReionEfficiency, 1);
    H5LTset_attribute_double(file_id, "/", "ReionNionPhotPerBary", &ReionNionPhotPerBary, 1);
    H5LTset_attribute_double(file_id, "/", "UnitLength_in_cm", &UnitLength_in_cm, 1);
    H5LTset_attribute_double(file_id, "/", "UnitMass_in_g", &UnitMass_in_g, 1);
    H5LTset_attribute_double(file_id, "/", "UnitTime_in_s", &UnitTime_in_s, 1);
    H5LTset_attribute_double(file_id, "/", "ReionRBubbleMax", &ReionRBubbleMax, 1);
    H5LTset_attribute_double(file_id, "/", "ReionRBubbleMin", &ReionRBubbleMin, 1);
    H5LTset_attribute_double(file_id, "/", "ReionDeltaRFactor", &ReionDeltaRFactor, 1);
    H5LTset_attribute_double(file_id, "/", "ReionGammaHaloBias", &ReionGammaHaloBias, 1);
    H5LTset_attribute_double(file_id, "/", "ReionAlphaUV", &ReionAlphaUV, 1);
    H5LTset_attribute_double(file_id, "/", "ReionEscapeFrac", &ReionEscapeFrac, 1);

    H5LTmake_dataset_float(file_id, "deltax", 1, (hsize_t []){slab_n_complex*2}, deltax);
    H5LTmake_dataset_float(file_id, "stars", 1, (hsize_t []){slab_n_complex*2}, stars);
    H5LTmake_dataset_float(file_id, "sfr", 1, (hsize_t []){slab_n_complex*2}, sfr);
    H5LTmake_dataset_float(file_id, "z_at_ionization", 1, (hsize_t []){slab_n_real}, z_at_ionization);
    H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, (hsize_t []){slab_n_real}, J_21_at_ionization);

    H5Fclose(file_id);
  }

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0;

  // Initialize arrays on the device
  cufftComplex *deltax_unfiltered_device  = (cufftComplex *)deltax;
  cufftComplex *stars_unfiltered_device   = (cufftComplex *)stars;
  cufftComplex *sfr_unfiltered_device     = (cufftComplex *)sfr;
  cufftComplex *deltax_filtered_device    = NULL;
  cufftComplex *stars_filtered_device     = NULL;
  cufftComplex *sfr_filtered_device       = NULL;
  float        *xH_device                 = NULL;
  float        *J_21_device               = NULL;
  float        *r_bubble_device           = NULL;
  float        *z_at_ionization_device    = NULL;
  float        *J_21_at_ionization_device = NULL;
  cudaMalloc((void**)&deltax_unfiltered_device,           sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&stars_unfiltered_device,            sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&sfr_unfiltered_device,              sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&deltax_filtered_device,             sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&stars_filtered_device,              sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&sfr_filtered_device,                sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&xH_device,                          sizeof(float)*slab_n_real);
  if(flag_ReionUVBFlag)
     cudaMalloc((void**)&J_21_device,                     sizeof(float)*slab_n_real);
  cudaMalloc((void**)&r_bubble_device,                    sizeof(float)*slab_n_real);
  cudaMalloc((void**)&z_at_ionization_device,             sizeof(float)*slab_n_real);
  cudaMalloc((void**)&J_21_at_ionization_device,          sizeof(float)*slab_n_real);

  // Perform host -> device transfer
  cudaMemcpy(deltax_unfiltered_device, deltax,            sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(stars_unfiltered_device,  stars,             sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(sfr_unfiltered_device,    sfr,               sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(z_at_ionization_device,   z_at_ionization,   sizeof(float)*slab_n_real,     cudaMemcpyHostToDevice);
  cudaMemcpy(J_21_at_ionization_device,J_21_at_ionization,sizeof(float)*slab_n_real,     cudaMemcpyHostToDevice);

  // Forward fourier transform to obtain k-space fields

  // Initialize cuFFT
  cufftHandle plan;
  cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_R2C);
  cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL);

  // Perform FFTs
  if (cufftExecR2C(plan,(cufftReal *)deltax_unfiltered_device,deltax_unfiltered_device) != CUFFT_SUCCESS ) {
    fprintf(stderr, "Cuda error 1.\n");
    return ;
  }
  if (cufftExecR2C(plan,(cufftReal *)stars_unfiltered_device,stars_unfiltered_device) != CUFFT_SUCCESS ) {
    fprintf(stderr, "Cuda error 2.\n");
    return ;
  }
  if (cufftExecR2C(plan,(cufftReal *)sfr_unfiltered_device,sfr_unfiltered_device) != CUFFT_SUCCESS ) {
    fprintf(stderr, "Cuda error 3.\n");
    return ;
  }

  // Clean-up 
  cufftDestroy(plan);

  // Make sure that the device has synchronized
  if(cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error 4.\n");
    return;
  }

  if (validation_output)
  {
    // prepare output file
    char fname[STRLEN];
    sprintf(fname, "validation_output-core%03d-z%.2f.h5", mpi_rank, redshift);
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
    cudaMemcpy(array_temp,deltax_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    H5LTmake_dataset_float(group, "deltax", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
    cudaMemcpy(array_temp,stars_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    H5LTmake_dataset_float(group, "stars", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
    cudaMemcpy(array_temp,sfr_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    H5LTmake_dataset_float(group, "sfr", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
    free(array_temp);
    H5Gclose(group);
    H5Fclose(file_id);
  }

  // Initialize GPU block and thread count
  int threads      = 256;  
  int grid_complex = (slab_n_complex+(threads-1))/threads;
  int grid_real    = (slab_n_real   +(threads-1))/threads;

  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  complex_vector_times_scalar<<<grid_complex, threads>>>(deltax_unfiltered_device,inv_total_n_cells,slab_n_complex);
  complex_vector_times_scalar<<<grid_complex, threads>>>(stars_unfiltered_device, inv_total_n_cells,slab_n_complex);
  complex_vector_times_scalar<<<grid_complex, threads>>>(sfr_unfiltered_device,   inv_total_n_cells,slab_n_complex);

  if(cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error 4a.\n");
    return;
  }

  // Initialize a few of the output grids
  set_array_gpu<<<grid_real,threads>>>(xH_device,      slab_n_real,1.f);
  if(cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error 4b.\n");
    return;
  }
  set_array_gpu<<<grid_real,threads>>>(r_bubble_device,slab_n_real,0.f);
  if(cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error 4c.\n");
    return;
  }
  if (flag_ReionUVBFlag){
     set_array_gpu<<<grid_real,threads>>>(J_21_device,slab_n_real,0.f);
     if(cudaThreadSynchronize() != cudaSuccess){
       fprintf(stderr, "Cuda error 4d.\n");
       return;
     }
  }

  // Loop through filter radii
  double R                     = fmin(ReionRBubbleMax, L_FACTOR * box_size); // Mpc/h
  bool   flag_last_filter_step = false;
  int    i_R=0; // for debugging
  while(!flag_last_filter_step)
  {
    i_R++;
    // check to see if this is our last filtering step
    if( ((R / ReionDeltaRFactor) <= (cell_length_factor * box_size / (double)ReionGridDim))
        || ((R / ReionDeltaRFactor) <= ReionRBubbleMin) )
    {
      flag_last_filter_step = true;
      R                     = cell_length_factor * box_size / (double)ReionGridDim;
    }

    mlog(".", MLOG_CONT);

    // copy the k-space grids
    cudaMemcpy(deltax_filtered_device,deltax_unfiltered_device,sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    cudaMemcpy(stars_filtered_device, stars_unfiltered_device, sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    cudaMemcpy(sfr_filtered_device,   sfr_unfiltered_device,   sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);

    // Perform convolution
    int local_ix_start = 0; // make this = (int)(slabs_ix_start[mpi_rank]);
    if(!flag_last_filter_step){
       filter_gpu<<<grid_complex,threads>>>(deltax_filtered_device,local_ix_start,local_nix,ReionGridDim,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType);
       filter_gpu<<<grid_complex,threads>>>(stars_filtered_device, local_ix_start,local_nix,ReionGridDim,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType);
       filter_gpu<<<grid_complex,threads>>>(sfr_filtered_device,   local_ix_start,local_nix,ReionGridDim,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType);
    }

    if (validation_output && i_R==1)
    {
        // prepare output file
        char fname[STRLEN];
        sprintf(fname, "validation_test-core%03d-z%.2f.h5", mpi_rank, redshift);
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "deltax_filtered_device", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "stars_filtered_device", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "sfr_filtered_device", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Gclose(group);
        H5Fclose(file_id);
    }

    // inverse fourier transform back to real space

    // Initialize cuFFT
    cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_C2R);
    cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL);
    
    // Perform FFTs
    if (cufftExecC2R(plan,(cufftComplex *)deltax_filtered_device, (cufftReal *)deltax_filtered_device) != CUFFT_SUCCESS ) {
      fprintf(stderr, "Cuda error 101.\n");
      return ;
    }
    if (cufftExecC2R(plan,(cufftComplex *)stars_filtered_device, (cufftReal *)stars_filtered_device) != CUFFT_SUCCESS ) {
      fprintf(stderr, "Cuda error 102.\n");
      return ;
    }
    if (cufftExecC2R(plan,(cufftComplex *)sfr_filtered_device, (cufftReal *)sfr_filtered_device) != CUFFT_SUCCESS ) {
      fprintf(stderr, "Cuda error 103.\n");
      return ;
    }

    // Clean-up device
    cufftDestroy(plan);

    if (validation_output && i_R==1)
    {
        // prepare output file
        char fname[STRLEN];
        sprintf(fname, "validation_test-core%03d-z%.2f.h5", mpi_rank, redshift);
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax_filtered_device", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "stars_filtered_device", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "sfr_filtered_device", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Fclose(file_id);
    }

    // Perform sanity checks to account for aliasing effects
    sanity_check_aliasing<<<grid_real,threads>>>(deltax_filtered_device,ReionGridDim,slab_n_real,-1.f + REL_TOL);
    sanity_check_aliasing<<<grid_real,threads>>>(stars_filtered_device, ReionGridDim,slab_n_real,0.);
    sanity_check_aliasing<<<grid_real,threads>>>(sfr_filtered_device,   ReionGridDim,slab_n_real,0.);

    if (validation_output && i_R==1)
    {
        // prepare output file
        char fname[STRLEN];
        sprintf(fname, "validation_test-core%03d-z%.2f.h5", mpi_rank, redshift);
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax_filtered_device_checked", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "stars_filtered_device_checked", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "sfr_filtered_device_checked", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Fclose(file_id);
    }

    // Main loop through the box...
    const double J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI)
      * ReionAlphaUV * PLANCK
      * 1e21 * ReionEscapeFrac
      * R *UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS
      * UnitMass_in_g / pow(UnitLength_in_cm, 3) / UnitTime_in_s;
    const double inv_pixel_volume = 1.f/pixel_volume;

    find_HII_bubbles_gpu_main_loop<<<grid_real,threads>>>(
        redshift, //
        slab_n_real, //
        flag_last_filter_step, //
        flag_ReionUVBFlag, //
        ReionGridDim, //
        R, //
        RtoM(R), //
        ReionEfficiency, //
        inv_pixel_volume, //
        J_21_aux_constant, //
        ReionGammaHaloBias, //
        xH_device,
        J_21_device,
        r_bubble_device,
        J_21_at_ionization_device,
        z_at_ionization_device,
        deltax_filtered_device,
        stars_filtered_device,
        sfr_filtered_device);

    if (validation_output && i_R==1)
    {
        // prepare output file
        char fname[STRLEN];
        sprintf(fname, "validation_test-core%03d-z%.2f.h5", mpi_rank, redshift);
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*slab_n_real);
        cudaMemcpy(array_temp,xH_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "xH", 1, (hsize_t []){slab_n_real}, array_temp);
        cudaMemcpy(array_temp,r_bubble_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "r_bubble", 1, (hsize_t []){slab_n_real}, array_temp);
        if(flag_ReionUVBFlag){
            cudaMemcpy(array_temp,J_21_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
            H5LTmake_dataset_float(file_id, "J_21", 1, (hsize_t []){slab_n_real}, array_temp);
        }
        free(array_temp);
        H5Fclose(file_id);
    }

    R /= ReionDeltaRFactor;
  }

  // Perform device -> host transfer
  cudaMemcpy(z_at_ionization,   z_at_ionization_device,   sizeof(float) * 2 * slab_n_complex,cudaMemcpyDeviceToHost);
  cudaMemcpy(J_21_at_ionization,J_21_at_ionization_device,sizeof(float) * 2 * slab_n_complex,cudaMemcpyDeviceToHost);
  cudaMemcpy(xH,                xH_device,                sizeof(float) * slab_n_real,cudaMemcpyDeviceToHost);
  cudaMemcpy(deltax,            deltax_filtered_device,   sizeof(float) * slab_n_real,cudaMemcpyDeviceToHost);
  if(flag_ReionUVBFlag)
     cudaMemcpy(J_21,           J_21_device,              sizeof(float) * slab_n_real,cudaMemcpyDeviceToHost);

  // Clean-up device
  cudaFree(z_at_ionization_device);
  cudaFree(J_21_at_ionization_device);
  cudaFree(xH_device);
  if(flag_ReionUVBFlag)
     cudaFree(J_21_device);
  cudaFree(deltax_unfiltered_device);
  cudaFree(stars_unfiltered_device);
  cudaFree(sfr_unfiltered_device);
  cudaFree(deltax_filtered_device);
  cudaFree(stars_filtered_device);
  cudaFree(sfr_filtered_device);

  // Find the volume and mass weighted neutral fractions
  // TODO: The deltax grid will have rounding errors from forward and reverse
  //       FFT. Should cache deltax slabs prior to ffts and reuse here.
  *volume_weighted_global_xH = 0.0;
  *mass_weighted_global_xH   = 0.0;
  double mass_weight         = 0.0;

  int ix,iy,iz;
  for (ix = 0; ix < local_nix; ix++)
    for (iy = 0; iy < ReionGridDim; iy++)
      for (iz = 0; iz < ReionGridDim; iz++)
      {
        const int i_real   = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
        const int i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
        double density_over_mean    = 1.0 + (double)((float *)deltax)[i_padded];
        *volume_weighted_global_xH += (double)xH[i_real];
        *mass_weighted_global_xH   += (double)(xH[i_real]) * density_over_mean;
        mass_weight                += density_over_mean;
      }

  *volume_weighted_global_xH                        *= inv_total_n_cells;
  *mass_weighted_global_xH                          /= mass_weight;

  if (validation_output)
  {
    // prepare output file
    char fname[STRLEN];
    sprintf(fname, "validation_output-core%03d-z%.2f.h5", mpi_rank, redshift);
    hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

    H5LTmake_dataset_float(file_id, "xH", 1, (hsize_t []){slab_n_real}, xH);
    H5LTmake_dataset_float(file_id, "z_at_ionization", 1, (hsize_t []){slab_n_real}, z_at_ionization);
    H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, (hsize_t []){slab_n_real}, J_21_at_ionization);

    H5LTset_attribute_double(file_id, "/", "volume_weighted_global_xH", volume_weighted_global_xH, 1);
    H5LTset_attribute_double(file_id, "/", "mass_weighted_global_xH", mass_weighted_global_xH, 1);

    H5Fclose(file_id);
  }
}

