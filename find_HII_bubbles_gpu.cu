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

__global__
void complex_vector_times_scalar(Complex *vector,double scalar,int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n){
        vector[i].x*=scalar;
        vector[i].y*=scalar;
    }
}

__global__
void sanity_check_aliasing(Complex *grid,int n,float val){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n)
          ((float *)grid)[i] = fmaxf(((float *)grid)[i], val);
}

__device__ void inline index2indices_FFT_k(const int dim,int index,int *i_k);
__device__ void inline index2indices_FFT_k(const int dim,int index,int *i_k){ // should match mode=INDEX_COMPLEX_HERM in grid_index ... check this!
  int i_d,j_d;
  int remainder;
  for(j_d=2,remainder=index;j_d>=0;j_d--){
    i_d=j_d;
    i_k[i_d]  =remainder%dim;
    remainder-=i_k[i_d];
    remainder/=dim;
  }
}

__device__ float k_mag_of_index(const int dim,int index);
__device__ float k_mag_of_index(const int dim,int index){
    int idxs[3];
    index2indices_FFT_k(dim,index,idxs);
    float k_mag = 0.f;
    return(k_mag);
}

__global__
void filter_gpu(Complex *grid,int dim,int n,float R){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n){
        float kR = R*k_mag_of_index(dim,i);
        float scalar =0.;
        int   support=false;
int filter_type = 0;
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

          case 2:        // Gaussian
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
        if(support){
            grid[i].x*=scalar;
            grid[i].y*=scalar;
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
    Complex *deltax_filtered_in,  // complex
    Complex *stars_filtered_in,  // complex
    Complex *sfr_filtered_in,  // complex

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

    H5Fclose(file_id);
  }

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0;

  // Init J_21
  int ii;
  if (flag_ReionUVBFlag)
    for(ii = 0; ii < slab_n_real; ii++)
      J_21[ii] = 0.0;

  // Init xH
  for(ii = 0; ii < slab_n_real; ii++)
    xH[ii] = 1.0;

  // Init r_bubble
  for(ii = 0; ii < slab_n_real; ii++)
    r_bubble[ii] = 0.0;

  // Forward fourier transform to obtain k-space fields
  // Initialize cuFFT
  cufftHandle plan;
  cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_R2C);
  cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL);

  // Initialize arrays on the device
  cufftComplex *deltax_unfiltered = (cufftComplex *)deltax;
  cufftComplex *stars_unfiltered  = (cufftComplex *)stars;
  cufftComplex *sfr_unfiltered    = (cufftComplex *)sfr;
  cufftComplex *deltax_filtered   = NULL;
  cufftComplex *stars_filtered    = NULL;
  cufftComplex *sfr_filtered      = NULL;
  cudaMalloc((void**)&deltax_unfiltered,sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&stars_unfiltered, sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&sfr_unfiltered,   sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&deltax_filtered,  sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&stars_filtered,   sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&sfr_filtered,     sizeof(cufftComplex)*slab_n_complex);
  cudaMemcpy(deltax,deltax_unfiltered,  sizeof(cufftComplex)*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(stars, stars_unfiltered,   sizeof(cufftComplex)*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(sfr,   sfr_unfiltered,     sizeof(cufftComplex)*slab_n_complex,cudaMemcpyHostToDevice);

  // Perform FFTs
  if (cufftExecR2C(plan,(cufftReal *)deltax_unfiltered,deltax_unfiltered) != CUFFT_SUCCESS ) {
    fprintf(stderr, "Cuda error 1.\n");
    return ;
  }
  if (cufftExecR2C(plan,(cufftReal *)stars_unfiltered,stars_unfiltered) != CUFFT_SUCCESS ) {
    fprintf(stderr, "Cuda error 2.\n");
    return ;
  }
  if (cufftExecR2C(plan,(cufftReal *)sfr_unfiltered,sfr_unfiltered) != CUFFT_SUCCESS ) {
    fprintf(stderr, "Cuda error 3.\n");
    return ;
  }

  // Make sure that the device has synchronized
  if(cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error 4.\n");
    return;
  }

  // Clean-up the device
  cufftDestroy(plan);

  if (validation_output)
  {
    // prepare output file
    char fname[STRLEN];
    sprintf(fname, "validation_output-core%03d-z%.2f.h5", mpi_rank, redshift);
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5LTmake_dataset_float(group, "deltax", 1, (hsize_t []){slab_n_complex * 2}, deltax);
    H5LTmake_dataset_float(group, "stars", 1, (hsize_t []){slab_n_complex * 2}, stars);
    H5LTmake_dataset_float(group, "sfr", 1, (hsize_t []){slab_n_complex * 2}, sfr);

    H5Gclose(group);
    H5Fclose(file_id);
  }

  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  int threads = 256;  
  int grid    = (slab_n_complex+255)/256;
  complex_vector_times_scalar<<<grid, threads>>>(deltax_unfiltered,inv_total_n_cells,slab_n_complex);
  complex_vector_times_scalar<<<grid, threads>>>(stars_unfiltered, inv_total_n_cells,slab_n_complex);
  complex_vector_times_scalar<<<grid, threads>>>(sfr_unfiltered,   inv_total_n_cells,slab_n_complex);

  // Loop through filter radii
  double R                     = fmin(ReionRBubbleMax, L_FACTOR * box_size); // Mpc/h

  bool  flag_last_filter_step = false;

  while(!flag_last_filter_step)
  {
    // check to see if this is our last filtering step
    if( ((R / ReionDeltaRFactor) <= (cell_length_factor * box_size / (double)ReionGridDim))
        || ((R / ReionDeltaRFactor) <= ReionRBubbleMin) )
    {
      flag_last_filter_step = true;
      R                     = cell_length_factor * box_size / (double)ReionGridDim;
    }

    mlog(".", MLOG_CONT);

    // copy the k-space grids
    cudaMemcpy(deltax_filtered,deltax_unfiltered,sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    cudaMemcpy(stars_filtered, stars_unfiltered, sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    cudaMemcpy(sfr_filtered,   sfr_unfiltered,   sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    if(!flag_last_filter_step){
       filter_gpu<<<grid,threads>>>(deltax_filtered,ReionGridDim,slab_n_complex,(float)R);
       filter_gpu<<<grid,threads>>>(stars_filtered, ReionGridDim,slab_n_complex,(float)R);
       filter_gpu<<<grid,threads>>>(sfr_filtered,   ReionGridDim,slab_n_complex,(float)R);
    }

    // inverse fourier transform back to real space
    // Initialize cuFFT
    cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_C2R);
    cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL);
    
    // Perform FFTs
    if (cufftExecC2R(plan,(cufftComplex *)deltax_filtered, (cufftReal *)deltax_filtered) != CUFFT_SUCCESS ) {
      fprintf(stderr, "Cuda error 101.\n");
      return ;
    }
    if (cufftExecC2R(plan,(cufftComplex *)stars_filtered, (cufftReal *)stars_filtered) != CUFFT_SUCCESS ) {
      fprintf(stderr, "Cuda error 102.\n");
      return ;
    }
    if (cufftExecC2R(plan,(cufftComplex *)sfr_filtered, (cufftReal *)sfr_filtered) != CUFFT_SUCCESS ) {
      fprintf(stderr, "Cuda error 103.\n");
      return ;
    }

    // Clean-up device
    cufftDestroy(plan);

    // Perform sanity checks to account for aliasing effects
    sanity_check_aliasing<<<grid,threads>>>(deltax_filtered,slab_n_complex,-1.f + REL_TOL);
    sanity_check_aliasing<<<grid,threads>>>(stars_filtered, slab_n_complex,0.);
    sanity_check_aliasing<<<grid,threads>>>(sfr_filtered,   slab_n_complex,0.);

    /*
     * Main loop through the box...
     */

    const double J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI)
      * ReionAlphaUV * PLANCK
      * 1e21 * ReionEscapeFrac
      * R *UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS
      * UnitMass_in_g / pow(UnitLength_in_cm, 3) / UnitTime_in_s;

#ifdef __NVCC__
#else
    const double inv_pixel_volume = 1.f/pixel_volume;
    for (ix = 0; ix < local_nix; ix++)
      for (iy = 0; iy < ReionGridDim; iy++)
        for (iz = 0; iz < ReionGridDim; iz++)
        {
          const int i_real   = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
          const int i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

          double density_over_mean = 1.0 + (double)((float *)deltax_filtered)[i_padded];

          double f_coll_stars      =  (double)((float *)stars_filtered)[i_padded] / (RtoM(R) * density_over_mean)
                               * (4.0 / 3.0) * M_PI * pow(R,3.0)  * inv_pixel_volume;

          double sfr_density       = (double)((float *)sfr_filtered)[i_padded] * inv_pixel_volume; // In internal units

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
    // iz
#endif

    R /= ReionDeltaRFactor;
  }


  // Find the volume and mass weighted neutral fractions
  // TODO: The deltax grid will have rounding errors from forward and reverse
  //       FFT. Should cache deltax slabs prior to ffts and reuse here.
  *volume_weighted_global_xH = 0.0;
  *mass_weighted_global_xH   = 0.0;
  double mass_weight         = 0.0;

#ifdef __NVCC__
#else
  for (ix = 0; ix < local_nix; ix++)
    for (iy = 0; iy < ReionGridDim; iy++)
      for (iz = 0; iz < ReionGridDim; iz++)
      {
        const int i_real   = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
        const int i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
        double density_over_mean    = 1.0 + (double)((float *)deltax_filtered)[i_padded];
        *volume_weighted_global_xH += (double)xH[i_real];
        *mass_weighted_global_xH   += (double)(xH[i_real]) * density_over_mean;
        mass_weight                += density_over_mean;
      }
#endif

  MPI_Allreduce(MPI_IN_PLACE, &volume_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &mass_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &mass_weight, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

  *volume_weighted_global_xH                        *= inv_total_n_cells;
  *mass_weighted_global_xH                          /= mass_weight;

  // Clean-up
  cudaFree(deltax_unfiltered);
  cudaFree(stars_unfiltered);
  cudaFree(sfr_unfiltered);
  cudaFree(deltax_filtered);
  cudaFree(stars_filtered);
  cudaFree(sfr_filtered);

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

void find_HII_bubbles_driver(
    double redshift,
    void  (*_find_HII_bubbles_passed)(
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
        Complex *deltax_filtered_in,  // complex
        Complex *stars_filtered_in,  // complex
        Complex *sfr_filtered_in,  // complex
    
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
        ),
    const char *reference_directory,
    timer_info *timer)
{
    int mpi_rank=run_globals.mpi_rank;
    int n_rank  =run_globals.mpi_size;
    if(n_rank!=1){
        if(mpi_rank==0) fprintf(stderr,"n_rank=%d but only n_rank==1 supported at this point.\n",n_rank);
        exit(1);
    }
    
    // Open inputs file
    char fname[STRLEN];
    sprintf(fname, "%s/validation_input-core%03d-z%.2f.h5",reference_directory,mpi_rank, redshift);
    hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    // Read input attributes
    double box_size;
    int    ReionGridDim;
    int    local_nix;
    int    flag_ReionUVBFlag;
    double ReionEfficiency;
    double ReionNionPhotPerBary;
    double UnitLength_in_cm;
    double UnitMass_in_g;
    double UnitTime_in_s;
    double ReionRBubbleMax;
    double ReionRBubbleMin;
    double ReionDeltaRFactor;
    double ReionGammaHaloBias;
    double ReionAlphaUV;
    double ReionEscapeFrac;
    H5LTget_attribute_double(file_id, "/", "redshift", &redshift);
    H5LTget_attribute_int   (file_id, "/", "mpi_rank", &mpi_rank);
    H5LTget_attribute_double(file_id, "/", "box_size", &box_size);
    H5LTget_attribute_int   (file_id, "/", "ReionGridDim", &ReionGridDim);
    H5LTget_attribute_int   (file_id, "/", "local_nix", &local_nix);
    H5LTget_attribute_int   (file_id, "/", "flag_ReionUVBFlag", &flag_ReionUVBFlag);
    H5LTget_attribute_double(file_id, "/", "ReionEfficiency", &ReionEfficiency);
    H5LTget_attribute_double(file_id, "/", "ReionNionPhotPerBary", &ReionNionPhotPerBary);
    H5LTget_attribute_double(file_id, "/", "UnitLength_in_cm", &UnitLength_in_cm);
    H5LTget_attribute_double(file_id, "/", "UnitMass_in_g", &UnitMass_in_g);
    H5LTget_attribute_double(file_id, "/", "UnitTime_in_s", &UnitTime_in_s);
    H5LTget_attribute_double(file_id, "/", "ReionRBubbleMax", &ReionRBubbleMax);
    H5LTget_attribute_double(file_id, "/", "ReionRBubbleMin", &ReionRBubbleMin);
    H5LTget_attribute_double(file_id, "/", "ReionDeltaRFactor", &ReionDeltaRFactor);
    H5LTget_attribute_double(file_id, "/", "ReionGammaHaloBias", &ReionGammaHaloBias);
    H5LTget_attribute_double(file_id, "/", "ReionAlphaUV", &ReionAlphaUV);
    H5LTget_attribute_double(file_id, "/", "ReionEscapeFrac", &ReionEscapeFrac);

    // Initialize fftw here (we need to do this because we need 'slab_n_complex' to set the size of the input datasets
    ptrdiff_t local_nix_2;
    ptrdiff_t local_ix_start;
    int       local_n_complex = fftwf_mpi_local_size_3d(ReionGridDim, ReionGridDim, ReionGridDim / 2 + 1, run_globals.mpi_comm, &local_nix_2, &local_ix_start);
    if(local_nix!=local_nix_2){
        printf("Error: local_nix!=local_nix_2 (ie. %d!=%d)\n",local_nix,local_nix_2);
        exit(1);
    }
    int       slab_n_complex    =  local_n_complex;
    int       slab_n_real       =  local_nix*ReionGridDim*ReionGridDim;
    ptrdiff_t slabs_n_complex[] = {slab_n_complex}; // works for 1 core only
    ptrdiff_t slabs_ix_start[]  = {0};              // works for 1 core only

    // Read input datasets
    float *deltax = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    float *stars  = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    float *sfr    = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    H5LTread_dataset_float(file_id, "deltax", deltax);
    H5LTread_dataset_float(file_id, "sfr",  sfr);
    H5LTread_dataset_float(file_id, "stars", stars);
    H5Fclose(file_id);

    // Initialize outputs    
    float   *J_21                     =(float   *)malloc(sizeof(float)  *slab_n_real);
    float   *r_bubble                 =(float   *)malloc(sizeof(float)  *slab_n_real); 
    Complex *deltax_filtered          =(Complex *)malloc(sizeof(Complex)*slab_n_complex);
    Complex *stars_filtered           =(Complex *)malloc(sizeof(Complex)*slab_n_complex);  
    Complex *sfr_filtered             =(Complex *)malloc(sizeof(float)  *slab_n_complex);
    float   *xH                       =(float   *)malloc(sizeof(float)  *slab_n_real);
    float   *z_at_ionization          =(float   *)calloc(slab_n_real,sizeof(float));
    float   *J_21_at_ionization       =(float   *)calloc(slab_n_real,sizeof(float));
    double   volume_weighted_global_xH;
    double   mass_weighted_global_xH;
    
    // Call the version of find_HII_bubbles we've been passed
    timer_start(timer);
    _find_HII_bubbles_passed(
        redshift, //
        run_globals.mpi_comm,
        mpi_rank,
        box_size, //
        ReionGridDim, //
        local_nix, //
        flag_ReionUVBFlag, //
        ReionEfficiency, //
        ReionNionPhotPerBary, //
        UnitLength_in_cm, //
        UnitMass_in_g, //
        UnitTime_in_s, //
        ReionRBubbleMax, //
        ReionRBubbleMin, //
        ReionDeltaRFactor, //
        ReionGammaHaloBias, //
        ReionAlphaUV, //
        ReionEscapeFrac, //
        true, //
        J_21,  
        r_bubble, 
        deltax, //
        stars, //
        sfr, //
        deltax_filtered,  
        stars_filtered,  
        sfr_filtered,  
        slabs_n_complex,
        slabs_ix_start,
        xH, 
        z_at_ionization,
        J_21_at_ionization,
        &volume_weighted_global_xH,
        &mass_weighted_global_xH);
    timer_stop(timer);
    
    // Clean-up
    free(J_21);
    free(r_bubble);
    free(deltax);
    free(stars);
    free(sfr);
    free(deltax_filtered);
    free(stars_filtered);
    free(sfr_filtered);
    free(xH);
    free(z_at_ionization);
    free(J_21_at_ionization);
}

