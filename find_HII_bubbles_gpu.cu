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
    float *J_21,     // real
    float *r_bubble, // real

    // input grids
    float *deltax,  // real & padded
    float *stars,   // real & padded
    float *sfr,     // real & padded

    // preallocated
    Complex *deltax_filtered,  // complex
    Complex *stars_filtered,   // complex
    Complex *sfr_filtered,     // complex

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

  if (validation_output)
  {
    // prepare output file
    char fname[STRLEN];
    sprintf(fname, "validation_input-core%03d-z%.2f.h5", mpi_rank, redshift);
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // write all of the input values
    H5LTset_attribute_double(file_id, "/", "redshift",             &redshift, 1);
    H5LTset_attribute_int   (file_id, "/", "mpi_rank",             &mpi_rank, 1);
    H5LTset_attribute_double(file_id, "/", "box_size",             &box_size, 1);
    H5LTset_attribute_int   (file_id, "/", "ReionGridDim",         &ReionGridDim, 1);
    H5LTset_attribute_int   (file_id, "/", "local_nix",            &local_nix, 1);
    H5LTset_attribute_int   (file_id, "/", "flag_ReionUVBFlag",    &flag_ReionUVBFlag, 1);
    H5LTset_attribute_double(file_id, "/", "ReionEfficiency",      &ReionEfficiency, 1);
    H5LTset_attribute_double(file_id, "/", "ReionNionPhotPerBary", &ReionNionPhotPerBary, 1);
    H5LTset_attribute_double(file_id, "/", "UnitLength_in_cm",     &UnitLength_in_cm, 1);
    H5LTset_attribute_double(file_id, "/", "UnitMass_in_g",        &UnitMass_in_g, 1);
    H5LTset_attribute_double(file_id, "/", "UnitTime_in_s",        &UnitTime_in_s, 1);
    H5LTset_attribute_double(file_id, "/", "ReionRBubbleMax",      &ReionRBubbleMax, 1);
    H5LTset_attribute_double(file_id, "/", "ReionRBubbleMin",      &ReionRBubbleMin, 1);
    H5LTset_attribute_double(file_id, "/", "ReionDeltaRFactor",    &ReionDeltaRFactor, 1);
    H5LTset_attribute_double(file_id, "/", "ReionGammaHaloBias",   &ReionGammaHaloBias, 1);
    H5LTset_attribute_double(file_id, "/", "ReionAlphaUV",         &ReionAlphaUV, 1);
    H5LTset_attribute_double(file_id, "/", "ReionEscapeFrac",      &ReionEscapeFrac, 1);

    H5LTmake_dataset_float(file_id, "deltax",             1, (hsize_t []){slab_n_complex*2}, deltax);
    H5LTmake_dataset_float(file_id, "stars",              1, (hsize_t []){slab_n_complex*2}, stars);
    H5LTmake_dataset_float(file_id, "sfr",                1, (hsize_t []){slab_n_complex*2}, sfr);
    H5LTmake_dataset_float(file_id, "z_at_ionization",    1, (hsize_t []){slab_n_real},      z_at_ionization);
    H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, (hsize_t []){slab_n_real},      J_21_at_ionization);

    H5Fclose(file_id);
  }

  // Initialize device arrays
  cufftComplex *deltax_unfiltered_device  = NULL;
  cufftComplex *stars_unfiltered_device   = NULL;
  cufftComplex *sfr_unfiltered_device     = NULL;
  cufftComplex *deltax_filtered_device    = NULL;
  cufftComplex *stars_filtered_device     = NULL;
  cufftComplex *sfr_filtered_device       = NULL;
  float        *xH_device                 = NULL;
  float        *r_bubble_device           = NULL;
  float        *z_at_ionization_device    = NULL;
  float        *J_21_at_ionization_device = NULL;
  float        *J_21_device               = NULL;
  cudaMalloc((void**)&deltax_unfiltered_device, sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&stars_unfiltered_device,  sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&sfr_unfiltered_device,    sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&deltax_filtered_device,   sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&stars_filtered_device,    sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&sfr_filtered_device,      sizeof(cufftComplex)*slab_n_complex);
  cudaMalloc((void**)&xH_device,                sizeof(float)*slab_n_real);
  cudaMalloc((void**)&r_bubble_device,          sizeof(float)*slab_n_real);
  cudaMalloc((void**)&z_at_ionization_device,   sizeof(float)*slab_n_real);
  cudaMalloc((void**)&J_21_at_ionization_device,sizeof(float)*slab_n_real);
  if(flag_ReionUVBFlag)
     cudaMalloc((void**)&J_21_device,           sizeof(float)*slab_n_real);

#ifdef GPU_HYBRID
  // Forward fourier transform to obtain k-space fields
  fftwf_complex *deltax_unfiltered = (fftwf_complex *)deltax;  // WATCH OUT!
  fftwf_plan     plan              = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax, deltax_unfiltered, mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  fftwf_complex *stars_unfiltered = (fftwf_complex *)stars;  // WATCH OUT!
  plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars, stars_unfiltered, mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  fftwf_complex *sfr_unfiltered = (fftwf_complex *)sfr;  // WATCH OUT!
  plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr, sfr_unfiltered, mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
#endif

  // Perform host -> device transfer
  cudaMemcpy(deltax_unfiltered_device, deltax,            sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(stars_unfiltered_device,  stars,             sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(sfr_unfiltered_device,    sfr,               sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
  cudaMemcpy(z_at_ionization_device,   z_at_ionization,   sizeof(float)*slab_n_real,     cudaMemcpyHostToDevice);
  cudaMemcpy(J_21_at_ionization_device,J_21_at_ionization,sizeof(float)*slab_n_real,     cudaMemcpyHostToDevice);

#ifndef GPU_HYBRID
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
#endif

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
  if(false){
    float *array_temp = (float *)malloc(sizeof(float)*(slab_n_complex*2));
    char fname[STRLEN];
    sprintf(fname, "%s/validation_output-core%03d-z%.2f.h5", "../test/", mpi_rank, redshift);
    hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group   = H5Gopen(file_id, "kspace",H5P_DEFAULT);
    H5LTread_dataset_float(group,"/kspace/deltax",array_temp); cudaMemcpy(deltax_unfiltered_device,array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
    H5LTread_dataset_float(group,"/kspace/stars", array_temp); cudaMemcpy(stars_unfiltered_device, array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
    H5LTread_dataset_float(group,"/kspace/sfr",   array_temp); cudaMemcpy(sfr_unfiltered_device,   array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
    H5Fclose(file_id);
    free(array_temp);
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

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  double cell_length_factor = L_FACTOR;
  if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0;

  // Initialize inverse FFTs
#ifdef GPU_HYBRID
  fftwf_plan plan_deltax = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, mpi_comm, FFTW_ESTIMATE);
  fftwf_plan plan_stars  = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)stars_filtered,  (float *)stars_filtered,  mpi_comm, FFTW_ESTIMATE);
  fftwf_plan plan_sfr    = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)sfr_filtered,    (float *)sfr_filtered,    mpi_comm, FFTW_ESTIMATE);
#else
  cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_C2R);
  cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL);
#endif

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

    char fname[STRLEN];
    char fname_full_dump[STRLEN];
    char fname_ref[STRLEN];
    sprintf(fname, "validation_test-core%03d-z%.2f_%03d.h5", mpi_rank, redshift,i_R);
    if(redshift>10.){
        sprintf(fname_full_dump,"validation_test-core%03d-z%.2f_%03d.h5", mpi_rank,10.11,i_R);
        sprintf(fname_ref,      "%s/validation_test-core%03d-z%.2f_%03d.h5", "../test/",mpi_rank,10.11,i_R);
    }
    else if(redshift>6.){
        sprintf(fname_full_dump,"validation_test-core%03d-z%.2f_%03d.h5", mpi_rank, 9.03,i_R);
        sprintf(fname_ref,      "%s/validation_test-core%03d-z%.2f_%03d.h5", "../test/",mpi_rank, 9.03,i_R);
    }
    else{
        sprintf(fname_full_dump,"validation_test-core%03d-z%.2f_%03d.h5", mpi_rank, 5.95,i_R);
        sprintf(fname_ref,      "%s/validation_test-core%03d-z%.2f_%03d.h5", "../test/",mpi_rank, 5.95,i_R);
    }
    if (validation_output && i_R==1 || !strcmp(fname,fname_full_dump))
    {
        // prepare output file
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "deltax_unfiltered", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "stars_unfiltered", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "sfr_unfiltered", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Gclose(group);
        H5Fclose(file_id);
    }
    if(false){ 
      float *array_temp = (float *)malloc(sizeof(float)*(slab_n_complex*2));
      hid_t file_id = H5Fopen(fname_ref, H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t group   = H5Gopen(file_id, "kspace",H5P_DEFAULT);
      H5LTread_dataset_float(group,"deltax_unfiltered",array_temp); cudaMemcpy(deltax_unfiltered_device,array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(group,"stars_unfiltered", array_temp); cudaMemcpy(stars_unfiltered_device, array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(group,"sfr_unfiltered",   array_temp); cudaMemcpy(sfr_unfiltered_device,   array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5Gclose(group);
      H5Fclose(file_id);
      free(array_temp);
    }

    // create working copies of the k-space grids
    cudaMemcpy(deltax_filtered_device,deltax_unfiltered_device,sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    cudaMemcpy(stars_filtered_device, stars_unfiltered_device, sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);
    cudaMemcpy(sfr_filtered_device,   sfr_unfiltered_device,   sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice);

    // Perform convolution
    int local_ix_start = 0; // make this = (int)(slabs_ix_start[mpi_rank]);
    if(!flag_last_filter_step){
       filter_gpu<<<grid_complex,threads>>>(deltax_filtered_device,local_nix,ReionGridDim,local_ix_start,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType);
       filter_gpu<<<grid_complex,threads>>>(stars_filtered_device, local_nix,ReionGridDim,local_ix_start,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType);
       filter_gpu<<<grid_complex,threads>>>(sfr_filtered_device,   local_nix,ReionGridDim,local_ix_start,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType);
    }

    if (validation_output && i_R==1 || !strcmp(fname,fname_full_dump))
    {
        // prepare output file
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "deltax_filtered", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "stars_filtered", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "sfr_filtered", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Gclose(group);
        H5Fclose(file_id);
    }
    if(false){ // good from here down
      float *array_temp = (float *)malloc(sizeof(float)*(slab_n_complex*2));
      hid_t file_id = H5Fopen(fname_ref, H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t group   = H5Gopen(file_id, "kspace",H5P_DEFAULT);
      H5LTread_dataset_float(group,"deltax_filtered",array_temp); cudaMemcpy(deltax_filtered_device,array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(group,"stars_filtered", array_temp); cudaMemcpy(stars_filtered_device, array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(group,"sfr_filtered",   array_temp); cudaMemcpy(sfr_filtered_device,   array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5Gclose(group);
      H5Fclose(file_id);
      free(array_temp);
    }

    // inverse fourier transform back to real space
#ifdef GPU_HYBRID
    cudaMemcpy(deltax_filtered,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    fftwf_execute(plan_deltax);
    cudaMemcpy(deltax_filtered_device,deltax_filtered,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);

    cudaMemcpy(stars_filtered,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    fftwf_execute(plan_stars);
    cudaMemcpy(stars_filtered_device,stars_filtered,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);

    cudaMemcpy(sfr_filtered,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    fftwf_execute(plan_sfr);
    cudaMemcpy(sfr_filtered_device,sfr_filtered,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
#else
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
#endif

    if (validation_output && i_R==1 || !strcmp(fname,fname_full_dump))
    {
        // prepare output file
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax_filtered_ift", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "stars_filtered_ift", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "sfr_filtered_ift", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Fclose(file_id);
    }
    if(false){
      hid_t file_id = H5Fopen(fname_ref, H5F_ACC_RDONLY, H5P_DEFAULT);
      float *array_temp = (float *)malloc(sizeof(float)*(2*slab_n_complex));
      H5LTread_dataset_float(file_id,"deltax_filtered_ift",array_temp);cudaMemcpy(deltax_filtered_device,array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(file_id,"stars_filtered_ift", array_temp);cudaMemcpy(stars_filtered_device, array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(file_id,"sfr_filtered_ift",   array_temp);cudaMemcpy(sfr_filtered_device,   array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      free(array_temp);
      H5Fclose(file_id);
    }

    // Perform sanity checks to account for aliasing effects
    sanity_check_aliasing<<<grid_real,threads>>>(deltax_filtered_device,ReionGridDim,local_ix_start,slab_n_real,-1.f + REL_TOL);
    sanity_check_aliasing<<<grid_real,threads>>>(stars_filtered_device, ReionGridDim,local_ix_start,slab_n_real,0.);
    sanity_check_aliasing<<<grid_real,threads>>>(sfr_filtered_device,   ReionGridDim,local_ix_start,slab_n_real,0.);

    if (validation_output && i_R==1 || !strcmp(fname,fname_full_dump))
    {
        // prepare output file
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax_checked", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "stars_checked", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "sfr_checked", 1, (hsize_t []){slab_n_complex * 2}, array_temp);
        free(array_temp);
        H5Fclose(file_id);
    }
    if(false){ 
      hid_t file_id = H5Fopen(fname_ref, H5F_ACC_RDONLY, H5P_DEFAULT);
      float *array_temp = (float *)malloc(sizeof(float)*(2*slab_n_complex));
      H5LTread_dataset_float(file_id,"deltax_checked",array_temp);cudaMemcpy(deltax_filtered_device,array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(file_id,"stars_checked", array_temp);cudaMemcpy(stars_filtered_device, array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(file_id,"sfr_checked",   array_temp);cudaMemcpy(sfr_filtered_device,   array_temp,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice);
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
        local_ix_start,
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

    if (validation_output && i_R==1 || !strcmp(fname,fname_full_dump))
    {
        // prepare output file
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*slab_n_real);
        cudaMemcpy(array_temp,xH_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "xH", 1, (hsize_t []){slab_n_real}, array_temp);
        cudaMemcpy(array_temp,r_bubble_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "r_bubble", 1, (hsize_t []){slab_n_real}, array_temp);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax", 1, (hsize_t []){slab_n_real}, array_temp);
        if(flag_ReionUVBFlag){
            cudaMemcpy(array_temp,J_21_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
            H5LTmake_dataset_float(file_id, "J_21", 1, (hsize_t []){slab_n_real}, array_temp);
        }
        free(array_temp);
        array_temp = (float *)malloc(sizeof(float)*(2*slab_n_complex));
        cudaMemcpy(array_temp,J_21_at_ionization_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, (hsize_t []){slab_n_real}, array_temp);
        cudaMemcpy(array_temp,z_at_ionization_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "z_at_ionization", 1, (hsize_t []){slab_n_real}, array_temp);
        free(array_temp);
        H5Fclose(file_id);
    }
    if(false){
      hid_t file_id = H5Fopen(fname_ref, H5F_ACC_RDONLY, H5P_DEFAULT);
      float *array_temp = (float *)malloc(sizeof(float)*(slab_n_real));
      H5LTread_dataset_float(file_id,"xH",      array_temp);cudaMemcpy(xH_device,      array_temp,sizeof(float)*slab_n_real,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(file_id,"r_bubble",array_temp);cudaMemcpy(r_bubble_device,array_temp,sizeof(float)*slab_n_real,cudaMemcpyHostToDevice);
      if(flag_ReionUVBFlag){
        H5LTread_dataset_float(file_id,"J_21",  array_temp);cudaMemcpy(J_21_device,array_temp,sizeof(float)*slab_n_real,cudaMemcpyHostToDevice);
      }
      H5LTread_dataset_float(file_id,"J_21_at_ionization",array_temp);cudaMemcpy(J_21_at_ionization_device,array_temp,sizeof(float)*slab_n_real,cudaMemcpyHostToDevice);
      H5LTread_dataset_float(file_id,"z_at_ionization",   array_temp);cudaMemcpy(z_at_ionization_device,   array_temp,sizeof(float)*slab_n_real,cudaMemcpyHostToDevice);
      free(array_temp);
      H5Fclose(file_id);
    }

    R /= ReionDeltaRFactor;
  }

#ifdef GPU_HYBRID
    fftwf_destroy_plan(plan_deltax);
    fftwf_destroy_plan(plan_stars);
    fftwf_destroy_plan(plan_sfr);
#else
    cufftDestroy(plan);
#endif

  // Perform device -> host transfer
  cudaMemcpy(xH,                xH_device,                sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost);
  cudaMemcpy(r_bubble,          r_bubble_device,          sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost);
  if(flag_ReionUVBFlag)
     cudaMemcpy(J_21,           J_21_device,              sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost);
  cudaMemcpy(z_at_ionization,   z_at_ionization_device,   sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost);
  cudaMemcpy(J_21_at_ionization,J_21_at_ionization_device,sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost);
  cudaMemcpy(deltax,            deltax_filtered_device,   sizeof(float) * 2* slab_n_complex, cudaMemcpyDeviceToHost);

  // Clean-up device
  cudaFree(deltax_unfiltered_device);
  cudaFree(stars_unfiltered_device);
  cudaFree(sfr_unfiltered_device);
  cudaFree(deltax_filtered_device);
  cudaFree(stars_filtered_device);
  cudaFree(sfr_filtered_device);
  cudaFree(xH_device);
  cudaFree(r_bubble_device);
  cudaFree(z_at_ionization_device);
  cudaFree(J_21_at_ionization_device);
  if(flag_ReionUVBFlag)
     cudaFree(J_21_device);

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

    H5LTmake_dataset_float(file_id, "xH",                 1, (hsize_t []){slab_n_real}, xH);
    H5LTmake_dataset_float(file_id, "z_at_ionization",    1, (hsize_t []){slab_n_real}, z_at_ionization);
    H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, (hsize_t []){slab_n_real}, J_21_at_ionization);

    H5LTset_attribute_double(file_id, "/", "volume_weighted_global_xH", volume_weighted_global_xH, 1);
    H5LTset_attribute_double(file_id, "/", "mass_weighted_global_xH",   mass_weighted_global_xH, 1);

    H5Fclose(file_id);
  }
}

