#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "meraxes.h"

void _find_HII_bubbles_gpu(double redshift,const bool flag_write_validation_output){
  // Fetch needed things from run_globals
  const int    mpi_rank            = run_globals.mpi_rank;
  const double box_size            = run_globals.params.BoxSize;
  const int    ReionGridDim        = run_globals.params.ReionGridDim;
  const int    flag_ReionUVBFlag   = run_globals.params.ReionUVBFlag;
  const double ReionEfficiency     = run_globals.params.physics.ReionEfficiency;
  const double ReionNionPhotPerBary= run_globals.params.physics.ReionNionPhotPerBary;
  const double UnitLength_in_cm    = run_globals.units.UnitLength_in_cm;
  const double UnitMass_in_g       = run_globals.units.UnitMass_in_g;
  const double UnitTime_in_s       = run_globals.units.UnitTime_in_s;
  const double ReionRBubbleMax     = run_globals.params.physics.ReionRBubbleMax;
  const double ReionRBubbleMin     = run_globals.params.physics.ReionRBubbleMin;
  const double ReionDeltaRFactor   = run_globals.params.ReionDeltaRFactor;
  const double ReionGammaHaloBias  = run_globals.params.physics.ReionGammaHaloBias;
  const double ReionAlphaUV        = run_globals.params.physics.ReionAlphaUV;
  const double ReionEscapeFrac     = run_globals.params.physics.ReionEscapeFrac;
  // grid parameters
  const ptrdiff_t *slabs_nix       = run_globals.reion_grids.slab_nix;
  const ptrdiff_t *slabs_n_complex = run_globals.reion_grids.slab_n_complex;
  const ptrdiff_t *slabs_ix_start  = run_globals.reion_grids.slab_ix_start;
  const int local_nix              = (int)slabs_nix[mpi_rank];
  const int local_ix_start         = (int)slabs_ix_start[mpi_rank];
  const int slab_n_complex         = (int)(slabs_n_complex[mpi_rank]);
  const int slab_n_real            = local_nix * ReionGridDim * ReionGridDim;
  // input grids
  float *deltax   = run_globals.reion_grids.deltax;  // real & padded
  float *stars    = run_globals.reion_grids.stars;   // real & padded
  float *sfr      = run_globals.reion_grids.sfr;     // real & padded
  // preallocated grids
  float   *J_21            = run_globals.reion_grids.J_21;     // real
  float   *r_bubble        = run_globals.reion_grids.r_bubble; // real
  // output grids
  float *xH                = run_globals.reion_grids.xH;                 // real
  float *z_at_ionization   = run_globals.reion_grids.z_at_ionization;    // real
  float *J_21_at_ionization= run_globals.reion_grids.J_21_at_ionization; // real
  // output values
  double *volume_weighted_global_xH=&(run_globals.reion_grids.volume_weighted_global_xH);
  double *mass_weighted_global_xH  =&(run_globals.reion_grids.mass_weighted_global_xH);
  // a few needed constants
  const double pixel_volume         = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
  const double total_n_cells        = pow((double)ReionGridDim, 3);
  const double inv_total_n_cells    = 1.f/total_n_cells;

  const hsize_t dset_real_size={slab_n_real};
  const hsize_t dset_cplx_size={2*slab_n_complex};

  if (flag_write_validation_output)
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

    H5LTmake_dataset_float(file_id, "deltax",             1, &dset_cplx_size, deltax);
    H5LTmake_dataset_float(file_id, "stars",              1, &dset_cplx_size, stars);
    H5LTmake_dataset_float(file_id, "sfr",                1, &dset_cplx_size, sfr);
    H5LTmake_dataset_float(file_id, "z_at_ionization",    1, &dset_real_size, z_at_ionization);
    H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, &dset_real_size, J_21_at_ionization);

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
  try{
    throw_on_cuda_error(cudaMalloc((void**)&deltax_unfiltered_device, sizeof(cufftComplex)*slab_n_complex),meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&stars_unfiltered_device,  sizeof(cufftComplex)*slab_n_complex),meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&sfr_unfiltered_device,    sizeof(cufftComplex)*slab_n_complex),meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&deltax_filtered_device,   sizeof(cufftComplex)*slab_n_complex),meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&stars_filtered_device,    sizeof(cufftComplex)*slab_n_complex),meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&sfr_filtered_device,      sizeof(cufftComplex)*slab_n_complex),meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&xH_device,                sizeof(float)*slab_n_real),          meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&r_bubble_device,          sizeof(float)*slab_n_real),          meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&z_at_ionization_device,   sizeof(float)*slab_n_real),          meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&J_21_at_ionization_device,sizeof(float)*slab_n_real),          meraxes_cuda_exception::MALLOC);
    if(flag_ReionUVBFlag)
       throw_on_cuda_error(cudaMalloc((void**)&J_21_device,           sizeof(float)*slab_n_real),meraxes_cuda_exception::MALLOC);
    // Throw an error if another rank has thrown an error
    throw_on_global_error(); 
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }

  // If we're not using CUFFT, do the forward FFT first, before sending it to the device
#ifndef USE_CUFFT
  // The following are only needed if we are using FFTW
  const MPI_Comm mpi_comm = run_globals.mpi_comm;
  Complex *deltax_filtered = (Complex *)run_globals.reion_grids.deltax_filtered;// complex TODO: Check the consistancy of Complex instances
  Complex *stars_filtered  = (Complex *)run_globals.reion_grids.stars_filtered; // complex TODO: Check the consistancy of Complex instances
  Complex *sfr_filtered    = (Complex *)run_globals.reion_grids.sfr_filtered;   // complex TODO: Check the consistancy of Complex instances

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

  // Perform host -> device transfer of input grids (note that these grids are k-space if we are using FFTW
  //    but are real-space if we are using CUFFT.  They will be transformed once on the device if the latter.
  try{
    throw_on_cuda_error(cudaMemcpy(deltax_unfiltered_device, deltax,            sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy(stars_unfiltered_device,  stars,             sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy(sfr_unfiltered_device,    sfr,               sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy(z_at_ionization_device,   z_at_ionization,   sizeof(float)*slab_n_real,     cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy(J_21_at_ionization_device,J_21_at_ionization,sizeof(float)*slab_n_real,     cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);
    // Throw an error if another rank has thrown an error
    throw_on_global_error(); 
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }

  // If we're using CUFFT, perform the forward FFT now that the data is on the device
#ifdef USE_CUFFT
  // Initialize cuFFT
  cufftHandle plan;
  try{
    throw_on_cuFFT_error(cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_R2C),meraxes_cuda_exception::CUFFT_CREATE_PLAN);
    throw_on_cuFFT_error(cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL),           meraxes_cuda_exception::CUFFT_SET_COMPATIBILITY);

    // Perform FFTs
    throw_on_cuFFT_error(cufftExecR2C(plan,(cufftReal *)deltax_unfiltered_device,deltax_unfiltered_device),meraxes_cuda_exception::CUFFT_R2C);
    throw_on_cuFFT_error(cufftExecR2C(plan,(cufftReal *)stars_unfiltered_device,stars_unfiltered_device),  meraxes_cuda_exception::CUFFT_R2C);
    throw_on_cuFFT_error(cufftExecR2C(plan,(cufftReal *)sfr_unfiltered_device,sfr_unfiltered_device),      meraxes_cuda_exception::CUFFT_R2C);

    // Clean-up 
    throw_on_cuFFT_error(cufftDestroy(plan),meraxes_cuda_exception::CUFFT_PLAN_DESTROY);

    // Make sure that the device has synchronized
    throw_on_cuda_error(cudaThreadSynchronize(),meraxes_cuda_exception::SYNC);

    // Throw an error if another rank has thrown an error
    throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }
#endif

  if (flag_write_validation_output)
  {
    // prepare output file
    char fname[STRLEN];
    sprintf(fname, "validation_output-core%03d-z%.2f.h5", mpi_rank, redshift);
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
    cudaMemcpy(array_temp,deltax_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    H5LTmake_dataset_float(group, "deltax", 1, &dset_cplx_size, array_temp);
    cudaMemcpy(array_temp,stars_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    H5LTmake_dataset_float(group, "stars", 1, &dset_cplx_size, array_temp);
    cudaMemcpy(array_temp,sfr_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
    H5LTmake_dataset_float(group, "sfr", 1, &dset_cplx_size, array_temp);
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
  int threads      = run_globals.gpu->n_threads;  
  int grid_complex = (slab_n_complex+(threads-1))/threads;
  int grid_real    = (slab_n_real   +(threads-1))/threads;

  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  try{
    throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(deltax_unfiltered_device,inv_total_n_cells,slab_n_complex)),meraxes_cuda_exception::KERNEL_CMPLX_AX);
    throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(stars_unfiltered_device, inv_total_n_cells,slab_n_complex)),meraxes_cuda_exception::KERNEL_CMPLX_AX);
    throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(sfr_unfiltered_device,   inv_total_n_cells,slab_n_complex)),meraxes_cuda_exception::KERNEL_CMPLX_AX);
    check_thread_sync(meraxes_cuda_exception::KERNEL_CMPLX_AX);
    // Throw an error if another rank has thrown an error
    throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }

  // Initialize a few of the output grids
  try{
    throw_on_kernel_error((set_array_gpu<<<grid_real,threads>>>(xH_device,      slab_n_real,1.f)),meraxes_cuda_exception::KERNEL_SET_ARRAY);
    throw_on_kernel_error((set_array_gpu<<<grid_real,threads>>>(r_bubble_device,slab_n_real,0.f)),meraxes_cuda_exception::KERNEL_SET_ARRAY);
    if(flag_ReionUVBFlag)
        throw_on_kernel_error((set_array_gpu<<<grid_real,threads>>>(J_21_device,slab_n_real,0.f)),meraxes_cuda_exception::KERNEL_SET_ARRAY);
    check_thread_sync(meraxes_cuda_exception::KERNEL_SET_ARRAY);
    // Throw an error if another rank has thrown an error
    throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  double cell_length_factor = L_FACTOR;
  if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0;

  // Initialize inverse FFTs
#ifdef USE_CUFFT
  try{
      throw_on_cuFFT_error(cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_C2R),meraxes_cuda_exception::CUFFT_C2R);
      throw_on_cuFFT_error(cufftSetCompatibilityMode(plan,CUFFT_COMPATIBILITY_FFTW_ALL),           meraxes_cuda_exception::CUFFT_SET_COMPATIBILITY);
      // Throw an error if another rank has thrown an error
      throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }
#else
  fftwf_plan plan_deltax = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, mpi_comm, FFTW_ESTIMATE);
  fftwf_plan plan_stars  = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)stars_filtered,  (float *)stars_filtered,  mpi_comm, FFTW_ESTIMATE);
  fftwf_plan plan_sfr    = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)sfr_filtered,    (float *)sfr_filtered,    mpi_comm, FFTW_ESTIMATE);
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
    if (flag_write_validation_output && (i_R==1 || !strcmp(fname,fname_full_dump)))
    {
        // prepare output file
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "deltax_unfiltered", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,stars_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "stars_unfiltered", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,sfr_unfiltered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "sfr_unfiltered", 1, &dset_cplx_size, array_temp);
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
    try{
        throw_on_cuda_error(cudaMemcpy(deltax_filtered_device,deltax_unfiltered_device,sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice),meraxes_cuda_exception::MEMCPY);
        throw_on_cuda_error(cudaMemcpy(stars_filtered_device, stars_unfiltered_device, sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice),meraxes_cuda_exception::MEMCPY);
        throw_on_cuda_error(cudaMemcpy(sfr_filtered_device,   sfr_unfiltered_device,   sizeof(Complex) * slab_n_complex,cudaMemcpyDeviceToDevice),meraxes_cuda_exception::MEMCPY);
        // Throw an error if another rank has thrown an error
        throw_on_global_error();
    }
    catch(const meraxes_cuda_exception e){
        e.process_exception();
    }

    // Perform convolution
    if(!flag_last_filter_step){
        try{
            throw_on_kernel_error((filter_gpu<<<grid_complex,threads>>>(deltax_filtered_device,local_nix,ReionGridDim,local_ix_start,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType)),meraxes_cuda_exception::KERNEL_FILTER);
            throw_on_kernel_error((filter_gpu<<<grid_complex,threads>>>(stars_filtered_device, local_nix,ReionGridDim,local_ix_start,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType)),meraxes_cuda_exception::KERNEL_FILTER);
            throw_on_kernel_error((filter_gpu<<<grid_complex,threads>>>(sfr_filtered_device,   local_nix,ReionGridDim,local_ix_start,slab_n_complex,R,box_size,run_globals.params.ReionRtoMFilterType)),meraxes_cuda_exception::KERNEL_FILTER);
            check_thread_sync(meraxes_cuda_exception::KERNEL_FILTER);
            // Throw an error if another rank has thrown an error
            throw_on_global_error();
        }
        catch(const meraxes_cuda_exception e){
            e.process_exception();
        }
    }

    if (flag_write_validation_output && (i_R==1 || !strcmp(fname,fname_full_dump)))
    {
        // prepare output file
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group = H5Gcreate(file_id, "kspace", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "deltax_filtered", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "stars_filtered", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(group, "sfr_filtered", 1, &dset_cplx_size, array_temp);
        free(array_temp);
        H5Gclose(group);
        H5Fclose(file_id);
    }
    if(false){ 
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
    try{
#ifdef USE_CUFFT
        throw_on_cuFFT_error(cufftExecC2R(plan,(cufftComplex *)deltax_filtered_device, (cufftReal *)deltax_filtered_device),meraxes_cuda_exception::CUFFT_C2R);
        throw_on_cuFFT_error(cufftExecC2R(plan,(cufftComplex *)stars_filtered_device, (cufftReal *)stars_filtered_device),  meraxes_cuda_exception::CUFFT_C2R);
        throw_on_cuFFT_error(cufftExecC2R(plan,(cufftComplex *)sfr_filtered_device, (cufftReal *)sfr_filtered_device),      meraxes_cuda_exception::CUFFT_C2R);
#else
        throw_on_cuda_error(cudaMemcpy(deltax_filtered,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
        fftwf_execute(plan_deltax);
        throw_on_cuda_error(cudaMemcpy(deltax_filtered_device,deltax_filtered,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);

        throw_on_cuda_error(cudaMemcpy(stars_filtered,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
        fftwf_execute(plan_stars);
        throw_on_cuda_error(cudaMemcpy(stars_filtered_device,stars_filtered,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);

        throw_on_cuda_error(cudaMemcpy(sfr_filtered,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
        fftwf_execute(plan_sfr);
        throw_on_cuda_error(cudaMemcpy(sfr_filtered_device,sfr_filtered,sizeof(float)*2*slab_n_complex,cudaMemcpyHostToDevice),meraxes_cuda_exception::MEMCPY);
#endif
        // Throw an error if another rank has thrown an error
        throw_on_global_error();
    }
    catch(const meraxes_cuda_exception e){
        e.process_exception();
    }

    if (flag_write_validation_output && (i_R==1 || !strcmp(fname,fname_full_dump)))
    {
        // prepare output file
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax_filtered_ift", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "stars_filtered_ift", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "sfr_filtered_ift", 1, &dset_cplx_size, array_temp);
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
    try{
        throw_on_kernel_error((sanity_check_aliasing<<<grid_real,threads>>>(deltax_filtered_device,ReionGridDim,local_ix_start,slab_n_real,-1.f + REL_TOL)),meraxes_cuda_exception::KERNEL_CHECK);
        throw_on_kernel_error((sanity_check_aliasing<<<grid_real,threads>>>(stars_filtered_device, ReionGridDim,local_ix_start,slab_n_real,0.)),meraxes_cuda_exception::KERNEL_CHECK);
        throw_on_kernel_error((sanity_check_aliasing<<<grid_real,threads>>>(sfr_filtered_device,   ReionGridDim,local_ix_start,slab_n_real,0.)),meraxes_cuda_exception::KERNEL_CHECK);
        check_thread_sync(meraxes_cuda_exception::KERNEL_CHECK);
        // Throw an error if another rank has thrown an error
        throw_on_global_error();
    }
    catch(const meraxes_cuda_exception e){
        e.process_exception();
    }

    if (flag_write_validation_output && (i_R==1 || !strcmp(fname,fname_full_dump)))
    {
        // prepare output file
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*2*slab_n_complex);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax_checked", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,stars_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "stars_checked", 1, &dset_cplx_size, array_temp);
        cudaMemcpy(array_temp,sfr_filtered_device,sizeof(float)*2*slab_n_complex,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "sfr_checked", 1, &dset_cplx_size, array_temp);
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

    try{
        throw_on_kernel_error((find_HII_bubbles_gpu_main_loop<<<grid_real,threads>>>(
            redshift, 
            slab_n_real, 
            flag_last_filter_step, 
            flag_ReionUVBFlag, 
            ReionGridDim, 
            local_ix_start,
            R, 
            RtoM(R), 
            ReionEfficiency, 
            inv_pixel_volume, 
            J_21_aux_constant, 
            ReionGammaHaloBias, 
            xH_device,
            J_21_device,
            r_bubble_device,
            J_21_at_ionization_device,
            z_at_ionization_device,
            deltax_filtered_device,
            stars_filtered_device,
            sfr_filtered_device)),meraxes_cuda_exception::KERNEL_MAIN_LOOP);
        check_thread_sync(meraxes_cuda_exception::KERNEL_MAIN_LOOP);
        // Throw an error if another rank has thrown an error
        throw_on_global_error();
    }
    catch(const meraxes_cuda_exception e){
        e.process_exception();
    }

    if (flag_write_validation_output && (i_R==1 || !strcmp(fname,fname_full_dump)))
    {
        // prepare output file
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
        float *array_temp = (float *)malloc(sizeof(float)*slab_n_real);
        cudaMemcpy(array_temp,xH_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "xH", 1, &dset_real_size, array_temp);
        cudaMemcpy(array_temp,r_bubble_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "r_bubble", 1, &dset_real_size, array_temp);
        cudaMemcpy(array_temp,deltax_filtered_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "deltax", 1, &dset_real_size, array_temp);
        if(flag_ReionUVBFlag){
            cudaMemcpy(array_temp,J_21_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
            H5LTmake_dataset_float(file_id, "J_21", 1, &dset_real_size, array_temp);
        }
        free(array_temp);
        array_temp = (float *)malloc(sizeof(float)*(2*slab_n_complex));
        cudaMemcpy(array_temp,J_21_at_ionization_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "J_21_at_ionization", 1, &dset_real_size, array_temp);
        cudaMemcpy(array_temp,z_at_ionization_device,sizeof(float)*slab_n_real,cudaMemcpyDeviceToHost);
        H5LTmake_dataset_float(file_id, "z_at_ionization", 1, &dset_real_size, array_temp);
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

  // Clean-up FFT plan(s)
#ifdef USE_CUFFT
  try{
    throw_on_cuFFT_error(cufftDestroy(plan),meraxes_cuda_exception::CUFFT_PLAN_DESTROY);
    // Throw an error if another rank has thrown an error
    throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }
#else
  fftwf_destroy_plan(plan_deltax);
  fftwf_destroy_plan(plan_stars);
  fftwf_destroy_plan(plan_sfr);
#endif

  // Perform device -> host transfer
  try{
    throw_on_cuda_error(cudaMemcpy((void *)xH,                (void *)xH_device,                sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy((void *)r_bubble,          (void *)r_bubble_device,          sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
    if(flag_ReionUVBFlag)
         throw_on_cuda_error(cudaMemcpy((void *)J_21,           (void *)J_21_device,              sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy((void *)z_at_ionization,   (void *)z_at_ionization_device,   sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy((void *)J_21_at_ionization,(void *)J_21_at_ionization_device,sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(cudaMemcpy((void *)deltax,            (void *)deltax_filtered_device,   sizeof(float) * 2* slab_n_complex, cudaMemcpyDeviceToHost),meraxes_cuda_exception::MEMCPY);
    // Throw an error if another rank has thrown an error
    throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }

  // Clean-up device
  try{
    throw_on_cuda_error(cudaFree(deltax_unfiltered_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(stars_unfiltered_device),  meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(sfr_unfiltered_device),    meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(deltax_filtered_device),   meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(stars_filtered_device),    meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(sfr_filtered_device),      meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(xH_device),                meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(r_bubble_device),          meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(z_at_ionization_device),   meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(J_21_at_ionization_device),meraxes_cuda_exception::FREE);
    if(flag_ReionUVBFlag)
       throw_on_cuda_error(cudaFree(J_21_device),meraxes_cuda_exception::FREE);
    // Throw an error if another rank has thrown an error
    throw_on_global_error();
  }
  catch(const meraxes_cuda_exception e){
      e.process_exception();
  }

  // Find the volume and mass weighted neutral fractions
  // TODO: The deltax grid will have rounding errors from forward and reverse
  //       FFT. Should cache deltax slabs prior to ffts and reuse here.
  *volume_weighted_global_xH = 0.0;
  *mass_weighted_global_xH   = 0.0;
  double mass_weight         = 0.0;

  // Calculate neutral fractions
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
  *volume_weighted_global_xH *= inv_total_n_cells;
  *mass_weighted_global_xH   /= mass_weight;
}

