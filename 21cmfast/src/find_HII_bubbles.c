#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "read_grids.c"
#include "read_galaxies.c"

/*
USAGE: find_HII_bubbles [-p <NUM THREADS>] <snapshot> [<ionization efficiency factor zeta>]

Program FIND_HII_BUBBLES reads in the halo list file from
../Output_files/Halo_lists/updated_halos_z%06.2f_%i_%.0fMpc,
where the relavant parameters are taken from INIT_PARAMS.H.
For each box pixel, it cycles through various bubble radii
(taken from ANAL_PARAMS.H), until it finds the largest 
radius such that the enclosed collapsed mass fraction 
(obtained by summing masses from the halo list file of
halos whose centers are within the bubble, or by taking 
the mean collapsed mass from conditional press-schechter)
is larger than 1/HII_EFF_FACTOR (from ANAL_PARAMS.H)
or the optional command line efficiency parameter 

If optional efficiency parameter is not passed, then the 
value of HII_EFF_FACTOR in ANAL_PARAMS.H is used.

NOTE: the optional argument of thread number including the -p flag, MUST
be the first two arguments.  If these are omitted, num_threads defaults
to NUMCORES in INIT_PARAMS.H

See ANAL_PARAMS.H for updated parameters and algorithms!

Author: Andrei Mesinger
Date: 01/10/07
*/

static int write_xH(char *meraxes_fname, int snapshot, float *xH, double global_xH, float ION_EFF_FACTOR, float M_MIN, float MFP)
{
  hid_t   file_id;
  hid_t   group_id;
  hsize_t ds_dims[1];
  char    target_group[256];
  int     tmp_int;
  float   tmp_float;

  if (!(file_id = H5Fopen(meraxes_fname, H5F_ACC_RDWR, H5P_DEFAULT))){
    fprintf(stderr, "find_HII_bubbles:: ERROR: unable to open file %s for writting!\n", meraxes_fname);
    fprintf(LOG, "find_HII_bubbles:: ERROR: unable to open file %s for writting!\n", meraxes_fname);
    global_xH = -1;
  }else 
  {
    fprintf(LOG, "find_HII_bubbles:: Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, meraxes_fname);
    fprintf(stderr, "find_HII_bubbles:: Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, meraxes_fname);
    fflush(LOG);

    sprintf(target_group, "Snap%03d", snapshot);
    group_id = H5Gopen(file_id, target_group, H5P_DEFAULT);
    ds_dims[0] = HII_TOT_NUM_PIXELS;
    H5LTmake_dataset(group_id, "xH_grid", 1, ds_dims, H5T_NATIVE_FLOAT, xH);
    H5LTset_attribute_double(group_id, "xH_grid", "global_xH", &global_xH, 1);
    tmp_float = ION_EFF_FACTOR;
    H5LTset_attribute_float(group_id, "xH_grid", "ion_eff_factor", &tmp_float, 1);
    tmp_float = M_MIN;
    H5LTset_attribute_float(group_id, "xH_grid", "Mvir_min", &tmp_float, 1);
    tmp_float = MFP;
    H5LTset_attribute_float(group_id, "xH_grid", "MFP", &tmp_float, 1);
    tmp_float = BOX_LEN;
    H5LTset_attribute_float(group_id, "xH_grid", "box_len", &tmp_float, 1);
    tmp_int = HII_FILTER;
    H5LTset_attribute_int(group_id, "xH_grid", "HII_filter", &tmp_int, 1);
    tmp_int = HII_DIM;
    H5LTset_attribute_int(group_id, "xH_grid", "HII_dim", &tmp_int, 1);

    switch(FIND_BUBBLE_ALGORITHM)
    {
      case 2:
        H5LTset_attribute_string(group_id, "xH_grid", "bubble_algorithm", "center");
        break;
      default:
        H5LTset_attribute_string(group_id, "xH_grid", "bubble_algorithm", "sphere");
    }

    H5Gclose(group_id);
    H5Fclose(file_id);
  }

  return global_xH;

}


int find_HII_bubbles(tocf_params_struct *params)
{

  char                filename[300];
  FILE               *F;
  FILE               *pPipe;
  float               REDSHIFT;
  float               mass;
  float               R;
  float               growth_factor;
  float               pixel_mass;
  float               cell_length_factor;
  float               ave_N_min_cell;
  float               ION_EFF_FACTOR;
  float               M_MIN;
  int                 x;
  int                 y;
  int                 z;
  int                 N_min_cell;
  int                 LAST_FILTER_STEP;
  int                 num_th;
  int                 arg_offset;
  int                 i                  = 0;
  int                 j;
  int                 k;
  unsigned long long  ct;
  unsigned long long  ion_ct;
  unsigned long long  sample_ct;
  float               f_coll_crit;
  float               pixel_volume;
  float               density_over_mean;
  float               erfc_num;
  float               erfc_denom;
  float               erfc_denom_cell;
  float               res_xH;
  float              *xH;
  float               TVIR_MIN;
  float               MFP;
  float               xHI_from_xrays;
  float               std_xrays;
  fftwf_complex      *M_coll_unfiltered;
  fftwf_complex      *M_coll_filtered;
  fftwf_complex      *deltax_unfiltered;
  fftwf_complex      *deltax_filtered;
  fftwf_complex      *xe_unfiltered;
  fftwf_complex      *xe_filtered;
  fftwf_plan          plan;
  double              global_xH;
  double              ave_xHI_xrays;
  double              ave_den;
  double              ST_over_PS;
  double              mean_f_coll_st;
  double              mean_f_coll_ps;
  double              f_coll;
  double              ave_fcoll;
  const gsl_rng_type *T;
  gsl_rng            *r;
  int                 snapshot;
  int                 corrected_snapshot;
  int                 status;
  char               *meraxes_fname;

  hid_t   file_id;
  char    target_group[256];

  FILE               *LOG;
  unsigned long long  SAMPLING_INTERVAL = (((unsigned long long)(HII_TOT_NUM_PIXELS/1.0e6)) + 1); //used to sample density field to compute mean collapsed fraction

  num_th = params->num_th;
  if (num_th<0)
    num_th = NUMCORES;
  fprintf(stderr, "find_HII_bubbles:: threading with %d threads\n", num_th);

  ION_EFF_FACTOR = params->ion_eff_factor;
  if (ION_EFF_FACTOR<0)
    ION_EFF_FACTOR = HII_EFF_FACTOR;

  TVIR_MIN       = params->tvir_min;
  if (TVIR_MIN<0)
    TVIR_MIN = ION_Tvir_MIN;

  MFP            = params->mean_free_path;
  if (MFP<0)
    MFP = R_BUBBLE_MAX;

  snapshot       = params->snapshot;
  meraxes_fname  = params->meraxes_fname;

  if (fftwf_init_threads()==0){
    fprintf(stderr, "find_HII_bubbles:: ERROR: problem initializing fftwf threads\nAborting\n.");
    return -1;
  }
  omp_set_num_threads(num_th);
  fftwf_plan_with_nthreads(num_th);    
  init_ps();
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  // Get the corrected (unsampled) snapshot number
  file_id = H5Fopen(meraxes_fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(target_group, "Snap%03d", snapshot);
  H5LTget_attribute_float(file_id, target_group, "Redshift", &REDSHIFT);
  H5LTget_attribute_int(file_id, target_group, "CorrectedSnap", &corrected_snapshot);
  H5Fclose(file_id);
  fprintf(stderr, "find_HII_bubbles:: Snapshot = %d (%d)-> redshift = %.2f\n", snapshot, corrected_snapshot, REDSHIFT);

  growth_factor = dicke(REDSHIFT);
  pixel_volume = pow(BOX_LEN/(float)HII_DIM, 3);
  pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM); 

  // This is for DM only:
  // f_coll_crit = 1/ION_EFF_FACTOR;
  // if we feed stellar mass instead of halo mass to the program then:
  f_coll_crit = 0.1/ION_EFF_FACTOR;
  // assuming that the star formation efficiency (f_* of eqn. 83 Furlanetto+ 2006) was previously assumed to be 10%

  cell_length_factor = L_FACTOR;
  // this parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations.  it probably doesn't matter for larger cell sizes.
  if (USE_HALO_FIELD && (FIND_BUBBLE_ALGORITHM==2)
      && ((BOX_LEN/(float)HII_DIM) < 1)){ // fairly arbitrary length based on 2 runs i did
    cell_length_factor = 1; 
  }

  //set the minimum source mass
  if ((TVIR_MIN > 0) && (ION_M_MIN > 0)){
    fprintf(stderr, "find_HII_bubbles:: You have to \"turn-off\" either the ION_M_MIN or \
        the ION_Tvir_MIN option in ANAL_PARAMS.H\nAborting...\n");
    free_ps(); return -1;
  }
  if (TVIR_MIN > 0){ // use the virial temperature for Mmin
    if (TVIR_MIN < 9.99999e3) // neutral IGM
      M_MIN = TtoM(REDSHIFT, TVIR_MIN, 1.22);
    else // ionized IGM
      M_MIN = TtoM(REDSHIFT, TVIR_MIN, 0.6);
  }
  else if (TVIR_MIN < 0){ // use the mass
    M_MIN = ION_M_MIN;
  }
  // check for WDM
  if (P_CUTOFF && ( M_MIN < M_J_WDM())){
    fprintf(stderr, "find_HII_bubbles:: The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
    M_MIN = M_J_WDM();
    fprintf(stderr, "find_HII_bubbles:: Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
  }


  // open log file
  sprintf(filename, "%s/HII_bubble_log_file_%03d.log", params->logfile_dir, snapshot);
  LOG = fopen(filename, "w");
  if (!LOG){
    fprintf(stderr, "find_HII_bubbles:: Error opening log file\nAborting...\n");
    fftwf_cleanup_threads();
    free_ps(); return -1;
  }

  // allocate memory for the neutral fraction box
  xH = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!xH){
    fprintf(stderr, "find_HII_bubbles:: Error allocating memory for xH box\nAborting...\n");
    fprintf(LOG, "find_HII_bubbles:: Error allocating memory for xH box\nAborting...\n");
    fclose(LOG); fftwf_cleanup_threads();
    free_ps(); return -1;
  }
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){    xH[ct] = 1;  }

  // lets check if we are going to bother with computing the inhmogeneous field at all...
  mean_f_coll_st = FgtrM_st(REDSHIFT, M_MIN);
  mean_f_coll_ps = FgtrM(REDSHIFT, M_MIN);
  if ((mean_f_coll_st/f_coll_crit < HII_ROUND_ERR) || (REDSHIFT > Z_HEAT_MAX)){ // way too small to ionize anything...
    fprintf(stderr, "find_HII_bubbles:: The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n I will just declare everything to be neutral\n", mean_f_coll_st, f_coll_crit);
    fprintf(LOG, "find_HII_bubbles:: The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n I will just declare everything to be neutral\n", mean_f_coll_st, f_coll_crit);

    if (USE_TS_IN_21CM  && (REDSHIFT < Z_HEAT_MAX)){ // use the x_e box
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
      if (!(F = fopen(filename, "rb"))){
        fprintf(stderr, "find_HII_bubbles:: find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
        fprintf(LOG, "find_HII_bubbles:: Unable to open x_e file at %s\nAborting...\n", filename);
        fclose(LOG); fftwf_free(xH);  fftwf_cleanup_threads();
        free_ps(); return -1;
      }
      for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        if (fread(&xH[ct], sizeof(float), 1, F)!=1){
          fprintf(stderr, "find_HII_bubbles:: Read error occured while reading xe box.\nAborting...\n");
          fprintf(LOG, "find_HII_bubbles:: Read error occured while reading xe box.\nAborting\n");
          fclose(LOG); fftwf_free(xH);  fftwf_cleanup_threads();
          free_ps(); fftwf_free(xe_filtered); fftwf_free(xe_unfiltered); return -1;
        }
        xH[ct] = 1-xH[ct]; // convert from x_e to xH
        global_xH += xH[ct];
      }
      fclose(F);
      global_xH /= (double)HII_TOT_NUM_PIXELS;
    }
    else{
      // find the neutral fraction
      init_heat();
      global_xH = 1 - xion_RECFAST(REDSHIFT, 0);;
      destruct_heat();
      for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        xH[ct] = global_xH;
      }
    }

    // print out the xH box
    global_xH = write_xH(meraxes_fname, snapshot, xH, global_xH, ION_EFF_FACTOR, M_MIN, MFP);
    fclose(F); fclose(LOG); fftwf_free(xH); fftwf_cleanup_threads();
    free_ps(); return (int) (global_xH * 100);
  }


  // see if we want to weigh by the number of neutral atoms from x-ray preheating
  if (USE_TS_IN_21CM){
    xe_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    xe_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (!xe_filtered || !xe_unfiltered){
      fprintf(stderr, "find_HII_bubbles:: Error allocating memory for xe boxes\nAborting...\n");
      fprintf(LOG, "find_HII_bubbles:: Error allocating memory for xe boxes\nAborting...\n");
      fclose(LOG); fftwf_free(xH);  fftwf_cleanup_threads();
      free_ps(); return -1;
    }

    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
    if (!(F = fopen(filename, "rb"))){
      fprintf(stderr, "find_HII_bubbles:: Unable to open x_e file at %s\nAborting...\n", filename);
      fprintf(LOG, "find_HII_bubbles:: Unable to open x_e file at %s\nAborting...\n", filename);
      fclose(LOG); fftwf_free(xH);  fftwf_cleanup_threads();
      free_ps(); fftwf_free(xe_filtered); fftwf_free(xe_unfiltered); return -1;
    }
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
        for (k=0; k<HII_DIM; k++){
          if (fread((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
            fprintf(stderr, "find_HII_bubbles:: Read error occured while reading xe box.\nAborting...\n");
            fprintf(LOG, "find_HII_bubbles:: Read error occured while reading xe box.\nAborting\n");
            fclose(LOG); fftwf_free(xH);  fftwf_cleanup_threads();
            free_ps(); fftwf_free(xe_filtered); fftwf_free(xe_unfiltered); return -1;
          }
        }
      }
    }
    fclose(F);
  }



  if (USE_HALO_FIELD){
    // allocate memory for the smoothed halo field
    M_coll_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    M_coll_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (!M_coll_unfiltered || !M_coll_filtered){
      fprintf(stderr, "find_HII_bubbles:: Error allocating memory for M_coll boxes\nAborting...\n");
      fprintf(LOG, "find_HII_bubbles.c:: Error allocating memory for M_coll boxes\nAborting...\n");
      fclose(LOG); fftwf_free(xH);  fftwf_cleanup_threads();
      free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return -1;
    }
    for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS; ct++){    *((float *)M_coll_unfiltered + ct) = 0;  }

    // read in the galaxy stellar mass field
    status = generate_stellarmass_field(meraxes_fname, snapshot, M_MIN, (float *)M_coll_unfiltered);
    if (status!=0)
    {
      fftwf_free(xH); fftwf_free(M_coll_unfiltered); fftwf_free(M_coll_filtered); fclose(LOG);
      fftwf_cleanup_threads();
      free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return -1;
    }

  } // end of the USE_HALO_FIELD option

  // read in the perturbed density field
  deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!deltax_unfiltered || !deltax_filtered){
    fprintf(stderr, "find_HII_bubbles:: Error allocating memory for deltax boxes\nAborting...\n");
    fprintf(LOG, "find_HII_bubbles:: Error allocating memory for deltax boxes\nAborting...\n");
    fftwf_free(xH); fclose(LOG);
    if (USE_HALO_FIELD){ fftwf_free(M_coll_unfiltered); fftwf_free(M_coll_filtered); }
    fftwf_cleanup_threads();
    free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return -1;
  }

  status = read_nbody_grid(params, corrected_snapshot, 0, (float *)deltax_unfiltered);
  if(status!=0)
  {
    fprintf(LOG, "find_HII_bubbles:: Read error occured while reading n_body grid.\n");
    fftwf_free(xH); fclose(LOG); fftwf_free(deltax_unfiltered); fftwf_free(deltax_filtered);
    if (USE_HALO_FIELD){ fftwf_free(M_coll_unfiltered); fftwf_free(M_coll_filtered); }
    fftwf_cleanup_threads(); fclose(F);
    free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return -1;
  }


  // do the fft to get the k-space M_coll field and deltax field
  fprintf(LOG, "find_HII_bubbles:: begin initial ffts, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
  fflush(LOG);
  if (USE_HALO_FIELD){
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)M_coll_unfiltered, (fftwf_complex *)M_coll_unfiltered, FFTW_ESTIMATE);
    fftwf_execute(plan);
  }
  if (USE_TS_IN_21CM){
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)xe_unfiltered, (fftwf_complex *)xe_unfiltered, FFTW_ESTIMATE);
    fftwf_execute(plan);    
  }
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
  //  real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
    if (USE_HALO_FIELD){  M_coll_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);}
    if (USE_TS_IN_21CM){ xe_unfiltered[ct] /= (float)HII_TOT_NUM_PIXELS;}
    deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
  }
  fprintf(LOG, "find_HII_bubbles:: end initial ffts, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
  fflush(LOG);


  // loop through the filter radii (in Mpc)
  erfc_denom_cell=1; //dummy value
  R=fmin(MFP, L_FACTOR*BOX_LEN);
  LAST_FILTER_STEP = 0;
  while (!LAST_FILTER_STEP){//(R > (cell_length_factor*BOX_LEN/(HII_DIM+0.0))){
    if ((R/DELTA_R_HII_FACTOR) <= (cell_length_factor*BOX_LEN/(float)HII_DIM)){
      LAST_FILTER_STEP = 1;
    }

    fprintf(LOG, "find_HII_bubbles:: before memcpy, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    if (USE_HALO_FIELD){    memcpy(M_coll_filtered, M_coll_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
    if (USE_TS_IN_21CM){    memcpy(xe_filtered, xe_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    fprintf(LOG, "find_HII_bubbles:: begin filter, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    if (USE_HALO_FIELD) { HII_filter(M_coll_filtered, HII_FILTER, R);}
    if (USE_TS_IN_21CM) { HII_filter(xe_filtered, HII_FILTER, R);}
    HII_filter(deltax_filtered, HII_FILTER, R);
    fprintf(LOG, "find_HII_bubbles:: end filter, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);

    // do the FFT to get the M_coll
    fprintf(LOG, "find_HII_bubbles:: begin fft with R=%f, clock=%06.2f\n", R, (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    if (USE_HALO_FIELD) {
      plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)M_coll_filtered, (float *)M_coll_filtered, FFTW_ESTIMATE);
      fftwf_execute(plan);
    }
    if (USE_TS_IN_21CM) {
      plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)xe_filtered, (float *)xe_filtered, FFTW_ESTIMATE);
      fftwf_execute(plan);
    }
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    fprintf(LOG, "find_HII_bubbles:: end fft with R=%f, clock=%06.2f\n", R, (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);


    /* Check if this is the last filtering scale.  If so, we don't need deltax_unfiltered anymore.
       We will re-read it to get the real-space field, which we will use to se the residual
       neutral fraction */
    ST_over_PS = 0;
    f_coll = 0;
    if (LAST_FILTER_STEP){
      status = read_nbody_grid(params, corrected_snapshot, 0, (float *)deltax_unfiltered);
      if(status!=0)
      {
        fprintf(LOG, "find_HII_bubbles:: Read error occured while reading n_body grid.\n");
        fftwf_free(xH); fclose(LOG); fftwf_free(deltax_unfiltered); fftwf_free(deltax_filtered); fftwf_free(M_coll_unfiltered); fftwf_free(M_coll_filtered);  fftwf_cleanup_threads(); fclose(F);
        free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return -1;
      }


      if (USE_HALO_FIELD){
        for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS; ct++){    *((float *)M_coll_unfiltered + ct) = 0;  }

        // read in the galaxy stellar mass field
        status = generate_stellarmass_field(meraxes_fname, snapshot, M_MIN, (float *)M_coll_unfiltered);
        if (status!=0)
        {
          fftwf_free(xH); fftwf_free(M_coll_unfiltered); fftwf_free(M_coll_filtered); fclose(LOG);
          fftwf_cleanup_threads();
          free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return -1;
        }
      } // end halo field option
      else{
        erfc_denom_cell = sqrt( 2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(RtoM(cell_length_factor*BOX_LEN/(float)HII_DIM)), 2) ) );
        if (erfc_denom_cell < 0){  // our filtering scale has become too small
          break;
        }

        // renormalize the collapse fraction so that the mean matches ST, 
        // since we are using the evolved (non-linear) density field
        sample_ct=0;
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct+=SAMPLING_INTERVAL){
          density_over_mean = 1.0 + *((float *)deltax_unfiltered + ct);
          erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
          f_coll += splined_erfc(erfc_num/erfc_denom_cell);	  
          sample_ct++;
        }
        f_coll /= (double) sample_ct;
        ST_over_PS = mean_f_coll_st/f_coll;
      }


      // and the spin temperature option is set
      if (USE_TS_IN_21CM){
        sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
        F = fopen(filename, "rb");
        for (i=0; i<HII_DIM; i++){
          for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
              fread((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F);
            }
          }
        }
        fclose(F);
      } // end re-reading the x_e box
    } // end if last filter step conditional statement

    // not the last filter step, and we operating on the density field
    else if (!USE_HALO_FIELD){
      fprintf(LOG, "find_HII_bubbles:: begin f_coll normalization if, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
      erfc_denom = sqrt( 2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(RtoM(R)), 2) ) );
      if (erfc_denom < 0){  // our filtering scale has become too small
        break;
      }

      // renormalize the collapse fraction so that the mean matches ST, 
      // since we are using the evolved (non-linear) density field
      sample_ct=0;
      for (ct=0; ct<HII_TOT_NUM_PIXELS; ct+=SAMPLING_INTERVAL){
        density_over_mean = 1.0 + *((float *)deltax_filtered + ct);
        erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
        f_coll += splined_erfc(erfc_num/erfc_denom);
        sample_ct++;
      }      
      f_coll /= (double) sample_ct;
      ST_over_PS = mean_f_coll_st/f_coll;
      fprintf(LOG, "find_HII_bubbles:: end f_coll normalization if, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
    }
    // fprintf(stderr, "Last filter %i, R_filter=%f, fcoll=%f, ST_over_PS=%f, mean normalized fcoll=%f\n", LAST_FILTER_STEP, R, f_coll, ST_over_PS, f_coll*ST_over_PS);

    /************  MAIN LOOP THROUGH THE BOX **************/
    fprintf(LOG, "find_HII_bubbles:: start of main loop scroll, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fprintf(stderr, "find_HII_bubbles:: start of main loop scroll, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    // now lets scroll through the filtered box
    ave_xHI_xrays = ave_den = ave_fcoll = std_xrays = 0;
    ion_ct=0;
    for (x=0; x<HII_DIM; x++){
      for (y=0; y<HII_DIM; y++){
        for (z=0; z<HII_DIM; z++){
          if (LAST_FILTER_STEP)
            density_over_mean = 1.0 + *((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z));
          else
            density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));

          /* check for aliasing which can occur for small R and small cell sizes,
             since we are using the analytic form of the window function for speed and simplicity */
          if (density_over_mean <= 0){
            // fprintf(LOG, "WARNING: aliasing during filtering step produced density n/<n> of %06.2f at cell (%i, %i, %i)\n Setting to 0\n", density_over_mean, x,y,z);
            density_over_mean = FRACT_FLOAT_ERR;
          }


          if (USE_HALO_FIELD){
            if (LAST_FILTER_STEP)
              f_coll = *((float *)M_coll_unfiltered + HII_R_FFT_INDEX(x,y,z)) / (RtoM(R)*density_over_mean);
            else
              f_coll = *((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) / (RtoM(R)*density_over_mean);
            f_coll *= (4/3.0)*PI*pow(R,3) / pixel_volume;
          }
          else{
            erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
            if (LAST_FILTER_STEP)
              f_coll = ST_over_PS * splined_erfc(erfc_num/erfc_denom_cell);
            else
              f_coll = ST_over_PS * splined_erfc(erfc_num/erfc_denom);
            //	      FgtrM_bias(REDSHIFT, M_MIN, density_over_mean-1, sigma_z0(RtoM(R)));
          }

          // adjust the denominator of the collapse fraction for the residual electron fraction in the neutral medium
          xHI_from_xrays = 1;
          if (USE_TS_IN_21CM){
            if (LAST_FILTER_STEP)
              xHI_from_xrays = (1 - *((float *)xe_unfiltered + HII_R_FFT_INDEX(x,y,z)));
            else
              xHI_from_xrays = (1 - *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)));

            if (xHI_from_xrays>1)
              xHI_from_xrays = 1;
            else if (xHI_from_xrays<0.001) // lets use a floor of 10^-3 for residual neutral fraction
              xHI_from_xrays = 0.001;

            f_coll /= xHI_from_xrays;
          }

          // check if ionized!
          if (f_coll > f_coll_crit){ // ionized

            if (FIND_BUBBLE_ALGORITHM == 2) // center method
              xH[HII_R_INDEX(x, y, z)] = 0;

            else if (FIND_BUBBLE_ALGORITHM == 1) // sphere method
              update_in_sphere(xH, HII_DIM, R/BOX_LEN, x/(HII_DIM+0.0), y/(HII_DIM+0.0), z/(HII_DIM+0.0));

            else{
              fprintf(stderr, "find_HII_bubbles:: Incorrect choice of find bubble algorithm: %i\nAborting...", FIND_BUBBLE_ALGORITHM);
              fprintf(LOG, "find_HII_bubbles:: Incorrect choice of find bubble algorithm: %i\nAborting...", FIND_BUBBLE_ALGORITHM);
              fflush(NULL);
              z=HII_DIM;y=HII_DIM,x=HII_DIM;R=0;
            }
          }

          /* check if this is the last filtering step.
             if so, assign partial ionizations to those cells which aren't fully ionized */
             else if (LAST_FILTER_STEP && (xH[HII_R_INDEX(x, y, z)] > TINY)){
               if (!USE_HALO_FIELD){
                 f_coll = ST_over_PS * splined_erfc(erfc_num/erfc_denom_cell);
                 if (f_coll>1) f_coll=1;
                 ave_N_min_cell = f_coll * pixel_mass*density_over_mean / M_MIN; // ave # of M_MIN halos in cell
                 if (ave_N_min_cell < N_POISSON){ 
                   // the collapsed fraction is too small, lets add poisson scatter in the halo number
                   N_min_cell = (int) gsl_ran_poisson(r, ave_N_min_cell);
                   f_coll = N_min_cell * M_MIN / (pixel_mass*density_over_mean);
                 }
               }
               else{
                 f_coll = *((float *)M_coll_unfiltered + HII_R_FFT_INDEX(x,y,z)) / (pixel_mass*density_over_mean);
               }
               if (f_coll>1) f_coll=1;
               res_xH = xHI_from_xrays - f_coll * ION_EFF_FACTOR;
               // and make sure fraction doesn't blow up for underdense pixels
               if (res_xH < 0)
                 res_xH = 0;
               else if (res_xH > 1)
                 res_xH = 1;

               xH[HII_R_INDEX(x, y, z)] = res_xH;
             } // end partial ionizations at last filtering step

        } // k
      } // j
    } // i

    R /= DELTA_R_HII_FACTOR;
  }
  fprintf(stderr, "find_HII_bubbles:: finished main loop scroll, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);

  // find the neutral fraction
  global_xH = 0;
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
    global_xH += xH[ct];
  }
  global_xH /= (float)HII_TOT_NUM_PIXELS;

  // save the xH box
  global_xH = write_xH(meraxes_fname, snapshot, xH, global_xH, ION_EFF_FACTOR, M_MIN, MFP);

  // deallocate
  fftwf_cleanup_threads();
  gsl_rng_free (r);
  fftwf_free(xH);
  fclose(LOG);
  fftwf_free(deltax_unfiltered);
  fftwf_free(deltax_filtered);
  if (USE_HALO_FIELD){  fftwf_free(M_coll_unfiltered);  fftwf_free(M_coll_filtered);}
  free_ps(); if (USE_TS_IN_21CM){ fftwf_free(xe_filtered); fftwf_free(xe_unfiltered);} return (int) (global_xH * 100);
  }
