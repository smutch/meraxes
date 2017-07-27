#include "meraxes.h"
#include "meraxes_gpu.h"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>

#include <hdf5.h>
#include <hdf5_hl.h>

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
    fftwf_complex *deltax_filtered,  // complex
    fftwf_complex *stars_filtered,  // complex
    fftwf_complex *sfr_filtered,  // complex

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
  double       pixel_volume         = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
  double       cell_length_factor   = L_FACTOR;
  double       total_n_cells        = pow((double)ReionGridDim, 3);
  int          slab_n_real          = local_nix * ReionGridDim * ReionGridDim;
  float        J_21_aux;
  double       density_over_mean;
  double       sfr_density;
  double       f_coll_stars;
  int          i_real;
  int          i_padded;

  int slab_n_complex = (int)(slabs_n_complex[mpi_rank]);

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
  if (flag_ReionUVBFlag)
    for(int ii = 0; ii < slab_n_real; ii++)
      J_21[ii] = 0.0;

  // Init xH
  for(int ii = 0; ii < slab_n_real; ii++)
    xH[ii] = 1.0;

  // Init r_bubble
  for(int ii = 0; ii < slab_n_real; ii++)
    r_bubble[ii] = 0.0;

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
  for (int ii = 0; ii < slab_n_complex; ii++)
  {
    deltax_unfiltered[ii] /= total_n_cells;
    stars_unfiltered[ii]  /= total_n_cells;
    sfr_unfiltered[ii]    /= total_n_cells;
  }

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
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
    memcpy(stars_filtered, stars_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
    memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex);

    // do the filtering unless this is the last filter step
    int local_ix_start = (int)(slabs_ix_start[mpi_rank]);
    if(!flag_last_filter_step)
    {
      filter((fftwf_complex *)deltax_filtered, local_ix_start, local_nix, ReionGridDim, (float)R);
      filter((fftwf_complex *)stars_filtered,  local_ix_start, local_nix, ReionGridDim, (float)R);
      filter((fftwf_complex *)sfr_filtered,    local_ix_start, local_nix, ReionGridDim, (float)R);
    }

    // inverse fourier transform back to real space
    plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)stars_filtered, (float *)stars_filtered, mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, (fftwf_complex *)sfr_filtered, (float *)sfr_filtered, mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Perform sanity checks to account for aliasing effects
    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++)
        {
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
          ((float *)deltax_filtered)[i_padded] = fmaxf(((float *)deltax_filtered)[i_padded], -1 + REL_TOL);
          ((float *)stars_filtered)[i_padded]  = fmaxf(((float *)stars_filtered)[i_padded], 0.0);
          ((float *)sfr_filtered)[i_padded]    = fmaxf(((float *)sfr_filtered)[i_padded], 0.0);
        }

    /*
     * Main loop through the box...
     */

    double J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI)
      * ReionAlphaUV * PLANCK
      * 1e21 * ReionEscapeFrac
      * R *UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS
      * UnitMass_in_g / pow(UnitLength_in_cm, 3) / UnitTime_in_s;

    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++)
        {
          i_real   = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

          density_over_mean = 1.0 + (double)((float *)deltax_filtered)[i_padded];

          f_coll_stars      =  (double)((float *)stars_filtered)[i_padded] / (RtoM(R) * density_over_mean)
                               * (4.0 / 3.0) * M_PI * pow(R,3.0) / pixel_volume;

          sfr_density       = (double)((float *)sfr_filtered)[i_padded] / pixel_volume; // In internal units

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
//fprintf(stderr,"ix,iy,iz=%03d,%03d,%03d R=%le\n",ix,iy,iz,r_bubble[i_real]);
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
//if(iz<100) fprintf(stderr,"%le ",J_21_at_ionization[i_real]);if(iz==100) fprintf(stderr," XXXXXXX\n");
        }
    // iz

    R /= ReionDeltaRFactor;
  }


  // Find the volume and mass weighted neutral fractions
  // TODO: The deltax grid will have rounding errors from forward and reverse
  //       FFT. Should cache deltax slabs prior to ffts and reuse here.
  *volume_weighted_global_xH = 0.0;
  *mass_weighted_global_xH   = 0.0;
  double mass_weight         = 0.0;

  for (int ix = 0; ix < local_nix; ix++)
    for (int iy = 0; iy < ReionGridDim; iy++)
      for (int iz = 0; iz < ReionGridDim; iz++)
      {
        i_real   = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
        i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
        *volume_weighted_global_xH += (double)xH[i_real];
        density_over_mean           = 1.0 + (double)((float *)deltax_filtered)[i_padded];
        *mass_weighted_global_xH   += (double)(xH[i_real]) * density_over_mean;
        mass_weight                += density_over_mean;
      }

  MPI_Allreduce(MPI_IN_PLACE, &volume_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &mass_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &mass_weight, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

  *volume_weighted_global_xH                        /= total_n_cells;
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
        fftwf_complex *deltax_filtered,  // complex
        fftwf_complex *stars_filtered,  // complex
        fftwf_complex *sfr_filtered,  // complex
    
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
    float *deltax = fftwf_alloc_real(slab_n_complex * 2);
    float *stars  = fftwf_alloc_real(slab_n_complex * 2);
    float *sfr    = fftwf_alloc_real(slab_n_complex * 2);
    H5LTread_dataset_float(file_id, "deltax", deltax);
    H5LTread_dataset_float(file_id, "sfr",  sfr);
    H5LTread_dataset_float(file_id, "stars", stars);
    H5Fclose(file_id);

    // Initialize outputs    
    float         *J_21                     =fftwf_alloc_real(slab_n_real);
    float         *r_bubble                 =fftwf_alloc_real(slab_n_real); 
    fftwf_complex *deltax_filtered          =fftwf_alloc_complex(slab_n_complex);
    fftwf_complex *stars_filtered           =fftwf_alloc_complex(slab_n_complex);  
    fftwf_complex *sfr_filtered             =fftwf_alloc_complex(slab_n_complex);
    float         *xH                       =fftwf_alloc_real(slab_n_real);
    float         *z_at_ionization          =fftwf_alloc_real(slab_n_real);
    float         *J_21_at_ionization       =fftwf_alloc_real(slab_n_real);
    double         volume_weighted_global_xH;
    double         mass_weighted_global_xH;
    
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
    fftwf_free(J_21);
    fftwf_free(r_bubble);
    fftwf_free(deltax);
    fftwf_free(stars);
    fftwf_free(sfr);
    fftwf_free(deltax_filtered);
    fftwf_free(stars_filtered);
    fftwf_free(sfr_filtered);
    fftwf_free(xH);
    fftwf_free(z_at_ionization);
    fftwf_free(J_21_at_ionization);

}
