#include "meraxes.h"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>

// DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et 
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger 
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were 
 * subsequently made by Simon Mutch & Paul Geil.
 */

double RtoM(double R){
  // All in internal units
  int filter = run_globals.params.ReionRtoMFilterType;
  double OmegaM = run_globals.params.OmegaM;
  double RhoCrit = run_globals.RhoCrit;

  switch (filter)
  {
   case 0: //top hat M = (4/3) PI <rho> R^3
    return (4.0/3.0)*PI*pow(R,3)*(OmegaM*RhoCrit);
    break;
   case 1: //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow(2*PI, 1.5) * OmegaM*RhoCrit * pow(R, 3);
    break;
   default: // filter not defined
    SID_log_error("Unrecognised filter (%d). Aborting...", filter);
    ABORT(EXIT_FAILURE);
    break;
  }

  return -1;
}


static void filter(fftwf_complex *box, float R)
{
  int filter_type = run_globals.params.ReionFilterType;
  int slab_nx = (int)(run_globals.reion_grids.slab_nix[SID.My_rank]);
  int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[SID.My_rank]);
  int ReionGridDim = run_globals.params.ReionGridDim;
  int HII_middle = ReionGridDim / 2;
  float box_size = run_globals.params.BoxSize;
  float delta_k = M_PI / box_size;

  // Loop through k-box
  for (int n_x=0; n_x<slab_nx; n_x++)
  {
    float k_x;
    int n_x_global = n_x + local_ix_start;

    if (n_x_global > HII_middle)
      k_x = (n_x_global - ReionGridDim)*delta_k;
    else
      k_x = n_x_global*delta_k;

    for (int n_y=0; n_y<ReionGridDim; n_y++)
    {
      float k_y;

      if (n_y>HII_middle)
        k_y = (n_y-ReionGridDim)*delta_k;
      else
        k_y = n_y*delta_k;

      for (int n_z=0; n_z<=HII_middle; n_z++)
      { 
        float k_z = n_z*delta_k;

        float k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

        float kR = k_mag*R;   // Real space top-hat

        switch(filter_type)
        {
          case 0:   // Real space top-hat
            if (kR > 1e-4)
              box[grid_index(n_x, n_y, n_z, ReionGridDim, INDEX_COMPLEX_HERM)] *= 3.0 * (sinf(kR)/powf(kR, 3) - cosf(kR)/powf(kR, 2));
            break;

          case 1:   // k-space top hat
            kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
            if (kR > 1)
              box[grid_index(n_x, n_y, n_z, ReionGridDim, INDEX_COMPLEX_HERM)] = 0.0;
            break;
        
          case 2:   // Gaussian
            kR *= 0.643;   // Equates integrated volume to the real space top-hat
            box[grid_index(n_x, n_y, n_z, ReionGridDim, INDEX_COMPLEX_HERM)] *= powf(M_E, -kR*kR/2.0);
            break;

          default:
            if ( (n_x==0) && (n_y==0) && (n_z==0) )
            {
              SID_log_error("ReionFilterType.c: Warning, ReionFilterType type %d is undefined!", filter_type);
              ABORT(EXIT_FAILURE);
            }
            break;
        }

      }

    }
  }   // End looping through k box

}


float find_HII_bubbles(float redshift)
{
  // TODO: TAKE A VERY VERY CLOSE LOOK AT UNITS!!!!

  float box_size = run_globals.params.BoxSize;  // Mpc/h
  int ReionGridDim = run_globals.params.ReionGridDim;
  float pixel_volume = powf(box_size/(float)ReionGridDim, 3);  // (Mpc/h)^3
  float l_factor = 0.620350491;  // Factor relating cube length to filter radius = (4PI/3)^(-1/3)
  float cell_length_factor = l_factor;
  int total_n_cells = (int)pow(ReionGridDim, 3);

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  if ((box_size/(float)ReionGridDim) < 1.0)   // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0; 

  // Init J_21
  int flag_ReionUVBFlag = run_globals.params.ReionUVBFlag;
  int slab_n_real = (int)(run_globals.reion_grids.slab_nix[SID.My_rank]) * ReionGridDim * ReionGridDim;
  float *J_21 = run_globals.reion_grids.J_21;
  if (flag_ReionUVBFlag)
    for(int ii=0; ii < slab_n_real; ii++)
      J_21[ii] = 0.0;

  // Init xH
  float *xH = run_globals.reion_grids.xH;
  for(int ii=0; ii < slab_n_real; ii++)
    xH[ii] = 1.0;

  // Forward fourier transform to obtain k-space fields
  // TODO: Ensure that fftwf_mpi_init has been called and fftwf_mpi_cleanup will be called
  // TODO: Don't use estimate and calculate plan in code init
  float *deltax = run_globals.reion_grids.deltax;
  fftwf_complex *deltax_unfiltered = (fftwf_complex *)deltax;  // WATCH OUT!
  fftwf_complex *deltax_filtered = run_globals.reion_grids.deltax_filtered;
  fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax, deltax_unfiltered, SID_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  float *stars = run_globals.reion_grids.stars;
  fftwf_complex *stars_unfiltered = (fftwf_complex *)stars;  // WATCH OUT!
  fftwf_complex *stars_filtered = run_globals.reion_grids.stars_filtered;
  plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars, stars_unfiltered, SID_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  
  float *sfr = run_globals.reion_grids.sfr;
  fftwf_complex *sfr_unfiltered = (fftwf_complex *)sfr;  // WATCH OUT!
  fftwf_complex *sfr_filtered = run_globals.reion_grids.sfr_filtered;
  plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr, sfr_unfiltered, SID_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  
  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  // TODO: Double check that looping over correct number of elements here
  int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[SID.My_rank]);
  for (int ii=0; ii<slab_n_complex; ii++)
  {
    deltax_unfiltered[ii] /= (float)total_n_cells;
    stars_unfiltered[ii] /= (float)total_n_cells;
    sfr_unfiltered[ii] /= (float)total_n_cells;
  }

  // Loop through filter radii
  float ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMax;  // Mpc/h
  float ReionRBubbleMin = run_globals.params.physics.ReionRBubbleMin;  // Mpc/h
  float R = fmin(ReionRBubbleMax, l_factor * box_size);  // Mpc/h
  float ReionDeltaRFactor = run_globals.params.ReionDeltaRFactor;
  float ReionGammaHaloBias = run_globals.params.physics.ReionGammaHaloBias;
  
  bool flag_last_filter_step = false;

  while(!flag_last_filter_step)
  {
    // check to see if this is our last filtering step
    if( ((R/ReionDeltaRFactor) <= (cell_length_factor*box_size/(float)ReionGridDim))
        || ((R/ReionDeltaRFactor) <= ReionRBubbleMin) )
    {
      flag_last_filter_step = true;
      R = cell_length_factor * box_size / (float)ReionGridDim;
    }

    // copy the k-space grids
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*slab_n_complex);
    memcpy(stars_filtered, stars_unfiltered, sizeof(fftwf_complex)*slab_n_complex);
    memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex)*slab_n_complex);

    // do the filtering unless this is the last filter step
    if(!flag_last_filter_step)
    {
      filter(deltax_filtered, R);
      filter(stars_filtered, R);
      filter(sfr_filtered, R);
    }

    // inverse fourier transform back to real space
    plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax_filtered, (float *)deltax_filtered, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars_filtered, (float *)stars_filtered, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr_filtered, (float *)sfr_filtered, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Perform sanity checks to account for aliasing effects
    int local_nix = (int)(run_globals.reion_grids.slab_nix[SID.My_rank]);
    for (int ix=0; ix<local_nix; ix++)
      for (int iy=0; iy<ReionGridDim; iy++)
        for (int iz=0; iz<ReionGridDim; iz++)
        {   
          ((float *)deltax_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = fmaxf(((float *)deltax_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)], -1 + REL_TOL);

          ((float *)stars_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = fmaxf(((float *)stars_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)], 0.0);

          ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = fmaxf(((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] , 0.0);
        }

    /*
     * Main loop through the box...
     */

    float ReionEfficiency = run_globals.params.physics.ReionEfficiency;
    float ReionNionPhotPerBary = run_globals.params.physics.ReionNionPhotPerBary;
    run_units_t *units = &(run_globals.units);

    float J_21_aux_constant = 1e21 * (1.0+redshift)*(1.0+redshift)/(4.0 * M_PI)
      * run_globals.params.physics.ReionAlphaUV * PLANCK * run_globals.params.physics.ReionEscapeFrac
      * R * units->UnitLength_in_cm
      * ReionNionPhotPerBary / PROTONMASS
      * units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3) / units->UnitTime_in_s;
  
    // DEBUG
    // for (int ix=0; ix<local_nix; ix++)
    //   for (int iy=0; iy<ReionGridDim; iy++)
    //     for (int iz=0; iz<ReionGridDim; iz++)
    //       if (((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] > 0)
    //       {
    //         SID_log("J_21_aux ==========", SID_LOG_OPEN);
    //         SID_log("redshift = %.2f", SID_LOG_COMMENT, redshift);
    //         SID_log("alpha = %.2e", SID_LOG_COMMENT, run_globals.params.physics.ReionAlphaUV);
    //         SID_log("PLANCK = %.2e", SID_LOG_COMMENT, PLANCK);
    //         SID_log("ReionEscapeFrac = %.2f", SID_LOG_COMMENT, ReionEscapeFrac);
    //         SID_log("ReionNionPhotPerBary = %.1f", SID_LOG_COMMENT, ReionNionPhotPerBary);
    //         SID_log("UnitMass_in_g = %.2e", SID_LOG_COMMENT, units->UnitMass_in_g);
    //         SID_log("UnitLength_in_cm = %.2e", SID_LOG_COMMENT, units->UnitLength_in_cm);
    //         SID_log("PROTONMASS = %.2e", SID_LOG_COMMENT, PROTONMASS);
    //         SID_log("-> J_21_aux_constant = %.2e", SID_LOG_COMMENT, J_21_aux_constant);
    //         SID_log("pixel_volume = %.2e", SID_LOG_COMMENT, pixel_volume);
    //         SID_log("sfr_density[%d, %d, %d] = %.2e", SID_LOG_COMMENT, ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume);
    //         SID_log("-> J_21_aux = %.2e", SID_LOG_COMMENT, ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume * J_21_aux_constant);
    //         SID_log("==========", SID_LOG_CLOSE);
    //         if (SID.My_rank == 0)
    //           ABORT(EXIT_SUCCESS);
    //       }

    for (int ix=0; ix<local_nix; ix++)
      for (int iy=0; iy<ReionGridDim; iy++)
        for (int iz=0; iz<ReionGridDim; iz++)
        {
          float density_over_mean = 1.0 + ((float *)deltax_filtered)[grid_index(ix,iy,iz, ReionGridDim, INDEX_PADDED)];

          float f_coll_stars =  ((float *)stars_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]/ (RtoM(R)*density_over_mean);
          f_coll_stars *= (4.0/3.0)*PI*pow(R,3.0) / pixel_volume;

// #ifdef DEBUG
//           debug("%d, %g, %g, %g, %g, %g, %g, %g, %g\n", SID.My_rank,
//               R, density_over_mean, f_coll_stars,
//               ((float *)stars_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)],
//               RtoM(R), (4.0/3.0)*PI*pow(R,3.0), pixel_volume, 1.0/ReionEfficiency);
// #endif

          float sfr_density = ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume;   // In internal units

          // TODO: I fixed an incorrect factor in this equation (1-Y_He instead
          // of 1-0.75*Y_He) that will need to be reverted when comparing with
          // old results. 
          float J_21_aux;
          if (flag_ReionUVBFlag)
            J_21_aux = sfr_density * J_21_aux_constant;

          // Check if ionised!
          if (f_coll_stars > 1.0/ReionEfficiency)   // IONISED!!!!
          {   
            // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
            if (xH[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] > REL_TOL)
              if(flag_ReionUVBFlag)
                J_21[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] = J_21_aux;


            // Mark as ionised
            xH[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] = 0;

          }
          // Check if this is the last filtering step.
          // If so, assign partial ionisations to those cells which aren't fully ionised 
          else if (flag_last_filter_step && (xH[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] > REL_TOL))
          {
            float res_xH = 1.0 - f_coll_stars * ReionEfficiency;
            xH[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] = res_xH;
          }

          // Check if new ionisation
          float *z_in = run_globals.reion_grids.z_at_ionization;
          if ( (xH[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] < REL_TOL) && (z_in[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] < 0) )   // New ionisation!
          {
            z_in[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] = redshift;
            if (flag_ReionUVBFlag)
              run_globals.reion_grids.J_21_at_ionization[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)] = J_21_aux * ReionGammaHaloBias;
          }
        } // iz

    R /= ReionDeltaRFactor;

  }


  // Find the neutral fraction
  float global_xH = 0.0;
  for (int ct=0; ct < slab_n_real; ct++)
    global_xH += xH[ct];
  SID_Allreduce(SID_IN_PLACE, &global_xH, 1, SID_FLOAT, SID_SUM, SID.COMM_WORLD);
  global_xH /= (float)pow(ReionGridDim, 3);

  return global_xH;
}
