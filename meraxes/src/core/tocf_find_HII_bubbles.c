#ifdef USE_TOCF

#include "meraxes.h"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et 
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger 
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were 
 * subsequently made by Simon Mutch & Paul Geil.
 */

// R in Mpc/h, M in 1e10 Msun/h 
double RtoM(double R){
  int filter = tocf_params.RtoM_filter;
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


void HII_filter(fftwf_complex *box, float R)
{
  int filter_type = tocf_params.HII_filter;
  int slab_nx = (int)(tocf_params.slab_nix[SID.My_rank]);
  int local_ix_start = (int)(tocf_params.slab_ix_start[SID.My_rank]);
  int HII_dim = tocf_params.HII_dim;
  int HII_middle = HII_dim / 2;
  float box_size = run_globals.params.BoxSize;
  float delta_k = M_PI / box_size;

  // Loop through k-box
  for (int n_x=0; n_x<slab_nx; n_x++)
  {
    float k_x;
    int n_x_global = n_x + local_ix_start;

    if (n_x_global > HII_middle)
      k_x = (n_x_global - HII_dim)*delta_k;
    else
      k_x = n_x_global*delta_k;

    for (int n_y=0; n_y<HII_dim; n_y++)
    {
      float k_y;

      if (n_y>HII_middle)
        k_y = (n_y-HII_dim)*delta_k;
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
              box[grid_index(n_x, n_y, n_z, HII_dim, INDEX_COMPLEX_HERM)] *= 3.0 * (sinf(kR)/powf(kR, 3) - cosf(kR)/powf(kR, 2));
            break;

          case 1:   // k-space top hat
            kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
            if (kR > 1)
              box[grid_index(n_x, n_y, n_z, HII_dim, INDEX_COMPLEX_HERM)] = 0.0;
            break;
        
          case 2:   // Gaussian
            kR *= 0.643;   // Equates integrated volume to the real space top-hat
            box[grid_index(n_x, n_y, n_z, HII_dim, INDEX_COMPLEX_HERM)] *= powf(M_E, -kR*kR/2.0);
            break;

          default:
            if ( (n_x==0) && (n_y==0) && (n_z==0) )
            {
              SID_log_error("HII_filter.c: Warning, HII_filter type %d is undefined!", filter_type);
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
  int HII_dim = tocf_params.HII_dim;
  float pixel_volume = powf(box_size/(float)HII_dim, 3);  // (Mpc/h)^3
  float l_factor = 0.620350491;  // Factor relating cube length to filter radius = (4PI/3)^(-1/3)
  float cell_length_factor = l_factor;
  double ion_tvir_min = tocf_params.ion_tvir_min;
  int total_n_cells = (int)pow(HII_dim, 3);

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  if ((box_size/(float)HII_dim) < 1.0)   // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0; 

  double M_min = 0.0;  // Minimum source mass  (1e10 Msol/h)
  if(tocf_params.ion_tvir_min > 0.0)
    M_min = Tvir_to_Mvir(ion_tvir_min, redshift);

  // Init J_21
  int flag_uvb_feedback = tocf_params.uvb_feedback;
  int slab_n_real = (int)(tocf_params.slab_nix[SID.My_rank]) * HII_dim * HII_dim;
  float *J_21 = run_globals.tocf_grids.J_21;
  for(int ii=0; ii < slab_n_real; ii++)
    J_21[ii] = 0.0;

  // Init xH
  float *xH = run_globals.tocf_grids.xH;
  for(int ii=0; ii < slab_n_real; ii++)
    xH[ii] = 0.0;

  // Forward fourier transform to obtain k-space fields
  // TODO: Ensure that fftwf_mpi_init has been called and fftwf_mpi_cleanup will be called
  float *deltax = run_globals.tocf_grids.deltax;
  fftwf_complex *deltax_unfiltered = (fftwf_complex *)deltax;  // WATCH OUT!
  fftwf_complex *deltax_filtered = run_globals.tocf_grids.deltax_filtered;
  fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(HII_dim, HII_dim, HII_dim, deltax, deltax_unfiltered, SID_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  float *stars = run_globals.tocf_grids.stars;
  fftwf_complex *stars_unfiltered = (fftwf_complex *)stars;  // WATCH OUT!
  fftwf_complex *stars_filtered = run_globals.tocf_grids.stars_filtered;
  plan = fftwf_mpi_plan_dft_r2c_3d(HII_dim, HII_dim, HII_dim, stars, stars_unfiltered, SID_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  
  float *sfr = run_globals.tocf_grids.sfr;
  fftwf_complex *sfr_unfiltered = (fftwf_complex *)sfr;  // WATCH OUT!
  fftwf_complex *sfr_filtered = run_globals.tocf_grids.sfr_filtered;
  plan = fftwf_mpi_plan_dft_r2c_3d(HII_dim, HII_dim, HII_dim, sfr, sfr_unfiltered, SID_COMM_WORLD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  
  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  // TODO: Double check that looping over correct number of elements here
  int slab_n_complex = (int)(tocf_params.slab_n_complex[SID.My_rank]);
  for (int ii=0; ii<slab_n_complex; ii++)
  {
    deltax_unfiltered[ii] /= (float)total_n_cells;
    stars_unfiltered[ii] /= (float)total_n_cells;
    sfr_unfiltered[ii] /= (float)total_n_cells;
  }

  // Loop through filter radii
  float r_bubble_max = tocf_params.r_bubble_max;  // Mpc/h
  float r_bubble_min = tocf_params.r_bubble_min;  // Mpc/h
  float R = fmin(r_bubble_max, l_factor * box_size);  // Mpc/h
  float delta_r_HII_factor = tocf_params.delta_r_HII_factor;
  float gamma_halo_bias = tocf_params.gamma_halo_bias;
  
  bool flag_last_filter_step = false;

  while(!flag_last_filter_step)
  {
    // check to see if this is our last filtering step
    if( ((R/delta_r_HII_factor) <= (cell_length_factor*box_size/(float)HII_dim))
        || ((R/delta_r_HII_factor) <= r_bubble_min) )
    {
      flag_last_filter_step = true;
      R = cell_length_factor * box_size / (float)HII_dim;
    }

    // copy the k-space grids
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*slab_n_complex);
    memcpy(stars_filtered, stars_unfiltered, sizeof(fftwf_complex)*slab_n_complex);
    memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex)*slab_n_complex);

    // do the filtering unless this is the last filter step
    if(!flag_last_filter_step)
    {
      HII_filter(deltax_filtered, R);
      HII_filter(stars_filtered, R);
      HII_filter(sfr_filtered, R);
    }

    // inverse fourier transform back to real space
    plan = fftwf_mpi_plan_dft_c2r_3d(HII_dim, HII_dim, HII_dim, deltax_filtered, (float *)deltax_filtered, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    plan = fftwf_mpi_plan_dft_c2r_3d(HII_dim, HII_dim, HII_dim, stars_filtered, (float *)stars_filtered, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    plan = fftwf_mpi_plan_dft_c2r_3d(HII_dim, HII_dim, HII_dim, sfr_filtered, (float *)sfr_filtered, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Perform sanity checks to account for aliasing effects
    int local_nix = (int)(tocf_params.slab_nix[SID.My_rank]);
    for (int ix=0; ix<local_nix; ix++)
      for (int iy=0; iy<HII_dim; iy++)
        for (int iz=0; iz<HII_dim; iz++)
        {   
          ((float *)deltax_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_PADDED)] = fmaxf(((float *)deltax_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)], -1 + REL_TOL);

          ((float *)stars_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_PADDED)] = fmaxf(((float *)stars_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)], 0.0);

          ((float *)sfr_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_PADDED)] = fmaxf(((float *)sfr_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] , 0.0);
        }

    /*
     * Main loop through the box...
     */

    int flag_uvb_feedback = tocf_params.uvb_feedback;
    float BaryonFrac = run_globals.params.BaryonFrac;
    float HII_eff_factor = tocf_params.HII_eff_factor;
    for (int ix=0; ix<HII_dim; ix++)
    {   
      for (int iy=0; iy<HII_dim; iy++)
      {
        for (int iz=0; iz<HII_dim; iz++)
        {
          float density_over_mean = 1.0 + ((float *)deltax_filtered)[grid_index(ix,iy,iz, HII_dim, INDEX_PADDED)];

          float f_coll_stars =  ((float *)stars_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_PADDED)]/ (RtoM(R)*density_over_mean);
          f_coll_stars *= (4.0/3.0)*PI*pow(R,3.0) / pixel_volume;

          float sfr_density = ((float *)sfr_filtered)[grid_index(ix, iy, iz, HII_dim, INDEX_PADDED)] / pixel_volume;   // In units of Msolar/s/cMpc/cMpc/cMpc

          // Adjust the denominator of the collapse fraction for the residual electron fraction in the neutral medium
          // Calculate the mfp of the ionising photons for this size 
          float J_21_aux;
          if (flag_uvb_feedback)
            J_21_aux = (1.0/(4.0 * M_PI)) * 1e21 * tocf_params.alpha_uv * PLANCK * ((1.0+redshift)*(1.0+redshift))*(R*MPC)
              * HII_eff_factor* BaryonFrac *(1.0-tocf_params.Y_He)
              * sfr_density*SOLAR_MASS/MPC/MPC/MPC / PROTONMASS;

          // Check if ionised!
          if (f_coll_stars > 1.0/HII_eff_factor)   // IONISED!!!!
          {   
            // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
            if (xH[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] > REL_TOL)
              if(flag_uvb_feedback)
                J_21[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] = J_21_aux*gamma_halo_bias;


            // Mark as ionised
            xH[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] = 0;

          }
          // Check if this is the last filtering step.
          // If so, assign partial ionisations to those cells which aren't fully ionised 
          else if (flag_last_filter_step && (xH[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] > REL_TOL))
          {
            float res_xH = 1.0 - f_coll_stars * HII_eff_factor;
            xH[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] = res_xH;
          }

          // Check if new ionisation
          float *z_in = run_globals.tocf_grids.z_at_ionization;
          if ( (xH[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] < REL_TOL) && (z_in[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] < 0) )   // New ionisation!
          {
            z_in[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] = redshift;
            if (flag_uvb_feedback)
              run_globals.tocf_grids.J_21_at_ionization[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)] = J_21_aux * gamma_halo_bias;
          }

        } // iz
      } // iy
    } // ix

    R /= delta_r_HII_factor;

  }


  // Find the neutral fraction
  float global_xH = 0.0;
  for (int ct=0; ct < slab_n_real; ct++)
    global_xH += xH[ct];
  SID_Allreduce(SID_IN_PLACE, &global_xH, 1, SID_FLOAT, SID_SUM, SID.COMM_WORLD);
  global_xH /= (float)pow(HII_dim, 3);

  // Renormalise the J_21 box
  if(flag_uvb_feedback)
    for(int ct=0; ct < slab_n_real; ct++)
      J_21[ct] /= gamma_halo_bias;

  // cleanup
  fftwf_mpi_cleanup();

  return global_xH;
}

#endif
