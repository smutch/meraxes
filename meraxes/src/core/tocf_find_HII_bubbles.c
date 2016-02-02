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
  int HII_dim = tocf_params.HII_dim;
  int HII_middle = HII_dim / 2;
  float box_size = run_globals.params.BoxSize;
  float delta_k = M_PI / box_size;
  float k_x, k_y, k_z, k_mag, kR;

  // Loop through k-box
  for (int n_x=0; n_x<HII_dim; n_x++)
  {
    if (n_x>HII_middle)
      k_x = (n_x-HII_dim)*delta_k;
    else
      k_x = n_x*delta_k;

    for (int n_y=0; n_y<HII_dim; n_y++)
    {
      if (n_y>HII_middle)
        k_y = (n_y-HII_dim)*delta_k;
      else
        k_y = n_y*delta_k;

      for (int n_z=0; n_z<=HII_middle; n_z++)
      { 
        k_z = n_z*delta_k;

        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

        kR = k_mag*R;   // Real space top-hat

        switch(filter_type)
        {
          case 0:   // Real space top-hat
            if (kR > 1e-4)
              box[HII_C_INDEX(n_x, n_y, n_z)] *= 3.0 * (sinf(kR)/powf(kR, 3) - cosf(kR)/powf(kR, 2));
            break;

          case 1:   // k-space top hat
            kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
            if (kR > 1)
              box[HII_C_INDEX(n_x, n_y, n_z)] = 0.0;
            break;
        
          case 2:   // Gaussian
            kR *= 0.643;   // Equates integrated volume to the real space top-hat
            box[HII_C_INDEX(n_x, n_y, n_z)] *= pow(M_E, -kR*kR/2.0);
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
  float pixel_mass = RtoM(L_FACTOR*box_size/(float)HII_dim);  // 1e10 Msol/h
  float cell_length_factor = L_FACTOR;
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
  int slab_n_real = tocf_params.slab_nix[SID.My_rank] * HII_dim * HII_dim;
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
  int slab_n_complex = tocf_params.slab_n_complex[SID.My_rank];
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

  }

}

#endif
