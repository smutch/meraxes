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

int delta_T_ps(
    int            snapshot,
    float         *average_T,
    float         *delta_T,
    float         **ps,
    int           *ps_nbins)
{
  fftwf_complex *deldel_T;
  fftwf_plan plan;

  bool flag_use_vel = false;
  bool flag_use_ts = false;

  unsigned long long ct, temp_ct;

  float pixel_x_HI, pixel_deltax;
  float maxi, maxj, maxk, maxdvdx, mini, minj, mink, mindvdx;
  float k_mag;
  float T_rad, pixel_Ts_factor;

  double dvdx, max_v_deriv;
  double ave_Ts, min_Ts, max_Ts, temp;

  float *xH = run_globals.tocf_grids.xH;
  float *deltax = run_globals.tocf_grids.deltax;

  // Begin initialisation
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  float min       = 1e3;
  float max       = -1e3;
  double ave       = 0.0;
  unsigned long long nonlin_ct = 0;
  int HII_dim = tocf_params.HII_dim;
  int local_nix = (int)(run_globals.tocf_grids.slab_nix[SID.My_rank]);

  // Set some redshift dependant values
  float redshift = run_globals.ZZ[snapshot];
  T_rad = TCMB*(1 + redshift);
  float H = (float)hubble_at_snapshot(snapshot);
  float Hubble_h = run_globals.params.Hubble_h;
  float OmegaM = run_globals.params.OmegaM;
  float OmegaB = OmegaM * run_globals.params.BaryonFrac;
  float const_factor = 27.0*(OmegaB * Hubble_h * Hubble_h / 0.023)*sqrt((0.15/OmegaM/Hubble_h/Hubble_h)*(1.0+redshift)/10.0);


  // delta_T grid (same size as the bubble box)
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  ave_Ts = 0.0;
  min_Ts = 1e5;
  max_Ts = 0.0;
  temp = 0.0;
  temp_ct = 0;

  for (int ii=0; ii<local_nix; ii++)
  {
    for (int jj=0; jj<HII_dim; jj++)
    {
      for (int kk=0; kk<HII_dim; kk++)
      {
        pixel_deltax = deltax[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)];

        int index = grid_index(ii, jj, kk, HII_dim, INDEX_REAL);
        pixel_x_HI = xH[index];

        if (pixel_x_HI > ABS_TOL)
        {
          temp = pixel_deltax;
          temp_ct++;
        }

        delta_T[index] = const_factor*pixel_x_HI*(1.0 + pixel_deltax);

        if (max < delta_T[index])
          max = delta_T[index];

        if (min > delta_T[index])
          min = delta_T[index];

        ave += delta_T[index];
      }
    }
  }

  int tot_num_pixels = (int)pow(HII_dim, 3);

  SID_Allreduce(SID_IN_PLACE, &ave, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  ave /= (double)tot_num_pixels;


  // Calculate power spectrum
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  float k_factor = 1.5;
  float delta_k = tocf_params.delta_k_ps;
  float k_first_bin_ceil = delta_k;
  float k_max = delta_k*HII_dim;

  // Initialise arrays
  // ghetto counting (lookup how to do logs of arbitrary bases in c...)
  int num_bins = 0;
  float k_floor = 0.0;
  float k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max)
  {
    num_bins++;
    k_floor = k_ceil;
    k_ceil *= k_factor;
  }

  //    printf("    NUM_BINS       = %d\n", NUM_BINS);

  double *p_box     = malloc(sizeof(double)*num_bins);
  double *k_ave     = malloc(sizeof(double)*num_bins);
  unsigned long long *in_bin_ct = malloc(sizeof(unsigned long long)*num_bins);

  for (int ii=0; ii < num_bins; ii++)
  {
    p_box[ii] = 0.0;
    k_ave[ii] = 0.0;
    in_bin_ct[ii] = 0;
  }

  deldel_T = fftwf_alloc_complex(run_globals.tocf_grids.slab_n_complex[SID.My_rank]);

  if (!deldel_T)
  {
    fprintf(stderr, "Unable to allocate memory for the deldel_T box!\n");
    fprintf(LOG, "Unable to allocate memory for the deldel_T box!\n");
    fclose(LOG);
    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    fftwf_cleanup_threads();
    return -1;
  }

  // Fill-up the real-space of the deldel box
  for (i=0; i<HII_DIM; i++)
  {
    for (j=0; j<HII_DIM; j++)
    {
      for (k=0; k<HII_DIM; k++)
      {
        *((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) = (delta_T[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(float)HII_TOT_NUM_PIXELS;
        // Note: we include the V/N factor for the scaling after the fft
      }
    }
  }

  // Transform to k-space
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();

  // Now construct the power spectrum
  for (n_x=0; n_x<HII_DIM; n_x++)
  {
    if (n_x>HII_MIDDLE)
    {
      k_x =(n_x-HII_DIM)*DELTA_K;   // Wrap around for FFT convention
    }
    else
    {
      k_x = n_x*DELTA_K;
    }

    for (n_y=0; n_y<HII_DIM; n_y++)
    {
      if (n_y>HII_MIDDLE)
      {
        k_y =(n_y-HII_DIM)*DELTA_K;
      }
      else
      {
        k_y = n_y*DELTA_K;
      }

      for (n_z=0; n_z<=HII_MIDDLE; n_z++)
      { 
        k_z = n_z*DELTA_K;

        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

        // Now go through the k bins and update
        ct = 0;
        k_floor = 0.0;
        k_ceil = k_first_bin_ceil;
        while (k_ceil < k_max)
        {
          // Check if we fal in this bin
          if ((k_mag >= k_floor) && (k_mag < k_ceil))
          {
            in_bin_ct[ct]++;
            p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
            // Note the 1/VOLUME factor, which turns this into a power density in k-space

            k_ave[ct] += k_mag;
            break;
          }

          ct++;
          k_floor = k_ceil;
          k_ceil *= k_factor;
        }
      }
    }
  }   // End looping through k box


  // Malloc and store the result
  *ps = (float *)calloc(3*NUM_BINS, sizeof(float));

  // NOTE - previous ct ran from 1 (not zero) to NUM_BINS
  for (ct=0; ct<NUM_BINS; ct++)
  {
    if (in_bin_ct[ct]>0)
    {
      (*ps)[0+3*ct] = k_ave[ct]/(float)in_bin_ct[ct];                              // Wavenumber
      (*ps)[1+3*ct] = p_box[ct]/(float)in_bin_ct[ct];                              // Power
      (*ps)[2+3*ct] = p_box[ct]/(float)in_bin_ct[ct]/sqrt((float)in_bin_ct[ct]);   // Error in power?
    }

  }
  *ps_nbins = NUM_BINS;
  *average_T = (float)ave;




  //    printf("\n");
  //    printf("21cm PS k bins\n");
  //    
  //    for (ct=1; ct<NUM_BINS; ct++)
  //    {
  //        printf("%3d   %g\n", ct, (*ps)[0+3*ct]);
  //    }
  //    printf("\n\n");




  // Deallocate
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  free(p_box);
  free(k_ave);
  free(in_bin_ct);
  fftwf_free(deldel_T);
  fclose(LOG);

  fftwf_cleanup_threads();

  return 0;
}
#endif
