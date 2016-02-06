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
    float         *xH,
    float         *deltax,
    float         *v,                 // Set to NULL if no velocity field
    float         *Ts,                // Set to NULL if no spin temperature field
    float         *average_T,
    float         *delta_T,
    float         **ps,
    int           *ps_nbins)
{
  fftwf_complex *deldel_T;
  fftwf_plan plan;

  bool flag_use_vel = false;
  bool flag_use_ts = false;

  int i, j, k, n_x, n_y, n_z, NUM_BINS;

  unsigned long long *in_bin_ct;
  unsigned long long ct, temp_ct;

  float pixel_x_HI, pixel_deltax;
  float maxi, maxj, maxk, maxdvdx, mini, minj, mink, mindvdx;
  float k_x, k_y, k_z, k_mag, k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
  float T_rad, pixel_Ts_factor;

  double *p_box, *k_ave;
  double dvdx, max_v_deriv;
  double ave_Ts, min_Ts, max_Ts, temp;


  // Begin initialisation
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  float min       = 1e3;
  float max       = -1e3;
  double ave       = 0.0;
  unsigned long long nonlin_ct = 0;

  // If no velocity field is passed then don't use velocities!
  if(v != NULL)
    flag_use_vel = true;

  // If no spin temperature field is passed then don't use spin temperature!
  if(Ts != NULL)
    flag_use_ts = true;

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

  for (i=0; i<HII_DIM; i++)
  {
    fprintf(LOG, "%i...", i);
    for (j=0; j<HII_DIM; j++)
    {
      for (k=0; k<HII_DIM; k++)
      {
        pixel_deltax = deltax[HII_R_FFT_INDEX(i,j,k)];   //*((float *)deltax + HII_R_FFT_INDEX(i,j,k));
        pixel_x_HI = xH[HII_R_INDEX(i,j,k)];

        if (pixel_x_HI > TINY)
        {
          temp = pixel_deltax;
          temp_ct++;
        }

        delta_T[HII_R_INDEX(i,j,k)] = const_factor*pixel_x_HI*(1.0 + pixel_deltax);

        if (USE_TS_IN_21CM)
        {
          pixel_Ts_factor = 1.0 - T_rad/Ts[HII_R_INDEX(i,j,k)];
          delta_T[HII_R_INDEX(i,j,k)] *= pixel_Ts_factor;
          ave_Ts += Ts[HII_R_INDEX(i,j,k)];

          if (min_Ts > Ts[HII_R_INDEX(i,j,k)]) {min_Ts = Ts[HII_R_INDEX(i,j,k)];}

          if (max_Ts < Ts[HII_R_INDEX(i,j,k)]) {max_Ts = Ts[HII_R_INDEX(i,j,k)];}

        }

        if (max < delta_T[HII_R_INDEX(i,j,k)]) {max = delta_T[HII_R_INDEX(i,j,k)];}

        if (min > delta_T[HII_R_INDEX(i,j,k)]) {min = delta_T[HII_R_INDEX(i,j,k)];}

        ave += delta_T[HII_R_INDEX(i,j,k)];
      }
    }
  }
  fprintf(LOG, "\n%e\n", temp/(double)temp_ct);
  ave /= (double)HII_TOT_NUM_PIXELS;
  fprintf(LOG, "Without velocities, max is %e, min is %e, ave is %e\n", max, min, ave);


  //    printf("    delta_T max    = %g\n", max);
  //    printf("    delta_T min    = %g\n", min);
  //    printf("    delta_T ave    = %g\n", ave);


  if (USE_TS_IN_21CM)
  {
    ave_Ts /= (double)(HII_TOT_NUM_PIXELS);
    fprintf(LOG, "Ts, min = %e, max = %e, ave = %e\n", min_Ts, max_Ts, ave_Ts);
    fprintf(LOG, "corresponding to (1-trad/Ts of), min = %e, max = %e, ave = %e\n", 1-T_rad/min_Ts, 1-T_rad/max_Ts, 1-T_rad/ave_Ts);
  }

  fflush(LOG);


  // Deal with velocities
  // ------------------------------------------------------------------------------------------------------

  if(T_USE_VELOCITIES)
  {
    min = 1e3;
    max = -1.0;
    ave = 0.0;

    // Take the derivative in k-space
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)v, (fftwf_complex *)v, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    for (n_x=0; n_x<HII_DIM; n_x++)
    {
      if (n_x>HII_MIDDLE)
      {
        k_x = (n_x-HII_DIM)*DELTA_K;  // Wrap around for FFT convention
      }
      else
      {
        k_x = n_x*DELTA_K;
      }

      for (n_y=0; n_y<HII_DIM; n_y++)
      {
        if (n_y>HII_MIDDLE)
        {
          k_y = (n_y-HII_DIM)*DELTA_K;
        }
        else
        {
          k_y = n_y*DELTA_K;
        }

        for (n_z=0; n_z<=HII_MIDDLE; n_z++)
        { 
          k_z = n_z*DELTA_K;

          // Take partial deriavative along the line of sight
          switch(VELOCITY_COMPONENT)
          {
            case 1:
              *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_x*I/(float)HII_TOT_NUM_PIXELS;
              break;
            case 3:
              *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_z*I/(float)HII_TOT_NUM_PIXELS;
              break;
            default:
              *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_y*I/(float)HII_TOT_NUM_PIXELS;
          }

        }
      }
    }
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)v, (float *)v, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();


    // Now add the velocity correction to the delta_T maps

    max_v_deriv = fabs(MAX_DVDR*H);

    for (i=0; i<HII_DIM; i++)
    {
      for (j=0; j<HII_DIM; j++)
      {
        for (k=0; k<HII_DIM; k++)
        {
          dvdx = v[HII_R_FFT_INDEX(i,j,k)];

          // Set maximum allowed gradient for this linear approximation
          if (fabs(dvdx) > max_v_deriv)
          {
            if (dvdx < 0)
            {
              dvdx = -max_v_deriv;
            }
            else
            {
              dvdx = max_v_deriv;
            }

            nonlin_ct++;
          }

          delta_T[HII_R_INDEX(i,j,k)] /= (dvdx/H + 1.0);

          if (max < delta_T[HII_R_INDEX(i,j,k)])
          {
            maxi = i;
            maxj = j;
            maxk = k;
            maxdvdx = dvdx;
            max = delta_T[HII_R_INDEX(i,j,k)];
          }

          if (min > delta_T[HII_R_INDEX(i,j,k)])
          {
            mini = i;
            minj = j;
            mink = k;
            mindvdx = dvdx;
            min = delta_T[HII_R_INDEX(i,j,k)];
          }

          ave += delta_T[HII_R_INDEX(i,j,k)];
        }
      }
    }
    ave /= (double)HII_TOT_NUM_PIXELS;

    fprintf(LOG, "With velocities:\nMax is %e\t dvdx is %e, ave is %e\n", max, maxdvdx, ave);
    fprintf(LOG, "%llu out of %llu voxels (fraction=%e) exceeded max allowed velocity gradient\n", nonlin_ct, HII_TOT_NUM_PIXELS, nonlin_ct/(double)HII_TOT_NUM_PIXELS);
    fprintf(LOG, "Min is %e\t dvdx is %e\n\n", min, mindvdx);
  }


  // Calculate power spectrum
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  k_factor = 1.5;
  k_first_bin_ceil = DELTA_K;
  k_max = DELTA_K*HII_DIM;

  // Initialise arrays
  // ghetto counting (lookup how to do logs of arbitrary bases in c...)
  NUM_BINS = 0;
  k_floor = 0.0;
  k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max)
  {
    NUM_BINS++;
    k_floor = k_ceil;
    k_ceil *= k_factor;
  }

  //    printf("    NUM_BINS       = %d\n", NUM_BINS);

  p_box     = (double *)malloc(sizeof(double)*NUM_BINS);
  k_ave     = (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);

  if (!p_box || !in_bin_ct || !k_ave) // a bit sloppy, but whatever..
  {
    fprintf(stderr, "delta_T.c: Error allocating memory.\nAborting...\n");
    fprintf(LOG, "delta_T.c: Error allocating memory.\nAborting...\n");
    fclose(LOG);
    fftwf_cleanup_threads();
    return -1;
  }
  for (ct=0; ct<NUM_BINS; ct++)
  {
    p_box[ct] = 0.0;
    k_ave[ct] = 0.0;
    in_bin_ct[ct] = 0;
  }

  deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

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
