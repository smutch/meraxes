/*
 * This code is an amalgamation of the requisite functions for LW background using
 * formulation of Qin2020a.
 * Written to work within Meraxes by Emanuele M. Ventura.
 */
 
 //>>>>>>Here you need to put all the includes>>>>>>//
 
#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>

#include <gsl/gsl_errno.h> //These are for integration
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>

#include "JLW.h"
#include "reionization.h" // I need this for the grid, not sure if this must be the same. Actually the grid is present also in meraxes.h 
#include "XRayHeatingFunctions.h" // Useful for zmax, drdz, dtdz and spectral_emissivity
#include "meraxes.h" // Main library. Here are defined some constants and all the run_globals ecc.
#include "misc_tools.h" // Random stuff, (e.g. INDEX_REAL and INDEX_PADDED).
#include "utils.h" // Timer stuff

 //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
 
 void _ComputeJLW(int snapshot)
 {
  double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc
  int ReionGridDim = run_globals.params.ReionGridDim; // Can I use the same grid I use for reionization? I would say yes, but maybe it's better to double check this
  double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc)^3 Remember that you need to convert this in cm^3
  double total_n_cells = pow((double)ReionGridDim, 3);
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]); //fft stuff
  int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]); //fft stuff
  run_units_t* units = &(run_globals.units); // meraxes.h
  
  double red = run_globals.ZZ[snapshot];
  double prev_red;
  if (snapshot == 0) {
    prev_red = run_globals.ZZ[snapshot];
  } else {
    prev_red = run_globals.ZZ[snapshot - 1];
  }
  
  int i_real, i_padded, i_smoothedSFR, R_ct, n_ct;
  
  double prev_zpp, prev_R, zpp, zp, R_factor, R, dzp, nuprime_LW;
  
  int TsNumFilterSteps = run_globals.params.TsNumFilterSteps; //Is the Number of Filter Step the same?? If not I need to define a new one (and maybe put definition in meraxes.h)
  
  double R_values[TsNumFilterSteps];
  double dt_dzpp_list[TsNumFilterSteps];
  
  double freq_int_pop2[TsNumFilterSteps]; 

  double result[1]; //Risultato dell'integrale. Per il momento te ne basta uno che hai solo le Pop2
  double SFR_POP2[TsNumFilterSteps]; // For now I'll just assume everything is Pop.II. What I want to do in the end is to differentiate between Pop.II and Pop.III. I could also differentiate 					         between AC and MC (like 21CMFAST works)
  
  float* JLW_box = run_globals.reion_grids.JLW_box;
  
  float* sfr = run_globals.reion_grids.sfr; // This is the SFR that we compute within Meraxes
  
  fftwf_complex* sfr_unfiltered = run_globals.reion_grids.sfr_unfiltered; // Filtering stuff
  fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
  fftwf_execute(run_globals.reion_grids.sfr_forward_plan);
  
  for (int ii = 0; ii < slab_n_complex; ii++) { // MPI stuff
    sfr_unfiltered[ii] /= (float)total_n_cells;
  }
  
  double* SMOOTHED_SFR_POP2 = run_globals.reion_grids.SMOOTHED_SFR_POP2; // Definition of this stuff in meraxes.h and reionization.c
  
  init_LW();
  
  double J_LW_ave;
  J_LW_ave = 0.0;
  
  zp = red; 
  dzp = zp - prev_red;
  dt_dzp = dtdz((float)(float)zp); 
  
  // Setup starting radius (minimum) and scaling to obtaining the maximum filtering radius for the LW background
  R = L_FACTOR * box_size / (float)ReionGridDim; // Take CARE that here you are doing the same than X-rays! Make a double check!
  R_factor = pow(R_XLy_MAX / R, 1 / (float)TsNumFilterSteps);

  // Smooth the density, stars and SFR fields over increasingly larger filtering radii (for evaluating the LW integrals)
  for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {

      R_values[R_ct] = R;

      memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex); // Copies fftwf_complex*slab_n complex characters from sfr_unfiltered to sfr_filtered.

      if (R_ct > 0) {
        int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

        filter(sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.TsHeatingFilterType); //Or ReionFilterType??
      }

      // inverse fourier transform back to real space
      fftwf_execute(run_globals.reion_grids.sfr_filtered_reverse_plan);

        for (int ix = 0; ix < local_nix; ix++)
          for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
              i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
              i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
              i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

              ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

              SMOOTHED_SFR_POP2[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) * // SFRD (Msolar/yr/Mpc^3)
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *  //SFRD (g/s/Mpc^3)
                                                pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS; //(g/s/cm^(-3)) (I think you should divide by PROTONMASS, did that in evolveLW)
      }
      R *= R_factor;
    }

      // Populate the initial LW tables
      for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {
  
        if (R_ct == 0) {
          prev_zpp = zp;
          prev_R = 0;
        } else {
          prev_zpp = zpp_edgee[R_ct - 1];
          prev_R = R_values[R_ct - 1];
        }

      zpp_edgee[R_ct] = prev_zpp - (R_values[R_ct] - prev_R) * MPC / (drdz((float)prev_zpp)); // cell size, drdz defined in XrayHeatingFunctions.c
      zpp = (zpp_edgee[R_ct] + prev_zpp) * 0.5; // average redshift value of shell: z'' + 0.5 * dz''

      dt_dzpp_list[R_ct] = dtdz((float)zpp);
      dt_dzpp = dt_dzpp_list[R_ct];	
      
      sum_lyn_LW[R_ct] = 0;
      for (n_ct = NSPEC_MAX; n_ct >= 2; n_ct--) {
        if (zpp > zmax((float)zp, n_ct))
          continue;

        nuprime_LW = nu_n(n_ct) * (1 + zpp) / (1.0 + zp);
        sum_lyn_LW[R_ct] += frecycle(n_ct) * spectral_emissivity_LW(nuprime_LW, 0, 2); //2 è per la Pop. //Check if it works
        freq_int_pop2[R_ct] = sum_lyn_LW[R_ct] * PLANCK; //erg * s
      }
    }
      
      for (int ix = 0; ix < local_nix; ix++)
          for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
              i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
              i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
              //i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);
              
              //SFR_POP2[R_ct] = SMOOTHED_SFR_POP2[i_smoothedSFR]; // Do I use this to move from Fourier Space to real space?
              //freq_int_pop2[R_ct] = 0.0;
      
      	      //for (n_ct = NSPEC_MAX; n_ct >= 2; n_ct--) {
              //  if (zpp > zmax((float)zp, n_ct)) //zmax defined in XrayHeatingFunctions.c
              //   continue;
              //   }          
                //freq_int_pop2[R_ct] += nu_integral(n_ct, zp, zpp, SFR_POP2[R_ct]);
                //freq_int_pop2[R_ct] += SFR_POP2[R_ct] * sum_lyn_LW[R_ct] * PLANCK_EV
                
              for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {
                i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

                SFR_POP2[R_ct] = SMOOTHED_SFR_POP2[i_smoothedSFR]/PROTONMASS; // g/(s*cm^-3) (DIVIDING BY PROTONMASS I make it for baryon
                dt_dzpp = dt_dzpp_list[R_ct]; //s
                }  
                 
              evolveLW((float)zp, freq_int_pop2, SFR_POP2, result); //baryon/(s*cm^-3)*s*(cm/s)*erg*s = baryon*cm^-2*Hz^-1*s^-1*erg*sr^-1
              
              JLW_box[i_real] = result[0]; 
              J_LW_ave += JLW_box[i_real];  
       }      
       MPI_Allreduce(MPI_IN_PLACE, &J_LW_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       J_LW_ave /= total_n_cells;
       run_globals.reion_grids.volume_ave_J_LW = J_LW_ave;
          
       destruct_LW(); 
       mlog("zp = %e J_LW_ave = %e", MLOG_MESG, zp, J_LW_ave);	  
}

/* typedef struct
{
  double StarFRateDens; // Not sure how to treat emiss since it's a function of nu
} nu_integral_params;
 
 double nu_integrand(double nu, void* params) 
 {
  
  double emissivity;
  
  nu_integral_params* p = (nu_integral_params*)params;
  emissivity = spectral_emissivity(nu , 0); // Function defined in XrayHeatingFunctions.h Careful! It may be nuprime instead of n! Check better spectral_emissivity
  
  return emissivity * PLANCK_EV * p->StarFRateDens;
   
 } 
 
 double nu_integral(int nn, float redd, float redprime, double SFRD) // This function should be ok!
 {
 
  double result, error;
  double rel_tol = 0.01;
  double nun, nun_next, nuprime;
  
  nun = NU_LL * (1.0 - pow(nn, -2));
  nun_next = NU_LL * (1.0 - pow(nn+1, -2));
  nuprime = nun * (1 + redprime) / (1 + redd);
   
  gsl_function F; // This data type defines a general function with parameters, to use this I need a structure with parameters of the function.
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000); // This function allocates a workspace sufficient to hold 1000 double precision intervals, their results and error. 

  nu_integral_params p;
  
  p.StarFRateDens = SFRD;
  
  F.params = &p;
  F.function = &nu_integrand;
  
  gsl_integration_qag (&F, fmax(NU_LW, nuprime), nun_next, 0, rel_tol, 1000, GSL_INTEG_GAUSS61, w, &result, &error); // Function to do the integral
  gsl_integration_workspace_free(w);
  
  return result;
  
 } */
 
 // Reads in and constructs table of the piecewise power-law fits to Pop 2 and Pop 3 stellar spectra, from Barkana, specifically designed for LW computation. Almost sure that is per baryon
 double spectral_emissivity_LW(double nu_norm, int flag, int Population)
 {
  static int n[NSPEC_MAX];
  static float nu_n[NSPEC_MAX], alpha_S_2[NSPEC_MAX];
  static float alpha_S_3[NSPEC_MAX], N0_2[NSPEC_MAX], N0_3[NSPEC_MAX];
  double n0_fac;
  double ans;
  double lower_limit;
  int i;
  FILE* F;
  
  char fname[STRLEN];
  
  if (flag == 1) {

    if (run_globals.mpi_rank == 0) {

      sprintf(fname, "%s/stellar_spectra.dat", run_globals.params.TablesForXHeatingDir);

      // Read in the data
      if (!(F = fopen(fname, "r"))) {
        mlog("spectral_emissivity: Unable to open file: stellar_spectra.dat at %s for reading\nAborting\n",
             MLOG_MESG,
             fname);
        return -1;
      }

      for (i = 1; i < NSPEC_MAX; i++) {
        fscanf(F, "%i %e %e %e %e", &n[i], &N0_2[i], &alpha_S_2[i], &N0_3[i], &alpha_S_3[i]);
      }
      fclose(F);

      for (i = 1; i < NSPEC_MAX; i++) {
        nu_n[i] = (float)(4.0 / 3.0 * (1.0 - 1.0 / pow(n[i], 2.0)));
      }

      for (i = 1; i < NSPEC_MAX; i++) {
        nu_n[i] = (float)(4.0 / 3.0 * (1.0 - 1.0 / pow(n[i], 2.0)));
      }

      for (i = 1; i < (NSPEC_MAX - 1); i++) {
        n0_fac = (pow(nu_n[i + 1], alpha_S_2[i] + 1) - pow(nu_n[i], alpha_S_2[i] + 1));
        N0_2[i] *= (alpha_S_2[i] + 1) / n0_fac * Pop2_ion;
        n0_fac = (pow(nu_n[i + 1], alpha_S_3[i] + 1) - pow(nu_n[i], alpha_S_3[i] + 1));
        N0_3[i] *= (alpha_S_3[i] + 1) / n0_fac * Pop3_ion;
      }
    }

    // broadcast the values to all cores
    MPI_Bcast(nu_n, sizeof(nu_n), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(alpha_S_2, sizeof(alpha_S_2), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(alpha_S_3, sizeof(alpha_S_3), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(N0_2, sizeof(N0_2), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(N0_3, sizeof(N0_3), MPI_BYTE, 0, run_globals.mpi_comm);

    return 0.0;
  }
  
  ans = 0.0;
  for (i = 1; i < (NSPEC_MAX - 1); i++) {
    if (nu_n[i] >= NU_LW)
      lower_limit = nu_n[i];
    else
      lower_limit = NU_LW; // Questa parte la puoi riscrivere come l'ha scritta Yuxiang
    if ((nu_norm >= lower_limit) && (nu_norm < nu_n[i + 1])) {
      // We are in the correct spectral region
      if (Population == 2)
        ans = N0_2[i] * pow(nu_norm, alpha_S_2[i]);
      else
        ans = N0_3[i] * pow(nu_norm, alpha_S_3[i]);

      return ans / Ly_alpha_HZ;
    }
  }

  i = NSPEC_MAX - 1;
  if (Pop == 2)
    return N0_2[i] * pow(nu_norm, alpha_S_2[i]) / Ly_alpha_HZ;
  else
    return N0_3[i] * pow(nu_norm, alpha_S_3[i]) / Ly_alpha_HZ;
}
 
 int init_LW() // Maybe you don't need this but just try to see if it works
 {
   size_t TsNumFilterSteps = (size_t)run_globals.params.TsNumFilterSteps;
   zpp_edgee = calloc(TsNumFilterSteps, sizeof(double));
   sum_lyn_LW = calloc(TsNumFilterSteps, sizeof(double));
   
   if (spectral_emissivity_LW(0, 1, 0) < 0)
    return -6;
 
   return 0;
 }
 
 void destruct_LW()
{
  free(zpp_edgee);
  spectral_emissivity_LW(0.0, 2, 0); //2 is the flag, frees memory.
  free(sum_lyn_LW);
}
 
 void evolveLW(float zp, const double integrand_POP2[], const double StarF_POP2[], double deriv[])
 {
 
  double dlw_dt_POP2;
  double zpp, dzpp;
  double zpp_integrand_POP2;
  int zpp_ct;
   
  dlw_dt_POP2 = 0.;
      
  for (zpp_ct = 0; zpp_ct < run_globals.params.TsNumFilterSteps; zpp_ct++) {
     // set redshift of half annulus; dz'' is negative since we flipped limits of integral
     if (zpp_ct == 0) {
       zpp = (zpp_edgee[0] + zp) * 0.5;
       dzpp = zp - zpp_edgee[0];
     } else {
       zpp = (zpp_edgee[zpp_ct] + zpp_edgee[zpp_ct - 1]) * 0.5;
       dzpp = zpp_edgee[zpp_ct - 1] - zpp_edgee[zpp_ct];
     }
      
     zpp_integrand_POP2 = integrand_POP2[zpp_ct];
     dlw_dt_POP2 += dt_dzpp * dzpp * StarF_POP2[zpp_ct] * zpp_integrand_POP2 * pow(1 + zp, 2) * (1 + zpp);     
  }
  
  dlw_dt_POP2 *= (SPEED_OF_LIGHT / (4. * M_PI)) / (PROTONMASS / SOLAR_MASS); //cm/s*s*(g/s*cm^-3)*erg*s*Hz (it should be cm*g*Hz*eV)
  
  deriv[0] = dlw_dt_POP2; // This is your final result, in the future you will disentangle between different components
}   
   
// This function makes sure that the right version of ComputeJLW() gets called.
// Note: Only the CPU version works for now
void ComputeJLW(int snapshot, timer_info* timer_total)
{
  // Call the version of ComputeTs we've been passed (and time it)
  timer_info timer;
  // double redshift = run_globals.ZZ[snapshot];
  // Run the Meraxes version of _ComputeJLW()
  mlog("Calling pure-CPU version of ComputeJLW() for snap=%d/z=%.2lf...",
       MLOG_OPEN | MLOG_TIMERSTART,
       snapshot,
       run_globals.ZZ[snapshot]);
  timer_start(&timer);
  _ComputeJLW(snapshot);
  timer_stop(&timer);
  timer_stop(timer_total);
  timer_gpu += timer_delta(timer);
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
  mlog("Total time spent in ComputeJLW vs. total run time (snapshot %d ): %.2f of %.2f s",
       MLOG_MESG,
       snapshot,
       timer_gpu,
       timer_delta(*timer_total));
}  
   
