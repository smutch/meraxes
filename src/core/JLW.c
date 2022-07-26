/*
 * This code is an amalgamation of the requisite functions for LW background using
 * formulation of Qin2020a.
 * Written to work within Meraxes by Emanuele M. Ventura.
 */
 
#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>

#include <gsl/gsl_errno.h> 
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>

#include "JLW.h"
#include "reionization.h"  
#include "XRayHeatingFunctions.h"
#include "meraxes.h"
#include "misc_tools.h" 
#include "utils.h"

 //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
 
 void _ComputeJLW(int snapshot)
 {
  double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; 
  int ReionGridDim = run_globals.params.ReionGridDim; 
  double pixel_volume = pow(box_size / (double)ReionGridDim, 3); 
  double total_n_cells = pow((double)ReionGridDim, 3);
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]); 
  int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
  run_units_t* units = &(run_globals.units); 
  
  double red = run_globals.ZZ[snapshot];
  double prev_red;
  if (snapshot == 0) {
    prev_red = run_globals.ZZ[snapshot];
  } else {
    prev_red = run_globals.ZZ[snapshot - 1];
  }
  
  int i_real, i_padded, i_smoothedSFR, R_ct, n_ct;
  
  double prev_zpp, prev_R, zpp, zp, R_factor, R, dzp, nuprime_LW;
  
  int TsNumFilterSteps = run_globals.params.TsNumFilterSteps; 
  
  double R_values[TsNumFilterSteps];
  double dt_dzpp_list[TsNumFilterSteps];
  
  double freq_int_pop2[TsNumFilterSteps]; 

  double result[1]; 
  double SFR_POP2[TsNumFilterSteps]; 
  
  float* JLW_box = run_globals.reion_grids.JLW_box;
  
  float* sfr = run_globals.reion_grids.sfr; 
  
  fftwf_complex* sfr_unfiltered = run_globals.reion_grids.sfr_unfiltered; 
  fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
  fftwf_execute(run_globals.reion_grids.sfr_forward_plan);
  
  for (int ii = 0; ii < slab_n_complex; ii++) { 
    sfr_unfiltered[ii] /= (float)total_n_cells;
  }
  
  double* SMOOTHED_SFR_POP2 = run_globals.reion_grids.SMOOTHED_SFR_POP2; 
  
  init_LW();
  
  double J_LW_ave;
  J_LW_ave = 0.0;
  
  zp = red; 
  dzp = zp - prev_red;
  dt_dzp = dtdz((float)(float)zp); 
  
  R = L_FACTOR * box_size / (float)ReionGridDim; 
  R_factor = pow(R_XLy_MAX / R, 1 / (float)TsNumFilterSteps);

  for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {

      R_values[R_ct] = R;

      memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex); 

      if (R_ct > 0) {
        int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

        filter(sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.TsHeatingFilterType); 
      }

      fftwf_execute(run_globals.reion_grids.sfr_filtered_reverse_plan);

        for (int ix = 0; ix < local_nix; ix++)
          for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
              i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
              i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
              i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

              ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

              //SMOOTHED_SFR_POP2[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) * 
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *  
                                                pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
              SMOOTHED_SFR_POP2[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) * 
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *  
                                                pow(units->UnitLength_in_cm, -3.);                                  
      }
      R *= R_factor;
    }

      for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {
  
        if (R_ct == 0) {
          prev_zpp = zp;
          prev_R = 0;
        } else {
          prev_zpp = zpp_edgee[R_ct - 1];
          prev_R = R_values[R_ct - 1];
        }

      zpp_edgee[R_ct] = prev_zpp - (R_values[R_ct] - prev_R) * MPC / (drdz((float)prev_zpp)); 
      zpp = (zpp_edgee[R_ct] + prev_zpp) * 0.5; 

      dt_dzpp_list[R_ct] = dtdz((float)zpp);
      dt_dzpp = dt_dzpp_list[R_ct];	
      
      sum_lyn_LW[R_ct] = 0;
      for (n_ct = NSPEC_MAX; n_ct >= 2; n_ct--) {
        if (zpp > zmax((float)zp, n_ct))
          continue;

        nuprime_LW = nu_n(n_ct) * (1 + zpp) / (1.0 + zp); // This is in Lyman-alpha units!
        sum_lyn_LW[R_ct] += spectral_emissivity_LW(nuprime_LW, 0, 2);
        //freq_int_pop2[R_ct] = sum_lyn_LW[R_ct] * PLANCK;
      }
      freq_int_pop2[R_ct] = sum_lyn_LW[R_ct] * PLANCK;
    }
      
      for (int ix = 0; ix < local_nix; ix++)
          for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
              i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
              i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                
              for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {
                i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

                SFR_POP2[R_ct] = SMOOTHED_SFR_POP2[i_smoothedSFR] / PROTONMASS;
                dt_dzpp = dt_dzpp_list[R_ct]; //s
                }  
                 
              evolveLW((float)zp, freq_int_pop2, SFR_POP2, result); 
              
              JLW_box[i_real] = result[0] * 1e21; // Compute LW in 1e-21 units 
              J_LW_ave += JLW_box[i_real];  
       }      
       MPI_Allreduce(MPI_IN_PLACE, &J_LW_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       J_LW_ave /= total_n_cells;
       run_globals.reion_grids.volume_ave_J_LW = J_LW_ave;
          
       destruct_LW(); 
       mlog("zp = %e J_LW_ave = %e", MLOG_MESG, zp, J_LW_ave);	  
}

 double spectral_emissivity_LW(double nu_norm, int flag, int Population) // if it works incorporate this inside spectral_emissivity
 {
   static int n[NSPEC_MAX];
   static float nu_n[NSPEC_MAX];
   static float alpha_S_3[NSPEC_MAX], N0_2[NSPEC_MAX], N0_3[NSPEC_MAX], alpha_S_2[NSPEC_MAX];
   double n0_fac;
   double ans;
   int i;
   double lower_limit;
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
    
    MPI_Bcast(nu_n, sizeof(nu_n), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(alpha_S_2, sizeof(alpha_S_2), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(alpha_S_3, sizeof(alpha_S_3), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(N0_2, sizeof(N0_2), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(N0_3, sizeof(N0_3), MPI_BYTE, 0, run_globals.mpi_comm);

    return 0.0;
  }
   
   for (i = 1; i < (NSPEC_MAX -1); i++) {
     if (nu_norm[i] >= NU_LW / Ly_alpha_HZ)
       lower_limit = nu_norm[i];
     else
       lower_limit = NU_LW / Ly_alpha_HZ;
     if ((nu_norm >= lower_limit) && (nu_norm < nu_n[i+1])) {
       if (Population == 2) {
         ans = N0_2[i]  / (alpha_S_2[i] + 1) * (pow(nu_n[i+1], alpha_S_2[i]+1) - pow(nu_norm, alpha_S_2[i]+1));
       }
       else if (Population == 3) {
         ans = N0_3[i] / (alpha_S_3[i] + 1) * (pow(nu_n[i+1], alpha_S_3[i]+1) - pow(nu_norm, alpha_S_3[i]+1));
       }
       else {
         mlog("Invalid value for Stellar Population", MLOG_MESG);
       }
       return ans;
     }
   }
 }
 
 int init_LW() 
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
  spectral_emissivity_LW(0.0, 2, 0); 
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
  
  //dlw_dt_POP2 *= (SPEED_OF_LIGHT / (4. * M_PI)) / (PROTONMASS / SOLAR_MASS); 
  dlw_dt_POP2 *= (SPEED_OF_LIGHT / (4. * M_PI));
  
  deriv[0] = dlw_dt_POP2; 
}   
   
void ComputeJLW(int snapshot, timer_info* timer_total)
{
  
  timer_info timer;
  
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
   
