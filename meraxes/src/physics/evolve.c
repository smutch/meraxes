#include <math.h>
#include "meraxes.h"

//! The formation history model physics function
static double physics_func(run_globals_t *run_globals, double prop, int snapshot)
{

  double *ZZ = run_globals->ZZ;
  
  double log_prop = log10(prop);

  double cur_peak        = run_globals->params.physics.peak * pow(1.0+ZZ[snapshot], run_globals->params.physics.peak_evo);
  double cur_stellarfrac = run_globals->params.physics.stellarfrac * pow(1.0+ZZ[snapshot], run_globals->params.physics.stellarfrac_evo);
  double cur_sigma       = run_globals->params.physics.sigma * pow(1.0+ZZ[snapshot], run_globals->params.physics.sigma_evo);

  return cur_stellarfrac * exp( -pow((log_prop-cur_peak)/cur_sigma ,2.) );

}

//! Evolve existing galaxies forward in time
int evolve_galaxies(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot, int NGal, int NFof)
{

  double         sfr             = 0.0;
  double         BaryonFrac      = run_globals->params.BaryonFrac;
  double         RecycleFraction = run_globals->params.RecycleFraction;
  galaxy_t *gal             = NULL;
  halo_t   *halo            = NULL;
  int            gal_counter     = 0;
  int            dead_gals       = 0;
  double         dMdt            = 0.0;
  double         burst_time      = 0.0;
  double         burst_mass      = 0.0;
  bool           cooling_flag    = true;

#ifdef DEBUG
  int            suppressed_cooling_count = 0;
#endif

  SID_log("Doing physics...", SID_LOG_OPEN|SID_LOG_TIMER);
  
  for(int i_fof=0; i_fof<NFof; i_fof++)
  {
    halo = fof_group[i_fof].FirstHalo;
    while (halo!=NULL) {
      gal = halo->Galaxy;

      while(gal!=NULL){
        if(gal->dM > 0)
          dMdt = (gal->dM)/(gal->dt);
        else
          dMdt = 0.;

#ifdef USE_TOCF
        if((run_globals->params.TOCF_Flag) && (tocf_params.uvb_feedback))
          cooling_flag = check_reionization_cooling(run_globals, halo);
#ifdef DEBUG
        if((!cooling_flag) && (gal->Type==0))
          suppressed_cooling_count++;
#endif
#endif

        if((gal->Mvir>0.0) && (gal->Type==0) && (dMdt>0.0) && (cooling_flag))
          switch (run_globals->params.physics.funcprop){
            case VMAX_PROP:
              sfr = BaryonFrac*dMdt * physics_func(run_globals, gal->Vmax, snapshot);
              break;
            case MVIR_PROP:
              sfr = BaryonFrac*dMdt * physics_func(run_globals, gal->Mvir*1.e10, snapshot);
              break;
            default:
              SID_log_error("Did not recognise physics_funcprop value!");
              ABORT(EXIT_FAILURE);
              break;
          }
        else
          sfr = 0.0;
      
        // update the star formation rate in the galaxy structure 
        gal->Sfr = sfr;

        // Instantaneous recycling approximation
        burst_mass = (1.0-RecycleFraction)*sfr*gal->dt;
        gal->StellarMass += burst_mass;

        // Add to the luminosities due to this stellar mass burst
        burst_time   = run_globals->LTTime[snapshot] + (0.5 * gal->dt);
        // double old_lum = gal->Lum[0][0];

        add_to_luminosities(run_globals, gal, burst_mass, 0.00666, burst_time);

        // DEBUG
        // if (gal->ID==0)
        // {
        //   SID_log("*** redshift = %.2f; burst_time = %.3e; dt = %.3e; burst_mass = %.3e; gal->Lum[1][0] = %.3e; Mag = %.3e; LTTime at output = %.3e", SID_LOG_COMMENT,
        //       run_globals->ZZ[snapshot], 
        //       burst_time* run_globals->units.UnitLength_in_cm / run_globals->units.UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR / run_globals->params.Hubble_h, 
        //       gal->dt, burst_mass, gal->Lum[1][0], lum_to_mag(gal->Lum[1][0]), run_globals->LTTime[run_globals->ListOutputSnaps[0]]);
        // }

        // DEBUG
        // SID_log("z=%.2f :: burst_mass=%.2f :: Lum[0][0] %.2f", SID_LOG_COMMENT, run_globals->ZZ[snapshot], burst_mass, gal->Lum[0][0]-old_lum);

        // If this is a type 2 then increment the merger clock
        if(gal->Type == 2)
          gal->MergTime -= gal->dt;

        gal_counter++;
        gal = gal->NextGalInHalo;
      }

      halo = halo->NextHaloInFOFGroup;
    }

    // Check for mergers
    halo = fof_group[i_fof].FirstHalo;
    while(halo!=NULL) {
      gal = halo->Galaxy;
      while(gal!=NULL) {
        if(gal->Type == 2)
        {

          // If the merger clock has run out or our target halo has already
          // merged then process a merger event.
          if((gal->MergTime <0) || (gal->MergerTarget->Type==3))
            merge_with_target(run_globals, gal, &dead_gals);

        }
        gal = gal->NextGalInHalo;
      }
      halo = halo->NextHaloInFOFGroup;
    }
  }
    
  if(gal_counter+(run_globals->NGhosts) != NGal)
  {
    SID_log_error("We have not processed the expected number of galaxies...");
    SID_log("gal_counter = %d but NGal = %d", SID_LOG_COMMENT, gal_counter, NGal);
    ABORT(EXIT_FAILURE);
  }

#ifdef DEBUG
  SID_log("Suppressed cooling in %d already present galaxies.", SID_LOG_COMMENT, suppressed_cooling_count);
#endif

  SID_log("...done", SID_LOG_CLOSE); 

  return gal_counter-dead_gals;
}
