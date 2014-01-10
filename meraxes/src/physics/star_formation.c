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


void form_stars_insitu(run_globals_t *run_globals, galaxy_t *gal, int snapshot, bool cooling_flag)
{

  double dMdt;
  double sfr             = 0.0;
  double BaryonFrac      = run_globals->params.BaryonFrac;
  double RecycleFraction = run_globals->params.RecycleFraction;
  double burst_time      = 0.0;
  double burst_mass      = 0.0;

  if(gal->dM > 0)
    dMdt = (gal->dM)/(gal->dt);
  else
    dMdt = 0.;

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

  add_to_luminosities(run_globals, gal, burst_mass, 0.00666, burst_time);
}
