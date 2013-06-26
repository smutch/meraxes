#include <math.h>
#include "meraxes.h"

//! The formation history model physics function
static double physics_func(run_globals_struct *run_globals, double prop, int snapshot)
{

  double *ZZ = run_globals->ZZ;
  
  double log_prop = log10(prop);

  double cur_peak        = run_globals->params.physics.peak * pow(1.0+ZZ[snapshot], run_globals->params.physics.peak_evo);
  double cur_stellarfrac = run_globals->params.physics.stellarfrac * pow(1.0+ZZ[snapshot], run_globals->params.physics.stellarfrac_evo);
  double cur_sigma       = run_globals->params.physics.sigma * pow(1.0+ZZ[snapshot], run_globals->params.physics.sigma_evo);

  return cur_stellarfrac * exp( -pow((log_prop-cur_peak)/cur_sigma ,2.) );

}

//! Evolve existing galaxies forward in time
int evolve_galaxies(run_globals_struct *run_globals, fof_group_struct *fof_group, int snapshot, int NGal, int NFof)
{

  double         sfr;
  double         BaryonFrac      = run_globals->params.BaryonFrac;
  double         RecycleFraction = run_globals->params.RecycleFraction;
  double         dt              = run_globals->LTTime[snapshot-1]-run_globals->LTTime[snapshot];
  galaxy_struct *gal             = NULL;
  galaxy_struct *parent          = NULL;
  halo_struct   *halo            = NULL;
  int            gal_counter     = 0;
  int            dead_gals       = 0;

  SID_log("Doing physics...", SID_LOG_OPEN|SID_LOG_TIMER);
  
  for(int i_fof=0; i_fof<NFof; i_fof++)
  {
    halo = fof_group[i_fof].FirstHalo;
    while (halo!=NULL) {
      gal = halo->Galaxy;

      while(gal!=NULL){
        if((gal->Mvir>0.0) && (gal->Type==0))
          switch (run_globals->params.physics.funcprop){
            case VMAX_PROP:
              sfr = BaryonFrac*gal->dMdt * physics_func(run_globals, gal->Vmax, snapshot);
              break;
            case MVIR_PROP:
              sfr = BaryonFrac*gal->dMdt * physics_func(run_globals, gal->Mvir*1.e10, snapshot);
              break;
            default:
              SID_log_error("Did not recognise physics_funcprop value!");
              ABORT(EXIT_FAILURE);
              break;
          }
        else
          sfr = 0.0;
      
        // update the star formation rate in the galaxy structure 
        for(int outputbin = 0; outputbin < NOUT; outputbin++)
        {
          if(snapshot == run_globals->ListOutputSnaps[outputbin])
          {
            gal->Sfr[outputbin] += sfr;
            break;
          }
        }

        // Instantaneous recycling approximation
        gal->StellarMass += (1.0-RecycleFraction)*sfr*dt;


        // If this is a type 2 then increment the merger clock
        if(gal->Type == 2)
          gal->MergTime -= dt;

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
          // merged then process a merger event.  Note that this relies on the
          // merger target coming before this galaxy in the linked list of halo
          // members.  This should be the case but I should confirm that it is
          // always true...
          if((gal->MergTime <0) || (gal->MergerTarget->Type==3))
          {
            // SID_log("Gal ID=%d has MergTime < 0", SID_LOG_COMMENT, gal->ID);

            // Merger!
            parent = gal->MergerTarget;

            while (parent->Type==3)
              parent = parent->MergerTarget;

            // Add galaxies together
            parent->StellarMass += gal->StellarMass;

            for(int outputbin = 0; outputbin < NOUT; outputbin++)
              parent->Sfr[outputbin] += gal->Sfr[outputbin];

            // Mark the merged galaxy as dead
            gal->Type          = 3;
            gal->HaloDescIndex = -1;
            dead_gals++;
          }
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

  SID_log("...done", SID_LOG_CLOSE); 

  return gal_counter-dead_gals;
}
