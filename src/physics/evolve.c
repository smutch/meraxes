#include "evolve.h"
#include "blackhole_feedback.h"
#include "cooling.h"
#include "core/stellar_feedback.h"
#include "core/virial_properties.h"
#include "infall.h"
#include "meraxes.h"
#include "mergers.h"
#include "reincorporation.h"
#include "star_formation.h"
#include "supernova_feedback.h"
#include <math.h>

//! Evolve existing galaxies forward in time
int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof, int* gal_counter_Pop3, int* gal_counter_Pop2, int* gal_counter_enriched)
{
  galaxy_t* gal = NULL;
  halo_t* halo = NULL;
  int gal_counter = 0;
  int dead_gals = 0;
  double infalling_gas = 0;
  double cooling_mass = 0;
  int NSteps = run_globals.params.NSteps;
  int MetalGridDim = run_globals.params.MetalGridDim;
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);

  mlog("Doing physics...", MLOG_OPEN | MLOG_TIMERSTART);
  // pre-calculate feedback tables for each lookback snapshot
  compute_stellar_feedback_tables(snapshot);

  for (int i_fof = 0; i_fof < NFof; i_fof++) {
    // First check to see if this FOF group is empty.  If it is then skip it.
    if (fof_group[i_fof].FirstOccupiedHalo == NULL)
      continue;

    infalling_gas = gas_infall(&(fof_group[i_fof]), snapshot);

    for (int i_step = 0; i_step < NSteps; i_step++) {
      halo = fof_group[i_fof].FirstHalo;
      while (halo != NULL) {
        gal = halo->Galaxy;

        while (gal != NULL) {
        
          if (Flag_Metals == true) { // Assign to newly formed galaxies metallicity of their cell according to a certain probability
            //mlog("PkProva is %f", MLOG_MESG, PowerSpectrum(0,0.1));
            //mlog("PkProva2 is %f", MLOG_MESG, PowerSpectrum(0,0.01));
            //double M75_0 = calculate_Mvir_2(7.5, 0);
            //double M8_0 = calculate_Mvir_2(8., 0);
            double M7_0 = calculate_Mvir_2(7.,0);
            double M7_10 = calculate_Mvir_2(7.,10);
            double M7_15 = calculate_Mvir_2(7.,15);
            double M7_20 = calculate_Mvir_2(7.,20);
            double M7_25 = calculate_Mvir_2(7.,25);
            //double pro = Sigma(0 , 2.752e4 * 0.6715); // Check S8
            //double pro2 = Sigma(10 , 2.752e4 * 0.6715); // Check S8
            //double pro3 = Sigma(10 , 2.752e3 * 0.6715); // Check S8
            //mlog("GF at z = 0 is %f", MLOG_MESG, pro);
            //mlog("GF at z = 1 is %f", MLOG_MESG, pro2);
            //mlog("GF at z = 2 is %f", MLOG_MESG, pro3);  
            //mlog("S8 at z = 0 is %f", MLOG_MESG, Sigma(0, M8_0));
            mlog("S7 at z = 0 is %f", MLOG_MESG, Sigma(0, M7_0));
            //mlog("S75 at z = 0 is %f", MLOG_MESG, Sigma(0, M75_0));
            mlog("S7 at z = 10 is %f", MLOG_MESG, Sigma(10, M7_10));
            mlog("S7 at z = 15 is %f", MLOG_MESG, Sigma(15, M7_15));
            mlog("S7 at z = 20 is %f", MLOG_MESG, Sigma(20, M7_20));
            mlog("S7 at z = 25 is %f", MLOG_MESG, Sigma(10, M7_25));
            //mlog("Sigma at z = 10 is %f", MLOG_MESG, pro3);
            if (gal->output_index == -1) { //Not sure if Ghostflag condition is needed
              double x;
              
              x = (double)rand() / RAND_MAX;
              
              if (x <= gal->Metal_Probability) {
                gal->MetalsHotGas = gal->HotGas * gal->Metallicity_IGM;
                gal->MetalsColdGas = gal->ColdGas * gal->Metallicity_IGM;
                gal->MetalsEjectedGas = gal->EjectedGas * gal->Metallicity_IGM;
                
                *gal_counter_enriched = *gal_counter_enriched + 1;
                if ((gal->Metallicity_IGM / 0.01) > 1e-4) {
                  *gal_counter_Pop2 = *gal_counter_Pop2 + 1;
                  gal->Galaxy_Population = 2;
                  }
                else
                  gal->Galaxy_Population = 3; // Enriched but not enough
              }
              
              else {
                gal->Galaxy_Population = 3;
                *gal_counter_Pop3 = *gal_counter_Pop3 + 1;
              }
            }
          }
          
          if (gal->Type == 0) {
            cooling_mass = gas_cooling(gal);

            add_infall_to_hot(gal, infalling_gas / ((double)NSteps));

            reincorporate_ejected_gas(gal);

            cool_gas_onto_galaxy(gal, cooling_mass);
          }

          if (gal->Type < 3) {
            if (!Flag_IRA)
              delayed_supernova_feedback(gal, snapshot);

            if (gal->BlackHoleAccretingColdMass > 0)
              previous_merger_driven_BH_growth(gal);

            insitu_star_formation(gal, snapshot);
            
            if (Flag_Metals == true) 
              calc_metal_bubble(gal, snapshot);
              
            // If this is a type 2 then decrement the merger clock
            if (gal->Type == 2)
              gal->MergTime -= gal->dt;
          }

          if (i_step == NSteps - 1)
            gal_counter++;

          gal = gal->NextGalInHalo;
        }

        halo = halo->NextHaloInFOFGroup;
      }

      // Check for mergers
      halo = fof_group[i_fof].FirstHalo;
      while (halo != NULL) {
        gal = halo->Galaxy;
        while (gal != NULL) {
          if (gal->Type == 2)
            // If the merger clock has run out or our target halo has already
            // merged then process a merger event.
            if ((gal->MergTime < 0) || (gal->MergerTarget->Type == 3))
              merge_with_target(gal, &dead_gals, snapshot);

          gal = gal->NextGalInHalo;
        }
        halo = halo->NextHaloInFOFGroup;
      }
    }
  }

  if (gal_counter + (run_globals.NGhosts) != NGal) {
    mlog_error("We have not processed the expected number of galaxies...");
    mlog("gal_counter = %d but NGal = %d", MLOG_MESG, gal_counter, NGal);
    ABORT(EXIT_FAILURE);
  }
  
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  return gal_counter - dead_gals;
}

void passively_evolve_ghost(galaxy_t* gal, int snapshot)
{
  // Passively evolve ghosts.
  // Currently, this just means evolving their stellar pops...

  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);

  if (!Flag_IRA)
    delayed_supernova_feedback(gal, snapshot);
    
  if (Flag_Metals == true) // You are updating this function to test why probability is decreasing in some cells
    calc_metal_bubble(gal, snapshot); 
}
