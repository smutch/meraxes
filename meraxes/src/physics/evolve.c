#include <math.h>
#include "meraxes.h"

// debug
static void print_galaxy(galaxy_t *gal)
{
  fprintf(stderr, "-------------------\n");
  fprintf(stderr, "ID = %d\n", gal->ID);
  fprintf(stderr, "Type = %d\n", gal->Type);
  fprintf(stderr, "Len = %d\n", gal->Len);
  fprintf(stderr, "Spin = %.3e\n", gal->Spin);
  fprintf(stderr, "HotGas = %.3e\n", gal->HotGas);
  fprintf(stderr, "MetalsHotGas = %.3e\n", gal->MetalsHotGas);
  fprintf(stderr, "ColdGas = %.3e\n", gal->ColdGas);
  fprintf(stderr, "MetalsColdGas = %.3e\n", gal->MetalsColdGas);
  fprintf(stderr, "StellarMass = %.3e\n", gal->StellarMass);
  fprintf(stderr, "MetalsStellarMass = %.3e\n", gal->MetalsStellarMass);
  fprintf(stderr, "EjectedGas = %.3e\n", gal->EjectedGas);
  fprintf(stderr, "MetalsEjectedGas = %.3e\n", gal->MetalsEjectedGas);
  fprintf(stderr, "DiskScaleLength = %.3e\n", gal->DiskScaleLength);
  fprintf(stderr, "-------------------\n");
}

//! Evolve existing galaxies forward in time
int evolve_galaxies(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot, int NGal, int NFof)
{

  galaxy_t *gal          = NULL;
  halo_t   *halo         = NULL;
  int       gal_counter  = 0;
  int       dead_gals    = 0;
  double    infalling_gas = 0;
  double    cooling_mass   = 0;

  SID_log("Doing physics...", SID_LOG_OPEN|SID_LOG_TIMER);


  for(int i_fof=0; i_fof<NFof; i_fof++)
  {

    if((fof_group[i_fof].FirstHalo->Galaxy != NULL) && (fof_group[i_fof].FirstHalo->Galaxy->id_MBP == DEBUG_MBP))
    {
      fprintf(stderr, "\nBefore infall calc:\n");
      print_galaxy(fof_group[i_fof].FirstHalo->Galaxy);
    }
    infalling_gas = gas_infall(run_globals, &(fof_group[i_fof]), snapshot);
    if((fof_group[i_fof].FirstHalo->Galaxy != NULL) && (fof_group[i_fof].FirstHalo->Galaxy->id_MBP == DEBUG_MBP))
    {
      fprintf(stderr, "After infall calc:\n");
      print_galaxy(fof_group[i_fof].FirstHalo->Galaxy);
    }

    halo = fof_group[i_fof].FirstHalo;
    while (halo!=NULL) {
      gal = halo->Galaxy;

      if(gal!=NULL)
      {
        if(gal->id_MBP == DEBUG_MBP)
        {
          fprintf(stderr, "Before cooling calc:\n");
          print_galaxy(gal);
        }
        cooling_mass = gas_cooling(run_globals, gal, snapshot);
        if(gal->id_MBP == DEBUG_MBP)
        {
          fprintf(stderr, "After cooling calc:\n");
          print_galaxy(gal);
        }
      }

      while(gal!=NULL){

        if(gal->Type == 0)
        {
          if(gal->id_MBP == DEBUG_MBP)
          {
            fprintf(stderr, "Before add infall:\n");
            print_galaxy(gal);
          }
          add_infall_to_hot(gal, infalling_gas);
          if(gal->id_MBP == DEBUG_MBP)
          {
            fprintf(stderr, "After add infall:\n");
            print_galaxy(gal);
          }
          if(gal->id_MBP == DEBUG_MBP)
          {
            fprintf(stderr, "Before reincorporation:\n");
            print_galaxy(gal);
          }
          reincorporate_ejected_gas(run_globals, gal);
          if(gal->id_MBP == DEBUG_MBP)
          {
            fprintf(stderr, "After reincorporation:\n");
            print_galaxy(gal);
          }
          if(gal->id_MBP == DEBUG_MBP)
          {
            fprintf(stderr, "Before cooling gas:\n");
            print_galaxy(gal);
          }
          cool_gas_onto_galaxy(gal, cooling_mass);
          if(gal->id_MBP == DEBUG_MBP)
          {
            fprintf(stderr, "After cooling gas:\n");
            print_galaxy(gal);
          }
        }

        if(gal->id_MBP == DEBUG_MBP)
        {
          fprintf(stderr, "Before insitu SF:\n");
          print_galaxy(gal);
        }
        insitu_star_formation(run_globals, gal, snapshot);
        if(gal->id_MBP == DEBUG_MBP)
        {
          fprintf(stderr, "After insitu SF:\n");
          print_galaxy(gal);
        }

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
          {
            if(gal->MergerTarget->id_MBP == DEBUG_MBP)
            {
              fprintf(stderr, "Before merging:\n");
              print_galaxy(gal->MergerTarget);
            }
            merge_with_target(run_globals, gal, &dead_gals, snapshot);
            if(gal->MergerTarget->id_MBP == DEBUG_MBP)
            {
              fprintf(stderr, "After merging:\n");
              print_galaxy(gal->MergerTarget);
            }
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
