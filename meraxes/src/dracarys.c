#include "meraxes.h"
#include "tree_flags.h"
#include <math.h>

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
static void evolve_galaxies(run_globals_struct *run_globals, fof_group_struct *fof_group, int snapshot, int NGal)
{

  double         sfr;
  double         BaryonFrac      = run_globals->params.BaryonFrac;
  double         RecycleFraction = run_globals->params.RecycleFraction;
  double         dt              = run_globals->Age[snapshot-1]-run_globals->Age[snapshot];
  galaxy_struct *gal             = NULL;
  galaxy_struct *parent          = NULL;
  halo_struct   *halo            = NULL;
  int            i_fof;
  int            gal_counter     = 0;
  double         mi;
  double         ma;
  double         mass_ratio;
  
  
  for(int i_fof=0; i_fof<NGal; i_fof++)
  {
    halo = fof_group[i_fof].FirstHalo;
    do {
      gal = halo->Galaxy;

      do {
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
      } while(gal!=NULL);

      halo = halo->NextHaloInFOFGroup;
    } while (halo!=NULL);

    // Check for mergers
    gal = fof_group[i_fof].FirstHalo->Galaxy;
    do {
      if(gal->Type == 2)
      {
        if(gal->MergTime <0)
        {
          // Merger!
          parent = gal->MergerTarget;

          // What if the parent (or the parent of the parent!?!) itself merges in this timestep?
          if((parent->Type ==2) && (parent->MergTime <0))
          {
            do{
              parent = parent->MergerTarget;
            } while ((parent->Type ==2) && (parent->MergTime <0));
          }

          // calculate mass ratio of merging galaxies 
          if(gal->StellarMass < parent->StellarMass)
          {
            mi = gal->StellarMass;
            ma = parent->StellarMass;
          }
          else
          {
            mi = parent->StellarMass;
            ma = gal->StellarMass;
          }
          if(ma > 0)
            mass_ratio = mi / ma;
          else
            mass_ratio = 1.0;

          // Add galaxies together
          parent->StellarMass += gal->StellarMass;

          for(int outputbin = 0; outputbin < NOUT; outputbin++)
            parent->Sfr[outputbin] += gal->Sfr[outputbin];

          // Mark the merged galaxy as dead
          gal->Type          = 3;
          gal->HaloDescIndex = -1;

        }
      }
      gal = gal->NextGalInHalo;
    } while(gal!=NULL);
  }
    
  if(gal_counter!=NGal)
  {
    SID_log_error("We have not processed the expected number of galaxies...");
    ABORT(EXIT_FAILURE);
  }

  // TODO: Updating of any final galaxy properties / indices
}


//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{

  trees_header_struct   trees_header;
  halo_struct         **halo;
  halo_struct          *cur_halo;
  fof_group_struct    **fof_group;
  galaxy_struct        *gal;
  galaxy_struct        *prev_gal;
  galaxy_struct        *cur_gal;
  int                   i_newhalo;
  int                   NGal         = 0;
  double                dt;

  for(int snapshot=0; snapshot<MAXSNAPS; snapshot++)
  {

    trees_header = read_halos(run_globals, snapshot, halo, fof_group);

    gal      = run_globals->FirstGal;
    prev_gal = NULL;
    dt       = run_globals->Age[snapshot-1]-run_globals->Age[snapshot];
    
    do {
      i_newhalo = gal->HaloDescIndex;

      if(i_newhalo>-1)
      {
        if( ((gal->TreeFlags & TREE_CASE_MERGER)==TREE_CASE_MERGER)
            && ((gal->TreeFlags & TREE_CASE_MAIN_PROGENITOR)!=TREE_CASE_MAIN_PROGENITOR) )
        {
          // Here we have a merger...  Mark it and deal with it below.
          gal->Type = 999;
          gal->Halo = &(*halo[i_newhalo]);
        } else if(gal->Type < 2)
        {
          copy_halo_to_galaxy(run_globals, &(*halo[i_newhalo]), gal);
          (*halo[i_newhalo]).Galaxy = gal;
        }
      } else
      {
        // This galaxy is done (merged, lost, whatever...) so get rid of it
        if(prev_gal!=NULL)
          prev_gal->Next = gal->Next;
        else
          run_globals->FirstGal = gal->Next;
        SID_free(SID_FARG gal);
        gal = prev_gal;
        NGal--;
      }

      prev_gal = gal;
      gal = gal->Next;
    } while (gal != NULL);

    // Incase we ended up removing the last galaxy, update the LastGal pointer
    run_globals->LastGal = prev_gal;

    // Find empty type 0 halos and place new galaxies in them
    for(int i_halo=0; i_halo<trees_header.n_subgroups; i_halo++)
    {
      if(((*halo[i_halo]).Type == 0) && ((*halo[i_halo]).Galaxy == NULL))
      {
        gal = SID_malloc(sizeof(galaxy_struct));
        copy_halo_to_galaxy(run_globals, &(*halo[i_halo]), gal);
        run_globals->LastGal->Next = gal;
        run_globals->LastGal = gal;
        NGal++;
      }
    }

    // Loop through each galaxy and deal with mergers now that all other galaxies have been 
    // correctly propogated forwards
    gal = run_globals->FirstGal;
    do {
      if(gal->Type == 999)
      {
        gal->Type = 2;
        cur_gal = gal->Halo->Galaxy;
        do {
          prev_gal = cur_gal;
          cur_gal = cur_gal->NextGalInHalo;
        } while (cur_gal!=NULL);
        prev_gal->NextGalInHalo = gal;
        
        gal->MergTime = calculate_merging_time(run_globals, gal, snapshot);
      }
      gal = gal->Next;
    } while (gal!=NULL);

    evolve_galaxies(run_globals, *fof_group, snapshot, NGal);
  
    SID_free(SID_FARG *halo);
    SID_free(SID_FARG *fof_group);
  }


}


