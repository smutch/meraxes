#include "meraxes.h"
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
static void evolve_galaxies(run_globals_struct *run_globals, galaxy_struct *Gal, int snapshot, int NGal)
{

  double sfr;
  double BaryonFrac = run_globals->params.BaryonFrac;
  double RecycleFraction = run_globals->params.RecycleFraction;
  double dt = run_globals->Age[snapshot-1]-run_globals->Age[snapshot];
  
  for(int i_gal=0; i_gal<NGal; i_gal++)
  {
    if((Gal[i_gal].Mvir>0.0) && (Gal[i_gal].Type==0))
      switch (run_globals->params.physics.funcprop){
        case VMAX_PROP:
          sfr = BaryonFrac*Gal[i_gal].dMdt * physics_func(run_globals, Gal[i_gal].Vmax, snapshot);
          break;
        case MVIR_PROP:
          sfr = BaryonFrac*Gal[i_gal].dMdt * physics_func(run_globals, Gal[i_gal].Mvir*1.e10, snapshot);
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
        Gal[i_gal].Sfr[outputbin] += sfr;
        break;
      }
    }

    // Instantaneous recycling approximation
    Gal[i_gal].StellarMass += (1.0-RecycleFraction)*sfr*dt;
  }

  // Check for mergers
  int i_central;
  double mi, ma, mass_ratio;
  for(int i_gal=0; i_gal<NGal; i_gal++)
  {
    if(Gal[i_gal].Type == 2)
    {
      Gal[gal_id].MergTime -= dt;

      if(Gal[gal_id].MergTime <0)
      {
        // Merger baby!
        i_central = Gal[i_gal].CentralGal;
        
        // calculate mass ratio of merging galaxies 
        if(Gal[i_gal].StellarMass < Gal[i_central].StellarMass)
        {
          mi = Gal[i_gal].StellarMass;
          ma = Gal[i_central].StellarMass;
        }
        else
        {
          mi = Gal[i_central].StellarMass;
          ma = Gal[i_gal].StellarMass;
        }
        if(ma > 0)
          mass_ratio = mi / ma;
        else
          mass_ratio = 1.0;

        // Add galaxies together
        Gal[i_gal].Type = 3;
        Gal[i_central].StellarMass += Gal[i_gal].StellarMass;
        Gal[i_central].BulgeMass += Gal[i_gal].StellarMass;

        for(int outputbin = 0; outputbin < NOUT; outputbin++)
          Gal[i_central].Sfr[outputbin] += Gal[i_gal].Sfr[outputbin];

      }
    }
  }

  // TODO: Updating of any final galaxy properties / indices
}


//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{

  trees_header_struct trees_header;
  halo_struct **halo;
  fof_group_struct **fof_group
  galaxy_struct *gal, *prev_gal;
  int i_newhalo;
  double dt;

  for(int snapshot=0; snapshot<MAXSNAPS; snapshot++)
  {

    trees_header = read_halos(run_globals, snapshot, halo, fof_group);

    gal = run_globals->FirstGal;
    prev_gal = NULL;
    dt = run_globals->Age[snapshot-1]-run_globals->Age[snapshot];
    
    do {
      i_newhalo = gal->HaloDescIndex;

      if(i_new_halo>-1)
      {
        if( ((gal->TreeFlags & TREE_CASE_MERGER)==TREE_CASE_MERGER)
            && ((gal->TreeFlags & TREE_CASE_MAIN_PROGENITOR)!=TREE_CASE_MAIN_PROGENITOR) )
        {
          // Here we have a merger...  Mark it and deal with it below.
          gal->Type = 2;
        } else
        {
          copy_halo_to_galaxy(run_globals, &(halo[i_new_halo]), gal);
          halo[i_new_halo].Galaxy = gal;
          if(halo[i_new_halo].Type == 0)
            halo[i_new_halo].CentralGal = gal;
        }
      } else
      {
        // This galaxy is done (merged, lost, whatever...) so get rid of it
        if(prev_gal!=NULL)
          prev_gal.Next = gal.Next;
        else
          run_globals->FirstGal = gal.Next;
        SID_free(SID_FARG gal);
      }

      prev_gal = gal;
      gal = gal.Next;
    } while (gal != NULL);

    // Incase we ended up removing the last galaxy, update the LastGal pointer
    run_globals->LastGal = prev_gal;

    // Find empty type 0 halos and place new galaxies in them
    for(int i_halo=0; i_halo<trees_header.n_subgroups; i_halo++)
    {
      if((*halo[i_halo].Type == 0) && (*halo[i_halo].Galaxy == NULL))
      {
        gal = SID_malloc(sizeof(galaxy_struct));
        copy_halo_to_galaxy(run_globals, &(halo[i_halo]), gal);
        run_globals->LastGal->Next = gal;
        run_globals->LastGal = gal;
      }
    }

    // Loop through all galaxies again and connect the FOF group members and deal with mergers
  
  SID_free(SID_FARG *halo);
  SID_free(SID_FARG *fof_group);
  }

}


