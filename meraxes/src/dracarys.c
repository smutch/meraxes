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

}

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{
  trees_header_struct  trees_header;
  halo_struct         *Halo;
  int                  i_newhalo;
  int                  NGal         = 0;
  double               dt;
  run_params_struct    params       = run_globals->params;
  galaxy_struct       *Gal          = NULL;

  for(int snapshot=0; snapshot<MAXSNAPS; snapshot++)
  {
    trees_header = read_halos(run_globals, snapshot, &Halo);

    // If this is the first read then use the n_halos_max parameter of the trees_header to malloc the galaxy array...
    if (Gal==NULL)
      init_galaxies(Gal, trees_header.n_halos_max);
    else
    {
      // otherwise, loop through each existing galaxy and update the properties appropriately.
      for(int i_gal=0; i_gal<NGal; i_gal++)
      {
        i_newhalo = Gal[i_gal].HaloDesc;
        dt = run_globals->Age[snapshot-1]-run_globals->Age[snapshot];

        if(i_newhalo==-1)
        {
          // Here we have a halo where we have lost tracking so we make the corresponding galaxy a type 2
          Gal[i_gal].Type = 2;

          // TODO: Start the merger clock etc.
        } else
        {
          copy_halo_to_galaxy(run_globals, &(Halo[i_newhalo]), &(Gal[i_gal]));

          Halo[i_newhalo].NGalaxies  = Gal[i_gal].HaloNGal;
          Gal[i_gal].CentralMvir     = Gal[Gal[i_gal].CentralGal].Mvir;
          Gal[i_gal].MergTime       -= dt;

          if ((Halo[i_newhalo].Mvir-Gal[i_gal].Mvir) >0.0)
          {
            Gal[i_gal].dM            = Halo[i_newhalo].Mvir-Gal[i_gal].Mvir;
            Gal[i_gal].dMdt          = Gal[i_gal].dM / dt;
          }
        }
      }
    }

    // Create new galaxies in empty type 0 halos
    for(int i_halo=0; i_halo<trees_header.n_subgroups; i_halo++)
    {
      if ( (Halo[i_halo].Type==0) && (Halo[i_halo].NGalaxies==0) )
      {
        copy_halo_to_galaxy(run_globals, &(Halo[i_halo]), &(Gal[NGal]));
        Gal[NGal].CentralGal  = NGal;
        Gal[NGal].CentralMvir = Halo[i_halo].Mvir;

        // Increment galaxy number counters
        Halo[i_halo].NGalaxies++;
        NGal++; 
      }
    }
    
    evolve_galaxies(run_globals, Gal, snapshot, NGal);

    // TODO: Save galaxies if this is an output snapshot
    
    free_halos(&Halo);
  }

  SID_free(SID_FARG Gal);

}


