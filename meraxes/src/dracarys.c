#include "meraxes.h"

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{
  trees_header_struct  trees_header;
  halo_struct         *Halo;
  int                  snapshot;
  int                  i_gal;
  int                  i_newhalo;
  double               dt;
 
  run_params_struct  params = run_globals->params;
  galaxy_struct     *Gal;

  // Initialise galaxy pointers and counters
  Gal = NULL;
  run_globals->NGal = 0;

  for (snapshot=0; snapshot<MAXSNAPS; snapshot++)
  {
    trees_header = read_halos(run_globals, snapshot, &Halo);

    // If this is the first read then use the n_halos_max parameter of the trees_header to malloc the galaxy array
    if (Gal==NULL)
      init_galaxies(Gal, trees_header.n_halos_max);
    else
    {
      // TODO: Copy over progenitor galaxy properties 
      for(i_gal=0; i_gal<run_globals->NGal; i_gal++)
      {
        i_newhalo = Gal[i_gal].HaloDescIndex;
        dt = (Age[snapshot-1]-Age[snapshot]);

        if(i_newhalo==-1)
        {
          // Here we have a halo where we have lost tracking so we make the corresponding galaxy a type 2...
          Gal[i_gal].Type = 2;

          // TODO: Start the merger clock etc.
        } else
        {
          Halo[i_newhalo].NGalaxies = Gal[i_gal].HaloNGal;

          Gal[i_gal].Type = Halo[i_newhalo].Type;
          Gal[i_gal].HaloDescIndex = Halo[i_newhalo].file_index;
          Gal[i_gal].CentralMvir   = Gal[Gal[i_gal].CentralGal].Mvir;
          Gal[i_gal].Mvir = Halo[i_newhalo].Mvir;
          Gal[i_gal].Rvir          = Halo[i_newhalo].Rvir;
          Gal[i_gal].Vvir          = Halo[i_newhalo].Vvir;
          Gal[i_gal].Vmax          = Halo[i_newhalo].Vmax;
          Gal[i_gal].MergTime      -= dt;
          if ((Halo[i_newhalo].Mvir-Gal[i_gal].Mvir) >0.0)
          {
            Gal[i_gal].dM            = Halo[i_newhalo].Mvir-Gal[i_gal].Mvir;
            Gal[i_gal].dMdt          = Gal[i_gal].dM / dt;
          }
        }
      }
    }

    // Create new galaxies in type 0 halos
    
    
    // TODO: Call physics

    // TODO: Save galaxies if this is an output snapshot
    
    free_halos(&Halo);
  }

  SID_free(SID_FARG Gal);

}
