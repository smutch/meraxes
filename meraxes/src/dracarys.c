#include "meraxes.h"

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{
  trees_header_struct  trees_header;
  halo_struct         *Halo;
  int                  snapshot;
  int                  i_newhalo;
  int                  NGal         = 0;
  double               dt;
  run_params_struct    params       = run_globals->params;
  galaxy_struct       *Gal          = NULL;

  for (snapshot=0; snapshot<MAXSNAPS; snapshot++)
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
    
    // TODO: Call physics

    // TODO: Save galaxies if this is an output snapshot
    
    free_halos(&Halo);
  }

  SID_free(SID_FARG Gal);

}
