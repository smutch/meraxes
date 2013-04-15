#include "meraxes.h"

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{
  trees_header_struct  trees_header;
  halo_struct         *halo;
  int                  snapshot;
 
  run_params_struct  params = run_globals->params;
  galaxy_struct     *Gal;

  // Initialise galaxy pointers and counters
  Gal = NULL;
  run_globals->Ngal = 0;

  for (snapshot=0; snapshot<MAXSNAPS; snapshot++)
  {
    trees_header = read_halos(run_globals, snapshot, &halo);

    // If this is the first read then use the n_halos_max parameter of the trees_header to malloc the galaxy array
    if (Gal==NULL)
      init_galaxies(Gal, trees_header.n_halos_max);
    else
    {
      // TODO: Copy over current halo properties to all galaxies
    }

    // TODO: Create new galaxies in type 0 halos
    
    // TODO: Call physics

    // TODO: Save galaxies if this is an output snapshot
    
    free_halos(&halo);
  }

  SID_free(SID_FARG Gal);

}
