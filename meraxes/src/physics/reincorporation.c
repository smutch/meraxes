#include "meraxes.h"


void reincorporate_ejected_gas(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot)
{

  halo_t *halo;
  galaxy_t *gal;
  galaxy_t *central;

  // loop through each galaxy in the halo and identify all type > 0 ones which
  // have ejected gas.  These galaxies must have just infallen into the FOF
  // group.
  halo = fof_group->FirstHalo;
  central = halo->Galaxy;

  if(central->Type != 0)
  {
    SID_log_error("There is something wrong here... fof_group->FirstHalo->Galaxy.Type != 0 !!!");
    ABORT(EXIT_FAILURE);
  }

  while(halo != NULL)
  {
    gal = halo->Galaxy;
    while(gal != NULL)
    {
      if((gal->Type > 0) && (gal->EjectedGas > 0))
      {
        // This galaxy is in a halo which has just infallen into this FOF group.
        // Give its ejected mass to the central...
        central->HotGas += gal->EjectedGas;
        gal->EjectedGas = 0;
      }
    }
  }

}
