#include "meraxes.h"

static void update_reservoirs_from_reincorporation(galaxy_t *gal, double reincorporated)
{

  double metals = reincorporated * calc_metallicity(gal->EjectedGas, gal->MetalsEjectedGas);

  gal->EjectedGas -= reincorporated;
  gal->MetalsEjectedGas -= metals;
  gal->HotGas += reincorporated;
  gal->MetalsHotGas += metals;

}


void reincorporate_ejected_gas(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot)
{

  halo_t *halo;
  galaxy_t *gal;
  galaxy_t *central;
  double t_dyn;
  double reincorporated;

  double ReincorporationEff = run_globals->params.physics.ReincorporationEff;

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
    if(gal == central)
      gal = gal->NextGalInHalo;
    while(gal != NULL)
    {
      if((gal->Type > 0) && (gal->EjectedGas > 0))
      {
        // This galaxy is in a halo which has just infallen into this FOF group.
        // Give its ejected mass to the central...
        central->HotGas       += gal->EjectedGas;
        central->MetalsHotGas += gal->MetalsEjectedGas;
        gal->EjectedGas        = 0;
        gal->MetalsEjectedGas  = 0;
      }
      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  if(central->EjectedGas > 0)
  {
    // now allow some of the ejected gas associated with the central to be
    // reincorporated following the prescription of Guo 2010 (which is actually
    // almost identical to SAGE).
    t_dyn = central->Rvir / central->Vvir;
    reincorporated = ReincorporationEff * central->Vvir / 220.0 * central->EjectedGas * (central->dt / t_dyn);

    // update the baryonic reservoirs
    update_reservoirs_from_reincorporation(central, reincorporated);
  }

}
