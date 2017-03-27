#include "meraxes.h"

static void update_reservoirs_from_reincorporation(galaxy_t *gal, double reincorporated)
{
  double metals = reincorporated * calc_metallicity(gal->EjectedGas, gal->MetalsEjectedGas);

  gal->EjectedGas       -= reincorporated;
  gal->MetalsEjectedGas -= metals;
  gal->HotGas           += reincorporated;
  gal->MetalsHotGas     += metals;
}


void reincorporate_ejected_gas(galaxy_t *gal)
{
  if (gal->EjectedGas > 0)
  {
    double       t_dyn;
    double       reincorporated;
    double       ReincorporationEff = run_globals.params.physics.ReincorporationEff;
    fof_group_t *fof_group          = gal->Halo->FOFGroup;

    // allow some of the ejected gas associated with the central to be
    // reincorporated following the prescription of Guo 2010 (which is actually
    // almost identical to SAGE).
    t_dyn          = fof_group->Rvir / fof_group->Vvir;
    // reincorporated = ReincorporationEff * gal->Vvir / 220.0 * gal->EjectedGas * (gal->dt / t_dyn);
    reincorporated = ReincorporationEff * gal->EjectedGas * (gal->dt / t_dyn);

    // ensure consistency
    if (reincorporated > gal->EjectedGas)
      reincorporated = gal->EjectedGas;

    // update the baryonic reservoirs
    update_reservoirs_from_reincorporation(gal, reincorporated);
  }
}
