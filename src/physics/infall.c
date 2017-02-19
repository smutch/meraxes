#include <math.h>
#include "meraxes.h"

double gas_infall(fof_group_t *FOFgroup, int snapshot)
{
  galaxy_t *gal;
  galaxy_t *central;
  halo_t   *halo;
  double    total_baryons = 0.;
  double    infall_mass = 0.;
  double    FOF_Mvir = FOFgroup->Mvir;
  double    fb_modifier, fb_modifier_hydro;

  double    total_stellarmass   = 0.0;
  double    total_hotgas        = 0.0;
  double    total_coldgas       = 0.0;
  double    total_ejectedgas    = 0.0;
  double    total_blackholemass = 0.0;

  // Calculate the total baryon mass in the FOF group
  halo    = FOFgroup->FirstHalo;
  central = FOFgroup->FirstOccupiedHalo->Galaxy;

  while (halo != NULL)
  {
    gal = halo->Galaxy;
    while (gal != NULL)
    {
      total_stellarmass   += gal->StellarMass;
      total_hotgas        += gal->HotGas;
      total_coldgas       += gal->ColdGas;
      total_ejectedgas    += gal->EjectedGas;
      total_blackholemass += gal->BlackHoleMass;

      if (gal != central)
      {
        central->HotGas       += gal->HotGas + gal->EjectedGas;
        central->MetalsHotGas += gal->MetalsHotGas + gal->MetalsEjectedGas;
        gal->HotGas            = 0.0;
        gal->MetalsHotGas      = 0.0;
        gal->EjectedGas        = 0.0;
        gal->MetalsEjectedGas  = 0.0;
      }

      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  total_baryons = total_stellarmass + total_hotgas + total_coldgas + total_ejectedgas + total_blackholemass;

  // Calculate the amount of fresh gas required to provide the baryon
  // fraction of this halo.
  fb_modifier       = reionization_modifier(central, FOF_Mvir, snapshot);
  fb_modifier_hydro = 1.0;
  if (run_globals.RequestedBaryonFracModifier == 1)
  {
    double logM;
    logM              = log10(FOF_Mvir * 1e10 / run_globals.params.Hubble_h);
    fb_modifier_hydro = interpolate_modifier(run_globals.baryon_frac_modifier, logM);
    //SID_log("%f\t%f\t%f", SID_LOG_COMMENT, FOF_Mvir, logM, fb_modifier_hydro);
  }
  infall_mass                 = fb_modifier * fb_modifier_hydro * run_globals.params.BaryonFrac * FOF_Mvir - total_baryons;

  // record the infall modifier
  central->BaryonFracModifier = fb_modifier;

  return infall_mass;
}


void add_infall_to_hot(galaxy_t *central, double infall_mass)
{
  // if we have mass to add then give it to the central
  if (infall_mass > 0)
    central->HotGas += infall_mass;
  else
  {
    double strip_mass = -infall_mass;
    // otherwise, strip the mass from the ejected
    if(central->EjectedGas > 0)
    {
      central->EjectedGas -= strip_mass;
      if(central->EjectedGas < 0)
      {
        strip_mass                = -central->EjectedGas;
        central->EjectedGas       = 0.0;
        central->MetalsEjectedGas = 0.0;
      }
      else
      {
        central->MetalsEjectedGas -= calc_metallicity(central->EjectedGas, central->MetalsEjectedGas) * strip_mass;
        strip_mass                 = 0.0;
      }
    }

    // if we still have mass left to remove after exhausting the mass of
    // the ejected component, the remove as much as we can from the hot gas
    if (strip_mass > 0)
    {
      central->HotGas -= strip_mass;
      if (central->HotGas < 0)
      {
        central->HotGas       = 0.0;
        central->MetalsHotGas = 0.0;
      }
      else
        central->MetalsHotGas -= calc_metallicity(central->HotGas, central->MetalsHotGas) * strip_mass;
    }
  }
}