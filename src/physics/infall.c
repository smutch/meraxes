#include <math.h>

#include "core/misc_tools.h"
#include "core/modifiers.h"
#include "infall.h"
#include "meraxes.h"
#include "reionization.h"
#include "tree_flags.h"

double gas_infall(fof_group_t* FOFgroup, int snapshot)
{
  halo_t* halo = FOFgroup->FirstHalo;

  galaxy_t* gal;
  galaxy_t* central;
  double total_baryons = 0.;
  double infall_mass = 0.;
  double FOF_Mvir = FOFgroup->Mvir;
  double FOFMvirModifier = FOFgroup->FOFMvirModifier;
  double fb_modifier;

  double total_stellarmass = 0.0;
  double total_hotgas = 0.0;
  double total_coldgas = 0.0;
  double total_ejectedgas = 0.0;
  double total_blackholemass = 0.0;
#if USE_MINI_HALOS
  double total_remnantmass = 0.0; // This come either from the BHs formed after Pop III stars that fail becoming SN and
                                  // directly collapse or from CCSN. Atm these don't accrete and they don't do anything.
#endif

  // Calculate the total baryon mass in the FOF group
  central = FOFgroup->FirstOccupiedHalo->Galaxy;

  while (halo != NULL) {
    gal = halo->Galaxy;
    while (gal != NULL) {
      total_stellarmass += gal->StellarMass;
      total_hotgas += gal->HotGas;
      total_coldgas += gal->ColdGas;
      total_ejectedgas += gal->EjectedGas;
      total_blackholemass += gal->BlackHoleMass + gal->BlackHoleAccretingColdMass;
#if USE_MINI_HALOS
      total_remnantmass += gal->Remnant_Mass;
#endif

      if (gal != central) {
        central->HotGas += gal->HotGas + gal->EjectedGas;
        central->MetalsHotGas += gal->MetalsHotGas + gal->MetalsEjectedGas;
        gal->HotGas = 0.0;
        gal->MetalsHotGas = 0.0;
        gal->EjectedGas = 0.0;
        gal->MetalsEjectedGas = 0.0;
      }

      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  if (check_for_flag(TREE_CASE_BELOW_VIRIAL_THRESHOLD, FOFgroup->FirstHalo->TreeFlags)) {
    // no infall and no hydrstatic hot halo in this case
    central->BaryonFracModifier = 0.0;
    central->EjectedGas += central->HotGas;
    central->MetalsEjectedGas += central->MetalsHotGas;
    central->HotGas = 0.0;
    central->MetalsHotGas = 0.0;
    return 0.0;
  }

  total_baryons = total_stellarmass + total_hotgas + total_coldgas + total_ejectedgas + total_blackholemass;
#if USE_MINI_HALOS
  total_baryons += total_remnantmass;
#endif

  // Calculate the amount of fresh gas required to provide the baryon
  // fraction of this halo.
  fb_modifier = reionization_modifier(central, FOF_Mvir, snapshot);
  if (run_globals.RequestedBaryonFracModifier == 1)
    fb_modifier *= interpolate_modifier(run_globals.baryon_frac_modifier,
                                        log10(FOF_Mvir / FOFMvirModifier / run_globals.params.Hubble_h) + 10.0);
  infall_mass = fb_modifier * run_globals.params.BaryonFrac * FOF_Mvir - total_baryons;

  // record the infall modifier
  central->BaryonFracModifier = fb_modifier;

  return infall_mass;
}

void add_infall_to_hot(galaxy_t* central, double infall_mass)
{
#if USE_MINI_HALOS
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);
#endif
  // if we have mass to add then give it to the central
  if (infall_mass > 0) {
    central->HotGas += infall_mass;
#if USE_MINI_HALOS
    if (Flag_Metals == true) {
      if (central->Flag_ExtMetEnr == 1) // If the halo is externally enriched, it will accrete polluted gas (metals).
        central->MetalsHotGas += infall_mass * central->Metallicity_IGM;
    }
#endif
  } else {
    double strip_mass = -infall_mass;
    // otherwise, strip the mass from the ejected
    if (central->EjectedGas > 0) {
      double metallicity = calc_metallicity(central->EjectedGas, central->MetalsEjectedGas);
      central->EjectedGas -= strip_mass;
      if (central->EjectedGas < 0) {
        strip_mass = -central->EjectedGas;
        central->EjectedGas = 0.0;
        central->MetalsEjectedGas = 0.0;
      } else {
        central->MetalsEjectedGas -= metallicity * strip_mass;
        strip_mass = 0.0;
      }
    }

    // if we still have mass left to remove after exhausting the mass of
    // the ejected component, the remove as much as we can from the hot gas
    if (strip_mass > 0) {
      double metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);
      central->HotGas -= strip_mass;
      if (central->HotGas < 0) {
        central->HotGas = 0.0;
        central->MetalsHotGas = 0.0;
      } else
        central->MetalsHotGas -= metallicity * strip_mass;
    }
  }
}
