#include <math.h>
#include "meraxes.h"

void gas_infall(run_globals_t *run_globals, fof_group_t *FOFgroup, int snapshot)
{

  halo_t *halo;
  galaxy_t *gal;
  double total_baryons = 0.;
  double infall_mass = 0.;
  double FOF_Mvir = FOFgroup->FirstHalo->Mvir;
  double fb_modifier;
  halo = FOFgroup->FirstHalo;

  // Calculate the total baryon mass in the FOF group
  halo = FOFgroup->FirstHalo;
  while(halo != NULL)
  {
    gal = halo->Galaxy;
    while(gal != NULL)
    {
      total_baryons += gal->StellarMass + gal->HotGas + gal->ColdGas + gal->EjectedGas;
      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  // Calculate the amount of fresh gas required to provide the baryon
  // fraction of this halo.
  if(run_globals->params.physics.ReionizationModifier)
    fb_modifier = reionization_modifier(run_globals, FOFgroup->FirstHalo, snapshot);
  else
    fb_modifier = 1.0;
  infall_mass = fb_modifier * run_globals->params.BaryonFrac * FOF_Mvir - total_baryons;

  // record the infall modifier
  gal = FOFgroup->FirstHalo->Galaxy;
  gal->BaryonFracModifier = fb_modifier;

  // Give this mass to the central
  if(infall_mass > 1e-10)
    if(gal != NULL)
      gal->HotGas += infall_mass;

}
