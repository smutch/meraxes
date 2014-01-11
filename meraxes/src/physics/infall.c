#include <math.h>
#include "meraxes.h"

void gas_infall(run_globals_t *run_globals, fof_group_t *FOFgroup)
{

  halo_t *halo;
  galaxy_t *gal;
  double total_baryons = 0.;
  double infall_mass = 0.;
  double FOF_Mvir = FOFgroup->FirstHalo->Mvir;
  double used_mass = 0.;
  double mass = 0.;

  // Calculate the total baryon mass in the FOF group
  halo = FOFgroup->FirstHalo;
  while(halo != NULL)
  {
    gal = halo->Galaxy;
    while(gal != NULL)
    {
      total_baryons += gal->StellarMass + gal->Gas;
      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  // Calculate the amount of fresh gas required to provide the universal baryon
  // fraction.
  infall_mass = run_globals->params.BaryonFrac * FOF_Mvir - total_baryons;

  // Split the infalling gas up umongst the subhalos with galaxies in them
  halo = FOFgroup->FirstHalo->NextHaloInFOFGroup;
  while(halo != NULL)
  {
    gal = halo->Galaxy;
    if(gal != NULL)
    {
      mass = (gal->Mvir/FOF_Mvir) * infall_mass;
      gal->Gas += mass;
      used_mass += mass;
      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  // The remaining gas goes to the central galaxy
  gal = FOFgroup->FirstHalo->Galaxy;
  gal->Gas += infall_mass-used_mass;

}
