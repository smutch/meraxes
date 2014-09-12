#include "meraxes.h"
#include <math.h>

double radio_mode_BH_heating(run_globals_t *run_globals, galaxy_t *gal, double cooling_mass)
{
  double accretion_rate;
  double eddington_rate;
  double accreted_mass;
  double heated_mass;
  double metallicity;

  fof_group_t *fof_group = gal->Halo->FOFGroup;

  run_units_t *units = &(run_globals->units);

  // if there is any hot gas
  if (gal->HotGas > 0.0)
  {
    // empirical accretion recipe of Croton et al. (2006)
    accretion_rate = run_globals->params.physics.RadioModeEff
                     / (units->UnitMass_in_g / units->UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
                     * (gal->BlackHoleMass / 0.01) * pow(fof_group->Vvir / 200.0, 3.0)
                     * ((gal->HotGas / fof_group->Mvir) / 0.1);

    // Eddington rate
    eddington_rate = 1.3e48 * gal->BlackHoleMass / (units->UnitEnergy_in_cgs / units->UnitTime_in_s) / 9e10;

    // limit accretion by the eddington rate
    if (accretion_rate > eddington_rate)
      accretion_rate = eddington_rate;

    accreted_mass = accretion_rate * gal->dt;

    // limit accretion by amount of hot gas available
    if (accreted_mass > gal->HotGas)
      accreted_mass = gal->HotGas;

    // mass heated by AGN following Croton et al. 2006
    // 1.34e5 = sqrt(2*eta*c^2), eta=0.1 (standard efficiency) and c in km/s
    heated_mass = (1.34e5 / fof_group->Vvir) * (1.34e5 / fof_group->Vvir) * accreted_mass;

    // limit the amount of heating to the amount of cooling
    if (heated_mass > cooling_mass)
    {
      accreted_mass = cooling_mass / heated_mass * accreted_mass;
      heated_mass   = cooling_mass;
    }

    // add the accreted mass to the black hole
    metallicity         = calc_metallicity(gal->HotGas, gal->MetalsHotGas);
    gal->BlackHoleMass += accreted_mass;
    gal->HotGas        -= accreted_mass;
    gal->MetalsHotGas  -= metallicity * accreted_mass;
  }
  else
  {
    // if there is no hot gas
    heated_mass = 0.0;
  }

  return heated_mass;
}


void merger_driven_BH_growth(run_globals_t *run_globals, galaxy_t *gal, double merger_ratio)
{
  if (gal->ColdGas > 0)
  {
    // If there is any cold gas to feed the black hole...
    double accreted_mass;
    double accreted_metals;
    double Vvir;

    // If this galaxy is the central of it's FOF group then use the FOF halo properties
    // TODO: This needs closer thought as to if this is the best thing to do...
    if (gal->Type == 0)
      Vvir = gal->Halo->FOFGroup->Vvir;
    else
      Vvir = gal->Vvir;

    accreted_mass = run_globals->params.physics.BlackHoleGrowthRate * merger_ratio /
                    (1.0 + (280.0 * 280.0 / Vvir / Vvir)) * gal->ColdGas;

    // limit accretion to what is available
    if (accreted_mass > gal->ColdGas)
      accreted_mass = gal->ColdGas;

    accreted_metals     = calc_metallicity(gal->ColdGas, gal->MetalsColdGas) * accreted_mass;
    gal->BlackHoleMass += accreted_mass;
    gal->ColdGas       -= accreted_mass;
    gal->MetalsColdGas -= accreted_metals;
  }
}
