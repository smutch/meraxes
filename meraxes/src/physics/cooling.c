#include "meraxes.h"
#include <math.h>

double gas_cooling(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{

  double cooling_mass;

  // we only need to do cooling if there is anything to cool!
  if(gal->HotGas > 1e-10)
  {

    double t_cool, max_cooling_mass, cooling_metals, Tvir;
    double logZ, lambda, x, rho_r_cool, r_cool, rho_at_Rvir;
    run_units_t *units = &(run_globals->units);

    // following Croton+ 2006, we set the maximum cooling time to be the
    // dynamical time of the host dark matter halo
    t_cool = gal->Rvir / gal->Vvir;  // internal units

    // calculate the halo virial temperature
    Tvir = 35.9 * gal->Vvir * gal->Vvir;  // internal units (Kelvin)

    // get the log10(metallicity) value
    if(gal->MetalsHotGas > 0)
      logZ = log10(calc_metallicity(gal->HotGas, gal->MetalsHotGas));
    else
      logZ = -10.0;

    // interpolate the temperature and metallicity dependant cooling rate (lambda)
    if(gal->id_MBP == DEBUG_MBP)
      lambda = interpolate_cooling_rate(log10(Tvir), logZ, 1);
    else
      lambda = interpolate_cooling_rate(log10(Tvir), logZ, 0);

    // following equation (3) of Croton+ 2006, calculate the hot gas density at
    // the radius r_cool (i.e. where the cooling time is equal to `t_cool`
    // above)
    x = PROTONMASS * BOLTZMANN * Tvir / lambda;              // now this has units sec g/cm^3
    x /= (units->UnitDensity_in_cgs * units->UnitTime_in_s); // now in internal units
    rho_r_cool = x / t_cool * 0.885;                         // 0.885 = 3/2 * mu, mu=0.59 for a fully ionized gas

    // under the assumption of an isothermal density profile extending to Rvir,
    // now calculate the cooling radius
    rho_at_Rvir = gal->HotGas / (4. * M_PI * gal->Rvir);
    r_cool = sqrt(rho_at_Rvir / rho_r_cool);

    // the maximum amount of gas we can possibly cool is limited by the amount
    // of mass within the free fall radius
    max_cooling_mass = gal->HotGas / t_cool * gal->dt;

    if(r_cool > gal->Rvir)
      // here we are in the rapid cooling regime and we accrete all gas within
      // the free-fall radius
      // cooling_mass = max_cooling_mass;
      cooling_mass = gal->HotGas;
    else
    {
      // here we are in the hot halo regime (but still limited by what's inside the free-fall radius)
      cooling_mass = 0.5 * gal->HotGas / gal->Rvir * r_cool / t_cool * gal->dt;
      // if(cooling_mass > max_cooling_mass)
      //   cooling_mass = max_cooling_mass;
    }

    // DEBUG
    if(gal->id_MBP == DEBUG_MBP)
    {
      fprintf(stderr, "COOLING DEBUG: (%d)\n", snapshot);
      fprintf(stderr, "r_cool = %.3e\n", r_cool);
      fprintf(stderr, "rho_at_Rvir = %.3e\n", rho_at_Rvir);
      fprintf(stderr, "HotGas = %.3e\n", gal->HotGas);
      fprintf(stderr, "Rvir = %.3e\n", gal->Rvir);
      fprintf(stderr, "Tvir = %.3e\n", Tvir);
      fprintf(stderr, "lambda = %.3e\n", lambda);
      fprintf(stderr, "x = %.3e\n", x);
      fprintf(stderr, "logZ = %.3e\n", logZ);
      fprintf(stderr, "t_cool = %.3e\n", t_cool);
      fprintf(stderr, "dt = %.3e\n", gal->dt);
      fprintf(stderr, "cooling_mass = %.3e\n", cooling_mass);
    }

    // do one last sanity check to ensure we aren't cooling more gas than is available etc.
    if(cooling_mass > gal->HotGas)
      cooling_mass = gal->HotGas;
    if(cooling_mass < 0)
      cooling_mass = 0.0;

  }
  else  // if there is no gas to cool...
    cooling_mass = 0.0;

  return cooling_mass;

}


void cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass)
{

  double cooling_metals;

  if(cooling_mass > gal->HotGas)
    cooling_mass = gal->HotGas;

  // what mass of metals is coming along with this cooling gas?
  cooling_metals = cooling_mass * calc_metallicity(gal->HotGas, gal->MetalsHotGas);

  // debug("%d %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
  //     gal->Type, gal->dt, gal->Rvir, gal->Vvir, gal->HotGas, gal->MetalsHotGas, t_cool, logZ, Tvir,
  //     lambda, rho_r_cool, rho_at_Rvir, r_cool, max_cooling_mass, cooling_mass);

  // save the cooling mass
  gal->Mcool = cooling_mass;

  // update the galaxy reservoirs
  gal->HotGas -= cooling_mass;
  gal->MetalsHotGas -= cooling_metals;
  gal->ColdGas += cooling_mass;
  gal->MetalsColdGas += cooling_metals;
}
