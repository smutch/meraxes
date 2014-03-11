#include <math.h>
#include "meraxes.h"
#include <gsl/gsl_sf_lambert.h>

static void update_reservoirs_from_sf(galaxy_t *gal, double new_stars)
{
  double metals = calc_metallicity(gal->ColdGas, gal->MetalsColdGas) * new_stars;
  gal->ColdGas -= new_stars;
  gal->MetalsColdGas -= metals;
  gal->StellarMass += new_stars;
  gal->MetalsStellarMass += metals;
}

void form_stars_insitu(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{

  // This is a new SF prescription which is based on Croton+ 2006, but is more
  // self consistent with the idea of a critical surface density threshold for
  // star formation.

  // there is no point doing anything if there is no cold gas!
  if(gal->ColdGas > 1e-8)
  {
    double r_d;  // disk scale radius
    double r_crit; // radius at which gas surface density exceeds critical for SF
    double total_sd;  // gas total surface density
    double m_crit;
    double m_gas;
    double b;
    double m_stars;
    double r_frac;

    double SfEfficiency = run_globals->params.physics.SfEfficiency;
    double SfRecycleFraction = run_globals->params.physics.SfRecycleFraction;
    double sqrt_2 = 1.414213562;

    // calculate disk scalelength using Mo, Mau & White (1998) eqn. 12
    r_d = gal->Spin * gal->Rvir / sqrt_2 * 1e3;  // kpc/h

    // calculate the total gas surface density
    total_sd = gal->ColdGas * 1e10 / (2.0 * M_PI * r_d * r_d); // Msol/h (pc/h)^-2

    // Assuming a critial surface density for star formation using the halo
    // dependant approximation of Kauffmann 1996 (eqn. 7), calculate the radius
    // out to which the cold gas surface density exceeds the critical gas surface
    // density.
    b = gal->Vvir * 0.59 / total_sd;
    debug("%d %d %e %e %e %e\n", snapshot, gal->ID, b, r_d, gal->ColdGas, gal->Vvir);
    r_crit = -r_d * gsl_sf_lambert_Wm1(-b/r_d);
    if(r_crit < 0)
    {
      gal->Sfr = 0.0;
      return;
    }

    // what is the critical mass within r_crit?
    m_crit = 0.59 * 2.0 * M_PI * gal->Vvir * r_crit / 1.e10; // 1e10 Msol/h

    // what is the cold gas mass within r_crit?
    r_frac = r_crit/r_d;
    m_gas = total_sd * 2.0 * M_PI * r_d*r_d * ( 1.0 - exp(-r_frac)*(r_frac + 1.0) ) / 1.e10;  // 1e10 Msol/h

    // now use a Croton+ 2006 style SF law to determine the SFR
    if(m_gas > m_crit)
      gal->Sfr = SfEfficiency * (m_gas - m_crit) / gal->Rvir * gal->Vvir;
    else
    {
      gal->Sfr = 0.0;
      return;
    }

    // make sure we aren't trying to create more stars than we have cold gas...
    m_stars = gal->Sfr * gal->dt;
    if (m_stars > gal->ColdGas)
    {
      gal->Sfr = (gal->ColdGas / m_stars) * gal->Sfr;
      m_stars = gal->ColdGas;
    }

    // instantaneous recycling approximation
    m_stars = (1 - SfRecycleFraction) * m_stars;

    // update the galaxy reservoirs
    update_reservoirs_from_sf(gal, m_stars);

  }
}
