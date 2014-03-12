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
    double r_disk;
    double m_crit;
    double m_stars;

    double SfEfficiency = run_globals->params.physics.SfEfficiency;
    double SfRecycleFraction = run_globals->params.physics.SfRecycleFraction;
    double sqrt_2 = 1.414213562;

    // calculate disk scalelength using Mo, Mau & White (1998) eqn. 12 and
    // multiply it by 3 to approximate the star forming region size (ala
    // Croton+ 2006).
    r_disk = gal->Spin * gal->Rvir / sqrt_2 * 3.0;

    // what is the critical mass within r_crit?
    m_crit = 0.19 * gal->Vvir * r_disk;

    if(gal->ColdGas > m_crit)
      gal->Sfr = SfEfficiency * (gal->ColdGas - m_crit) / r_disk * gal->Vvir;
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
