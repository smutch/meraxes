#include <math.h>
#include "meraxes.h"
#include <gsl/gsl_sf_lambert.h>

void update_reservoirs_from_sf(run_globals_t *run_globals, galaxy_t *gal, double new_stars)
{
  double metals;

  gal->Sfr = new_stars / gal->dt;

  // instantaneous recycling approximation
  new_stars = (1.0 - run_globals->params.physics.SfRecycleFraction) * new_stars;

  metals = calc_metallicity(gal->ColdGas, gal->MetalsColdGas) * new_stars;

  gal->ColdGas           -= new_stars;
  gal->MetalsColdGas     -= metals;
  gal->StellarMass       += new_stars;
  gal->MetalsStellarMass += metals;
}


void insitu_star_formation(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{

  // These star formation & supernova feedback prescriptions are taken directly
  // from Croton+ 2006

  // there is no point doing anything if there is no cold gas!
  if(gal->ColdGas > 1e-8)
  {
    double r_disk;
    double m_crit;
    double m_stars;

    double SfEfficiency = run_globals->params.physics.SfEfficiency;
    double sqrt_2 = 1.414213562;

    // calculate disk scalelength using Mo, Mau & White (1998) eqn. 12 and
    // multiply it by 3 to approximate the star forming region size (ala
    // Croton+ 2006).
    r_disk = gal->Spin * gal->Rvir / sqrt_2 * 3.0;

    // what is the critical mass within r_crit?
    m_crit = 0.19 * gal->Vvir * r_disk;

    if(gal->ColdGas > m_crit)
      m_stars = SfEfficiency * (gal->ColdGas - m_crit) / r_disk * gal->Vmax * gal->dt;
    else
      // no star formation
      return;

    // apply supernova feedback and update baryonic reservoirs
    supernova_feedback(run_globals, gal, m_stars);
  }
}
