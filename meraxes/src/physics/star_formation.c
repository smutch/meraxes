#include <math.h>
#include "meraxes.h"
#include <gsl/gsl_sf_lambert.h>

void update_reservoirs_from_sf(run_globals_t *run_globals, galaxy_t *gal, double new_stars, double m_recycled)
{
  double metallicity;
  double current_time;
  double remaining_stars;
  double new_metals;

  // update the galaxy's SFR value
  gal->Sfr += new_stars / gal->dt;

  // update the stellar mass history
  gal->NewStars[0] += new_stars;

  // instantaneous recycling approximation of stellar mass
  metallicity     = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);
  remaining_stars = new_stars - m_recycled;

  gal->ColdGas           -= remaining_stars;
  gal->MetalsColdGas     -= remaining_stars * metallicity;
  gal->StellarMass       += remaining_stars;
  gal->MetalsStellarMass += remaining_stars * metallicity;

  // update the luminosities
  current_time = gal->LTTime + 0.5 * gal->dt;
  add_to_luminosities(run_globals, gal, new_stars, metallicity, current_time);
}


void insitu_star_formation(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{
  // These star formation & supernova feedback prescriptions are taken directly
  // from Croton+ 2006

  // there is no point doing anything if there is no cold gas!
  if (gal->ColdGas > 1e-10)
  {
    double r_disk;
    double m_crit;
    double m_stars;
    double m_reheat;
    double m_eject;
    double new_metals;
    double m_recycled;

    double SfEfficiency = run_globals->params.physics.SfEfficiency;

    // calculate disk scalelength using Mo, Mau & White (1998) eqn. 12 and
    // multiply it by 3 to approximate the star forming region size (ala
    // Croton+ 2006).
    r_disk = gal->DiskScaleLength * 3.0;

    // what is the critical mass within r_crit?
    m_crit = 0.19 * gal->Vmax * r_disk;

    if (gal->ColdGas > m_crit)
      m_stars = SfEfficiency * (gal->ColdGas - m_crit) / r_disk * gal->Vmax * gal->dt;
    else
      // no star formation
      return;

    // apply supernova feedback and update baryonic reservoirs
    supernova_feedback(run_globals, gal, &m_stars, &m_reheat, &m_eject, &m_recycled, &new_metals, snapshot);

    // update the baryonic reservoirs (note the order makes a difference here!)
    update_reservoirs_from_sf(run_globals, gal, m_stars, m_recycled);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, new_metals);
  }
}
