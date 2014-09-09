#include <math.h>
#include "meraxes.h"
#include <gsl/gsl_sf_lambert.h>
#include <assert.h>

void update_reservoirs_from_sf(run_globals_t *run_globals, galaxy_t *gal, double new_stars)
{
  if (new_stars > 0)
  {
    double metallicity;
    double current_time;

    // update the galaxy's SFR value
    gal->Sfr += new_stars / gal->dt;
    assert(gal->Sfr >= 0);

    // update the stellar mass history
    gal->NewStars[0] += new_stars;

    // instantaneous recycling approximation of stellar mass
    metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

    gal->ColdGas           -= new_stars;
    gal->MetalsColdGas     -= new_stars * metallicity;
    gal->StellarMass       += new_stars;
    gal->GrossStellarMass  += new_stars;
    gal->MetalsStellarMass += new_stars * metallicity;

    // update the luminosities
    current_time = run_globals->LTTime[gal->LastIdentSnap] - 0.5 * gal->dt;
    add_to_luminosities(run_globals, gal, new_stars, metallicity, current_time);

    // Check the validity of the modified reservoir values.
    // Note that the ColdGas reservers *can* be negative at this point.  This
    // is because some fraction of the stars in this burst will go nova and
    // return mass to the ISM.  This will be accounted for when we update the
    // reservoirs due to supernova feedback.
    if (gal->StellarMass < 0)
      gal->StellarMass = 0.0;
    if (gal->MetalsStellarMass < 0)
      gal->MetalsStellarMass = 0.0;
  }
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
    double m_recycled;
    double new_metals;

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

    if (m_stars > gal->ColdGas)
      m_stars = gal->ColdGas;

    // calculate the total supernova feedback which would occur if this star
    // formation happened continuously and evenly throughout the snapshot
    contemporaneous_supernova_feedback(run_globals, gal, &m_stars, snapshot, &m_reheat, &m_eject, &m_recycled, &new_metals);

    // update the baryonic reservoirs (note that the order we do this in will change the result!)
    update_reservoirs_from_sf(run_globals, gal, m_stars);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, m_recycled, new_metals);
  }
}
