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

static void update_reservoirs_from_sn_feedback(galaxy_t *gal, double m_reheat, double m_eject)
{
  double metallicity;
  double m_to_hot;

  metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);
  m_to_hot = m_reheat - m_eject;
  if(m_to_hot < 0)
    m_to_hot = 0.0;

  gal->HotGas           += m_to_hot;
  gal->MetalsHotGas     += m_to_hot * metallicity;
  gal->EjectedGas       += m_eject;
  gal->MetalsEjectedGas += m_eject * metallicity;
  gal->ColdGas          -= m_reheat;
  gal->MetalsColdGas    -= m_reheat * metallicity;

}

void insitu_star_formation_and_feedback(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{

  // These star formation & supernova feedback prescriptions are taken directly
  // from Croton+ 2006

  // there is no point doing anything if there is no cold gas!
  if(gal->ColdGas > 1e-8)
  {
    double r_disk;
    double m_crit;
    double m_stars;
    double m_reheat;
    double m_eject;
    double factor;

    double SnReheatEff = run_globals->params.physics.SnReheatEff;
    double SnEjectionEff = run_globals->params.physics.SnEjectionEff;
    double SfEfficiency = run_globals->params.physics.SfEfficiency;
    double SfRecycleFraction = run_globals->params.physics.SfRecycleFraction;
    double sqrt_2 = 1.414213562;
    double sn_velocity = 630.0;  // km/s

    // calculate disk scalelength using Mo, Mau & White (1998) eqn. 12 and
    // multiply it by 3 to approximate the star forming region size (ala
    // Croton+ 2006).
    r_disk = gal->Spin * gal->Rvir / sqrt_2 * 3.0;

    // what is the critical mass within r_crit?
    m_crit = 0.19 * gal->Vvir * r_disk;

    if(gal->ColdGas > m_crit)
      gal->Sfr = SfEfficiency * (gal->ColdGas - m_crit) / r_disk * gal->Vmax;
    else
    {
      gal->Sfr = 0.0;
      return;
    }

    // given the timestep, what mass of stars does this correspond to?
    m_stars = gal->Sfr * gal->dt;

    // following the SN feedback model Croton+ 2006, what mass of cold gas will
    // we end up reheating due to this star formation episode?
    m_reheat = SnReheatEff * m_stars;

    // make sure we aren't trying to use more cold gas than is available...
    if ((m_stars + m_reheat) > gal->ColdGas)
    {
      factor = gal->ColdGas / (m_stars + m_reheat);
      m_stars *= factor;
      gal->Sfr *= factor;
      m_reheat *= factor;
    }

    // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
    m_eject = (SnEjectionEff * pow(sn_velocity/gal->Vvir, 2) - SnReheatEff) * m_stars;

    // make sure we are being consistent
    if(m_eject > m_reheat)
      m_eject = m_reheat;
    else if(m_eject < 0)
      m_eject = 0.0;

    // instantaneous recycling approximation
    m_stars = (1 - SfRecycleFraction) * m_stars;

    // update the baryonic reserviors (note the order makes a difference here!)
    update_reservoirs_from_sf(gal, m_stars);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject);

  }
}
