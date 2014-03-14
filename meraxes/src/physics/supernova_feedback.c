#include "meraxes.h"
#include <math.h>

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


void supernova_feedback(run_globals_t *run_globals, galaxy_t *gal, double m_stars)
{

    // NOTE: m_stars should be mass of stars formed **before** instantaneous
    // recycling approximation is applied

    double m_reheat;
    double m_eject;
    double factor;

    double SnReheatEff = run_globals->params.physics.SnReheatEff;
    double SnEjectionEff = run_globals->params.physics.SnEjectionEff;
    double sn_velocity = 630.0;  // km/s

    // following the SN feedback model Croton+ 2006, what mass of cold gas will
    // we end up reheating due to this star formation episode?
    m_reheat = SnReheatEff * m_stars;

    // make sure we aren't trying to use more cold gas than is available...
    if (m_stars + m_reheat > gal->ColdGas)
    {
      factor = gal->ColdGas / (m_stars+ m_reheat);
      m_stars *= factor;
      m_reheat *= factor;
    }

    // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
    m_eject = (SnEjectionEff * pow(sn_velocity/gal->Vvir, 2) - SnReheatEff) * m_stars;

    // make sure we are being consistent
    if(m_eject > m_reheat)
      m_eject = m_reheat;
    else if(m_eject < 0)
      m_eject = 0.0;

    // update the baryonic reservoirs (note the order makes a difference here!)
    update_reservoirs_from_sf(run_globals, gal, m_stars);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject);

}
