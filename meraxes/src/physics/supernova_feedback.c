#include "meraxes.h"
#include <math.h>

static void update_reservoirs_from_sn_feedback(galaxy_t *gal, double m_reheat, double m_eject)
{
  double metallicity;
  galaxy_t *central = gal->Halo->FOFGroup->FirstHalo->Galaxy;

  metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

  central->HotGas           += m_reheat;
  central->MetalsHotGas     += m_reheat * metallicity;
  gal->ColdGas              -= m_reheat;
  gal->MetalsColdGas        -= m_reheat * metallicity;

  metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);

  if(m_eject > central->HotGas)
    m_eject = central->HotGas;

  central->EjectedGas       += m_eject;
  central->MetalsEjectedGas += m_eject * metallicity;
  central->HotGas           -= m_eject;
  central->MetalsHotGas     -= m_eject * metallicity;

  if(central->id_MBP == DEBUG_MBP)
    fprintf(stderr, "gal->ID %d just ejected %.3e worth of gas into id_MBP = %d\n", gal->ID, m_eject, DEBUG_MBP);

}


void supernova_feedback(run_globals_t *run_globals, galaxy_t *gal, double m_stars, double merger_mass_ratio, int snapshot)
{

    // NOTE: m_stars should be mass of stars formed **before** instantaneous
    // recycling approximation is applied

    double m_reheat;
    double m_eject;
    double factor;
    galaxy_t *central = gal->Halo->FOFGroup->FirstHalo->Galaxy;

    double SnReheatEff = run_globals->params.physics.SnReheatEff;
    double SnEjectionEff = run_globals->params.physics.SnEjectionEff;
    double sn_velocity = 630.0;  // km/s

    // following the SN feedback model Croton+ 2006, what mass of cold gas will
    // we end up reheating due to this star formation episode?
    m_reheat = SnReheatEff * m_stars;

    // make sure we aren't trying to use more cold gas than is available...
    if ((m_stars + m_reheat) > gal->ColdGas)
    {
      factor = gal->ColdGas / (m_stars + m_reheat);
      m_stars *= factor;
      m_reheat *= factor;
    }

    // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
    m_eject = (SnEjectionEff * pow(sn_velocity/central->Vvir, 2) - SnReheatEff) * m_stars;

    // make sure we are being consistent
    if(m_eject < 0)
      m_eject = 0.0;

    if(gal->id_MBP == DEBUG_MBP)
    {
      fprintf(stderr, "DEBUG SN FEEDBACK (%d)\n", snapshot);
      fprintf(stderr, "m_stars = %.3e\n", m_stars);
      fprintf(stderr, "m_reheat = %.3e\n", m_reheat);
      fprintf(stderr, "m_eject = %.3e\n", m_eject);
      fprintf(stderr, "m_cold = %.3e\n", gal->ColdGas);
    }

    // update the baryonic reservoirs (note the order makes a difference here!)
    update_reservoirs_from_sf(run_globals, gal, m_stars, merger_mass_ratio, snapshot);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject);

}
