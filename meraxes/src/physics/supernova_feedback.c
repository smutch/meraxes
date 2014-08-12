#include "meraxes.h"
#include <math.h>
#include <assert.h>

void update_reservoirs_from_sn_feedback(galaxy_t *gal, double m_reheat, double m_eject, double new_metals)
{
  double metallicity;
  galaxy_t *central = gal->Halo->FOFGroup->FirstHalo->Galaxy;

  // assuming instantaneous recycling approximation and enrichment from SNII
  // only, work out the mass of metals returned to the ISM by this SF burst
  if (gal->ColdGas > 1e-10)
    gal->MetalsColdGas += new_metals;
  else
    central->MetalsHotGas += new_metals;

  metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

  central->HotGas       += m_reheat;
  central->MetalsHotGas += m_reheat * metallicity;
  gal->ColdGas          -= m_reheat;
  gal->MetalsColdGas    -= m_reheat * metallicity;

  metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);

  if (m_eject > central->HotGas)
    m_eject = central->HotGas;

  central->EjectedGas       += m_eject;
  central->MetalsEjectedGas += m_eject * metallicity;
  central->HotGas           -= m_eject;
  central->MetalsHotGas     -= m_eject * metallicity;
}


static inline double calc_ejected_mass(
    double m_reheat,
    double sn_energy,
    double Vvir)
{
  double m_eject = 0.0;

  if(m_reheat > 0)
  {
    double Vvir_sqrd = Vvir*Vvir;
    double reheated_energy = 0.5 * m_reheat * Vvir_sqrd;
    double specific_hot_halo_energy = 0.5 * Vvir_sqrd;
    m_eject = (sn_energy - reheated_energy)/specific_hot_halo_energy;
    if (m_eject < 0)
      m_eject = 0.0;
  }

  return m_eject;
}


static inline double calc_eta_sn(run_globals_t *run_globals, double m_high, double m_low)
{
  // work out the number of supernova per 1e10 Msol formed at the current time
  double exponent = run_globals->params.physics.IMFSlope + 1.0;  // should be -1.35 for Salpeter
  double const_phi = run_globals->params.physics.IMFNormConst;   // should be 0.1706 for Salpeter
  double eta_sn = const_phi * 1.0/exponent * 1e10 * (pow(m_high, exponent) - pow(m_low, exponent));
  assert(eta_sn >= 0);
  return eta_sn;
}

static inline double calc_sn_energy(run_globals_t *run_globals, double stars, double eta_sn)
{
  // work out the energy produced by the supernova and add it to our total at this snapshot
  double E_sn = run_globals->params.physics.EnergyPerSN / run_globals->units.UnitEnergy_in_cgs;
  double SnEjectionEff = run_globals->params.physics.SnEjectionEff;
  double sn_energy = 0.5 * SnEjectionEff * stars * E_sn * eta_sn;
  assert(sn_energy >= 0);
  return sn_energy;
}


static inline void calc_sn_frac(run_globals_t *run_globals, double m_high, double m_low, double *sf_frac, double *sn_frac)
{
  // calculate the mass ejected (from fraction of total SN-II that have gone off) from this burst
  double  m_frac_SNII   = 0.1442; // fraction of total stellar mass in stars more massive than 8Msol
  double  const_phi = run_globals->params.physics.IMFNormConst;   // should be 0.1706 for Salpeter
  double exponent = run_globals->params.physics.IMFSlope + 2.0;
  *sf_frac = const_phi * 1.0/exponent * (pow(m_high, exponent) - pow(m_low, exponent));
  *sn_frac = *sf_frac / m_frac_SNII;
  assert((*sn_frac <= 1) && (*sn_frac >= 0));
}


static double inline sn_m_low(double log_dt)
{
  // log_dt must be in units of log10(dt/Myr)
  // returned value is in units of Msol

  // This is a fit to the H+He core burning lifetimes of stars of varying
  // masses from Table 14 of Portinari, L., Chiosi, C. & Bressan, A.
  // Galactic chemical enrichment with new metallicity dependent stellar
  // yields.  Astronomy and Astrophysics 334, 505–539 (1998).

  double const_a     = 0.80680268;
  double const_b     = -2.81547136;
  double const_c     = -5.73419164;
  double const_d     = 0.475119568;
  double m_high      = 120.0;       // highest mass star produced in stellar mass burst (Msol)
  double m_low;

  m_low = pow(10.0, (const_a / log_dt) + const_b * exp(const_c / log_dt) + const_d);

  // the fitting formula for m_low is only valid until t=const_d
  if (m_low > m_high)
    m_low = m_high;
  // we are only including SN-II which corresponds to m_low >= 8Msol
  else if (m_low < 8.0)
    m_low = 8.0;

  return m_low;
}


void supernova_feedback(run_globals_t *run_globals, galaxy_t *gal, double *m_stars, double *m_reheat, double *m_eject, double *m_recycled, double *new_metals, int snapshot)
{
  run_units_t *units   = &(run_globals->units);
  double SnReheatEff   = run_globals->params.physics.SnReheatEff;
  double *LTTime       = run_globals->LTTime;
  fof_group_t *fof_group = gal->Halo->FOFGroup;

  double m_high = 120.0;  // Msol
  double m_low;
  double eta_sn;
  double sf_frac;
  double sn_frac;
  double sn_energy;
  double log_dt;

  // work out the lowest mass star which would have expended it's H & He core
  // fuel in this time
  // log_dt = log10((LTTime[snapshot-1] - LTTime[snapshot]) * units->UnitTime_in_Megayears / run_globals->params.Hubble_h);
  // m_low = sn_m_low(log_dt);  // Msol
  m_low = 8.0;

  // work out the number of supernova per unit stellar mass formed at the current time
  eta_sn = calc_eta_sn(run_globals, m_high, m_low);

  // now work out the energy produced by the supernova and add it to our total at this snapshot
  sn_energy = calc_sn_energy(run_globals, (*m_stars), eta_sn);

  // finally calculate the mass reheated (from fraction of total SN-II that have gone off) from this burst
  calc_sn_frac(run_globals, m_high, m_low, &sf_frac, &sn_frac);
  *m_reheat  = SnReheatEff * sn_frac * (*m_stars);

  // make sure we aren't trying to use more cold gas than is available...
  if (((*m_stars) + (*m_reheat)) > gal->ColdGas)
  {
    double factor = gal->ColdGas / ((*m_stars) + (*m_reheat));
    *m_stars  *= factor;
    *m_reheat *= factor;
  }

  // how much mass is recycled and ejected due to this star formation episode? (ala Croton+ 2006)
  // TODO: remove sf_frac = 0.43 for IRA
  *m_recycled = (*m_stars) * run_globals->params.physics.SfRecycleFraction;// * sf_frac;
  *m_eject = calc_ejected_mass((*m_reheat), sn_energy, fof_group->Vvir);
  *new_metals = run_globals->params.physics.Yield * (*m_stars);

  // make sure we are being consistent
  assert(*m_reheat >= 0);
  assert(*m_recycled >= 0);
  assert(*new_metals >= 0);
  assert(*m_eject >= 0);

}
