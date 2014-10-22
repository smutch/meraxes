#include "meraxes.h"
#include <math.h>
#include <assert.h>

void update_reservoirs_from_sn_feedback(galaxy_t *gal, double m_reheat, double m_eject, double m_recycled, double new_metals)
{
  double metallicity;
  galaxy_t *central;

  // If this is a ghost then it doesn't have an identified halo at this
  // snapshot.  We will therefore dump all of the reheated gas into the ghost's
  // hot halo, to be recollected and distributed when the ghost is reidentified
  // at a later time.
  if (gal->ghost_flag)
    central = gal;
  else
    central = gal->Halo->FOFGroup->FirstHalo->Galaxy;

  gal->StellarMass -= m_recycled;
  gal->ColdGas     += m_recycled;

  // assuming instantaneous recycling approximation and enrichment from SNII
  // only, work out the mass of metals returned to the ISM by this SF burst
  if (gal->ColdGas > 1e-10)
    gal->MetalsColdGas += new_metals;
  else
    central->MetalsHotGas += new_metals;

  // make sure we aren't trying to use more cold gas than is available...
  if (m_reheat > gal->ColdGas)
    m_reheat = gal->ColdGas;

  metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

  gal->ColdGas          -= m_reheat;
  gal->MetalsColdGas    -= m_reheat * metallicity;
  central->MetalsHotGas += m_reheat * metallicity;
  central->HotGas       += m_reheat;

  // If this is a ghost then we don't know what the ejected mass is as we don't
  // know the properties of the halo!
  if (!gal->ghost_flag)
  {
    metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);

    if (m_eject > central->HotGas)
      m_eject = central->HotGas;

    central->HotGas           -= m_eject;
    central->MetalsHotGas     -= m_eject * metallicity;
    central->EjectedGas       += m_eject;
    central->MetalsEjectedGas += m_eject * metallicity;
  }

  // Check the validity of the modified reservoir values
  if (central->HotGas < 0)
    central->HotGas = 0.0;
  if (central->MetalsHotGas < 0)
    central->MetalsHotGas = 0.0;
  if (gal->ColdGas < 0)
    gal->ColdGas = 0.0;
  if (gal->MetalsColdGas < 0)
    gal->MetalsColdGas = 0.0;
  if (gal->StellarMass < 0)
    gal->StellarMass = 0.0;
  if (central->EjectedGas < 0)
    central->EjectedGas = 0.0;
  if (central->MetalsEjectedGas < 0)
    central->MetalsEjectedGas = 0.0;
}


static inline double calc_ejected_mass(
  double m_reheat,
  double sn_energy,
  double Vvir)
{
  double m_eject = 0.0;

  if (m_reheat > 0)
  {
    double Vvir_sqrd                = Vvir * Vvir;
    double reheated_energy          = 0.5 * m_reheat * Vvir_sqrd;
    double specific_hot_halo_energy = 0.5 * Vvir_sqrd;
    m_eject = (sn_energy - reheated_energy) / specific_hot_halo_energy;
    if (m_eject < 0)
      m_eject = 0.0;
  }

  return m_eject;
}


static inline double calc_eta_sn(run_globals_t *run_globals, double m_high, double m_low, double *snII_frac)
{
  // work out the number of supernova per 1e10 Msol formed at the current time
  double exponent  = run_globals->params.physics.IMFSlope + 1.0; // should be -1.35 for Salpeter
  double const_phi = run_globals->params.physics.IMFNormConst;   // should be 0.1706 for Salpeter
  double eta_SNII  = 7.4319792e-3; // total number of type II SN per solar mass of burst

  double eta_sn = const_phi * 1.0 / exponent * (pow(m_high, exponent) - pow(m_low, exponent));

  *snII_frac = eta_sn / eta_SNII;
  eta_sn    *= 1.e10;

  if (*snII_frac > 1)
    *snII_frac = 1.0;

  assert((eta_sn >= 0) && (*snII_frac >= 0));
  return eta_sn;
}

static inline double calc_sn_energy(run_globals_t *run_globals, double stars, double Vmax, double eta_sn)
{
  // work out the energy produced by the supernova and add it to our total at this snapshot
  double E_sn          = run_globals->params.physics.EnergyPerSN / run_globals->units.UnitEnergy_in_cgs;
  double SnEjectionScaling = run_globals->params.physics.SnEjectionScaling;
  double SnEjectionNorm = run_globals->params.physics.SnEjectionNorm;
  double SnEjectionEff = run_globals->params.physics.SnEjectionEff;
  double sn_energy;

  if (SnEjectionScaling != 0)
  {
    SnEjectionEff *= 0.5 + pow(Vmax/SnEjectionNorm, -SnEjectionScaling);
    if (SnEjectionEff > 1.0)
      SnEjectionEff = 1.0;
  }

  sn_energy  = 0.5 * SnEjectionEff * stars * E_sn * eta_sn;
  assert(sn_energy >= 0);

  return sn_energy;
}


double calc_recycled_frac(run_globals_t *run_globals, double m_high, double m_low, double *burst_mass_frac)
{
  // calculate the mass ejected (from fraction of total SN-II that have gone off) from this burst
  double const_phi = run_globals->params.physics.IMFNormConst;    // should be 0.1706 for Salpeter
  double exponent  = run_globals->params.physics.IMFSlope + 2.0;

  double burst_recycled_frac = const_phi * 1.0 / exponent * (pow(m_high, exponent) - pow(m_low, exponent));
  double frac_mass_SSP_above_SNII = 0.14417;  // Fraction of SSP with M>8Msol

  if (burst_recycled_frac < 0)
    SID_log("WTF? %g (%g, %g)", SID_LOG_COMMENT, burst_recycled_frac, m_low, m_high);

  assert(burst_recycled_frac >= 0);

  *burst_mass_frac = burst_recycled_frac / frac_mass_SSP_above_SNII;
  if (*burst_mass_frac > 1.0)
    *burst_mass_frac = 1.0;

  assert(*burst_mass_frac >= 0);

  return burst_recycled_frac;
}


double sn_m_low(double log_dt)
{
  // log_dt must be in units of log10(dt/Myr)
  // returned value is in units of Msol

  // This is a fit to the H+He core burning lifetimes of stars of varying
  // masses from Table 14 of Portinari, L., Chiosi, C. & Bressan, A.
  // Galactic chemical enrichment with new metallicity dependent stellar
  // yields.  Astronomy and Astrophysics 334, 505–539 (1998).

  double const_a = 0.80680268;
  double const_b = -2.81547136;
  double const_c = -5.73419164;
  double const_d = 0.475119568;
  double m_high  = 120.0;           // highest mass star produced in stellar mass burst (Msol)
  double m_low;

  m_low = pow(10.0, (const_a / log_dt) + const_b * exp(const_c / log_dt) + const_d);

  // the fitting formula for m_low is only valid until t=const_d
  if (m_low > m_high)
    m_low = m_high;
  else if (m_low < 0.0)
    m_low = 0.0;

  return m_low;
}


void delayed_supernova_feedback(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{
  run_units_t *units = &(run_globals->units);
  double SnReheatScaling = run_globals->params.physics.SnReheatScaling;
  double SnReheatNorm = run_globals->params.physics.SnReheatNorm;
  double SnReheatEff = run_globals->params.physics.SnReheatEff;
  double *LTTime     = run_globals->LTTime;

  double m_high;
  double m_low;
  double eta_sn;
  double burst_recycled_frac;
  double burst_mass_frac;
  double snII_frac;
  double log_dt;
  double sn_energy  = 0.0;
  double m_reheat   = 0.0;
  double m_eject    = 0.0;
  double m_recycled = 0.0;
  double new_metals = 0.0;
  double m_stars;

  // If we are at snapshot < N_HISTORY_SNAPS-1 then only try to look back to snapshot 0
  int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;

  // scale the reheating efficiency
  if (SnReheatScaling != 0)
    SnReheatEff *= 0.5 + pow(gal->Vmax/SnReheatNorm, -SnReheatScaling);

  // Loop through each of the last `N_HISTORY_SNAPS` recorded stellar mass
  // bursts and calculate the amount of energy and mass that they will release
  // in the current time step.
  for (int i_burst = 1; i_burst < n_bursts; i_burst++)
  {
    m_stars = gal->NewStars[i_burst];

    // Only need to do this if any stars formed in this history bin
    if (m_stars > 1e-10)
    {
      // work out the higest mass star which would have expended it's H & He
      // core fuel in this time (i.e. the highest mass star that would not
      // already have left the main sequence since this burst occured).
      log_dt = log10(((LTTime[snapshot - i_burst - 1] + LTTime[snapshot - i_burst]) / 2.0 - LTTime[snapshot - 1]) * units->UnitTime_in_Megayears / run_globals->params.Hubble_h);
      m_high = sn_m_low(log_dt);

      // work out the lowest mass star which would have expended it's H & He core
      // fuel in this time
      log_dt = log10(((LTTime[snapshot - i_burst - 1] + LTTime[snapshot - i_burst]) / 2.0 - LTTime[snapshot]) * units->UnitTime_in_Megayears / run_globals->params.Hubble_h);
      m_low  = sn_m_low(log_dt); // Msol

      // calculate the mass recycled from this burst
      burst_recycled_frac = calc_recycled_frac(run_globals, m_high, m_low, &burst_mass_frac);
      m_recycled         += m_stars * burst_recycled_frac;

      // If m_high is > 8.0 Msol then we have already used all of the SN-II in
      // the previous recorded NewStars bins.  We therefore don't need to
      // calculate any supernova feedback quantites.
      if (m_high <= 8.0)
        continue;

      if (m_low < 8.0)
        m_low = 8.0;

      // work out the number of supernova per unit stellar mass formed at the
      // current time
      eta_sn = calc_eta_sn(run_globals, m_high, m_low, &snII_frac);

      // increment the total reheated and new metals masses
      m_reheat   += SnReheatEff * snII_frac * m_stars;
      new_metals += run_globals->params.physics.Yield * m_stars * burst_mass_frac;

      // now work out the energy produced by the supernova and add it to our total at this snapshot
      sn_energy += calc_sn_energy(run_globals, m_stars, gal->Vmax, eta_sn);
    }
  }

  assert(m_reheat >= 0);
  assert(m_recycled >= 0);
  assert(new_metals >= 0);

  if (!gal->ghost_flag)
  {
    // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
    // Note that we use the Vvir of the host group here, as we are assuming that
    // only the group holds a hot halo (which is stored by the central galaxy).
    m_eject = calc_ejected_mass(m_reheat, sn_energy, gal->Halo->FOFGroup->Vvir);
  }

  // update the baryonic reservoirs
  update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, m_recycled, new_metals);
}


static void backfill_ghost_NewStars(run_globals_t *run_globals, galaxy_t *gal, double m_stars, int snapshot)
{
  if ((snapshot - gal->LastIdentSnap) <= N_HISTORY_SNAPS)
  {
    double *LTTime    = run_globals->LTTime;
    double burst_time = LTTime[gal->LastIdentSnap] - gal->dt * 0.5;

    for (int ii = 1; ii < N_HISTORY_SNAPS; ii++)
      if (LTTime[snapshot - ii] > burst_time)
      {
        gal->NewStars[ii - 1] += m_stars;
        break;
      }
  }
}


void contemporaneous_supernova_feedback(
  run_globals_t *run_globals,
  galaxy_t      *gal,
  double        *m_stars,
  int            snapshot,
  double        *m_reheat,
  double        *m_eject,
  double        *m_recycled,
  double        *new_metals)
{
  // Here we approximate a constant SFR accross the timestep by a single burst
  // at t=0.5*dt.  This is a pretty good approximation (to within ~15% of the
  // true number of SN that would have gone of by the end of the timestep for a
  // constant SFR).

  run_units_t *units = &(run_globals->units);
  double SnReheatScaling = run_globals->params.physics.SnReheatScaling;
  double SnReheatNorm = run_globals->params.physics.SnReheatNorm;
  double SnReheatEff = run_globals->params.physics.SnReheatEff;
  bool Flag_IRA  = (bool)(run_globals->params.physics.Flag_IRA);

  double m_high = 120.0;  // Msol
  double m_low  = 8.0;
  double eta_sn;
  double burst_recycled_frac = run_globals->params.physics.SfRecycleFraction;
  double burst_mass_frac = 1.0;
  double snII_frac;
  double log_dt;
  double sn_energy = 0.0;

  // init (just in case!)
  *m_reheat = *m_recycled = *new_metals = *m_eject = 0.0;

  // scale the reheating efficiency
  if (SnReheatScaling != 0)
    SnReheatEff *= 0.5 + pow(gal->Vmax/SnReheatNorm, -SnReheatScaling);

  // work out the lowest mass star which would have expended it's H & He core
  // fuel in this time
  // N.B. If Flag_IRA is true then m_low and burst_recycled_frac will equal values in above declaration
  if (!Flag_IRA)
  {
    assert(snapshot > 0);
    log_dt = log10(gal->dt * 0.5 * units->UnitTime_in_Megayears / run_globals->params.Hubble_h);
    m_low  = sn_m_low(log_dt); // Msol

    // calculate the mass reheated (from fraction of total SN-II that have gone off) from this burst
    burst_recycled_frac = calc_recycled_frac(run_globals, m_high, m_low, &burst_mass_frac);
  }

  *m_recycled = *m_stars * burst_recycled_frac;

  if (m_low < 8.0)
    m_low = 8.0;

  // work out the number of supernova per unit stellar mass formed at the current time
  eta_sn = calc_eta_sn(run_globals, m_high, m_low, &snII_frac);

  // calculate the total reheated
  *m_reheat = SnReheatEff * snII_frac * *m_stars;

  // attenuate the star formation if necessary, so that we are being consistent
  // if (*m_reheat + *m_stars - *m_recycled > gal->ColdGas)
  if (*m_reheat + *m_stars > gal->ColdGas)
  {
    // double frac = gal->ColdGas / (*m_reheat + *m_stars - *m_recycled);
    double frac = gal->ColdGas / (*m_reheat + *m_stars);
    *m_reheat *= frac;
    *m_stars  *= frac;
    // *m_recycled *= frac;
  }

  // how much new metals will be created by this burst?
  *new_metals = run_globals->params.physics.Yield * *m_stars * burst_mass_frac;

  // now work out the energy produced by the supernova
  sn_energy = calc_sn_energy(run_globals, *m_stars, gal->Vmax, eta_sn);

  assert(*m_reheat >= 0);
  assert(*m_recycled >= 0);
  assert(*new_metals >= 0);

  // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
  // Note that we use the Vvir of the host group here, as we are assuming that
  // only the group holds a hot halo (which is stored by the central galaxy).
  *m_eject = calc_ejected_mass(*m_reheat, sn_energy, gal->Halo->FOFGroup->Vvir);

  // If this is a reidentified ghost, then back fill NewStars to reflect this
  // new SF burst.
  if (!Flag_IRA && (gal->LastIdentSnap < (snapshot - 1)))
    backfill_ghost_NewStars(run_globals, gal, *m_stars, snapshot);
}
