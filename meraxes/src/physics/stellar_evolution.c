#include "meraxes.h"
#include <math.h>
#include <assert.h>


void evolve_stellar_pops(galaxy_t *gal, int snapshot)
{
  if (gal->mwmsa_num > 0)
  {
    double log_dt;
    double m_stars = gal->mwmsa_denom;
    double mwmsa   = gal->mwmsa_num / gal->mwmsa_denom;
    double m_high, m_low, burst_recycled_frac, m_recycled;
    run_units_t *units = &(run_globals.units);
    double *LTTime     = run_globals.LTTime;
    double burst_mass_frac;

    // Use the pre-calculated mwmsa numerator and denomenator to get the mwmsa
    // of all stars formed before what is tracked in the NewStars property
    log_dt = log10((mwmsa - LTTime[snapshot - 1]) * units->UnitTime_in_Megayears / run_globals.params.Hubble_h);

    // work out the higest mass star which would have expended it's H & He
    // core fuel in this time (i.e. the highest mass star that would not
    // already have left the main sequence since this burst occured).
    m_high = sn_m_low(log_dt);

    // work out the lowest mass star which would have expended it's H & He core
    // fuel in this time
    log_dt = log10((mwmsa - LTTime[snapshot]) * units->UnitTime_in_Megayears / run_globals.params.Hubble_h);
    m_low  = sn_m_low(log_dt); // Msol

    burst_recycled_frac = calc_recycled_frac(m_high, m_low, &burst_mass_frac);
    m_recycled          = m_stars * burst_recycled_frac;
    update_reservoirs_from_sn_feedback(gal, 0.0, 0.0, m_recycled, 0.0);
  }
}
