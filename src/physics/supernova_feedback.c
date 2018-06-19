#include "meraxes.h"
#include <assert.h>
#include <math.h>

void update_reservoirs_from_sn_feedback(galaxy_t* gal, double m_reheat, double m_eject, 
                                        double m_recycled, double new_metals)
{
    double metallicity;
    galaxy_t* central;

    // If this is a ghost then it doesn't have an identified halo at this
    // snapshot.  We will therefore dump all of the reheated gas into the ghost's
    // hot halo, to be recollected and distributed when the ghost is reidentified
    // at a later time.
    if (gal->ghost_flag)
        central = gal;
    else
        central = gal->Halo->FOFGroup->FirstOccupiedHalo->Galaxy;

    gal->StellarMass -= m_recycled;
    //gal->MetalsStellarMass -= new_metals;
    gal->ColdGas += m_recycled;

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

    gal->ColdGas -= m_reheat;
    gal->MetalsColdGas -= m_reheat * metallicity;
    central->MetalsHotGas += m_reheat * metallicity;
    central->HotGas += m_reheat;

    // If this is a ghost then we don't know what the real ejected mass is as we
    // don't know the properties of the halo!
    if (!gal->ghost_flag) {
        metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);

        if (m_eject > central->HotGas)
            m_eject = central->HotGas;

        central->HotGas -= m_eject;
        central->MetalsHotGas -= m_eject * metallicity;
        central->EjectedGas += m_eject;
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
    if (gal->MetalsStellarMass < 0)
        gal->MetalsStellarMass = 0.0;
    if (central->EjectedGas < 0)
        central->EjectedGas = 0.0;
    if (central->MetalsEjectedGas < 0)
        central->MetalsEjectedGas = 0.0;
}

static inline double calc_ejected_mass(
    double* m_reheat,
    double sn_energy,
    double Vvir,
    double fof_Vvir)
{
    double m_eject = 0.0;

    if (*m_reheat > 0) {
        if (run_globals.params.physics.Flag_ReheatToFOFGroupTemp)
            Vvir = fof_Vvir;

        // Begin by calculating if we have enough energy to get m_reheat of gas to
        // Tvir of the host *subhalo*.
        double Vvir_sqrd = Vvir * Vvir;
        double reheated_energy = 0.5 * (*m_reheat) * Vvir_sqrd;
        double specific_hot_halo_energy = 0.5 * Vvir_sqrd;

        m_eject = (sn_energy - reheated_energy) / specific_hot_halo_energy;

        if (m_eject <= 0) {
            // If there is not enough energy to reheat all of the gas to Tvir of the
            // subhalo then how much can we reheat?
            m_eject = 0.0;
            *m_reheat = 2.0 * sn_energy / Vvir_sqrd;
        }
        else if (fof_Vvir > 0) {
            // If we were able to reheat all of the mass with energy left to spare,
            // is there enough energy to further eject gas from the host *FOF group*?
            Vvir_sqrd = fof_Vvir * fof_Vvir;
            reheated_energy = 0.5 * (*m_reheat) * Vvir_sqrd;
            specific_hot_halo_energy = 0.5 * Vvir_sqrd;

            m_eject = (sn_energy - reheated_energy) / specific_hot_halo_energy;

            if (m_eject < 0)
                m_eject = 0.0;
        }
    }

    return m_eject;
}

static inline double calc_eta_sn(double m_high, double m_low, double* snII_frac)
{
    // work out the number of supernova per 1e10 Msol formed at the current time
    double exponent = run_globals.params.physics.IMFSlope + 1.0; // should be -1.35 for Salpeter
    double const_phi = run_globals.params.physics.IMFNormConst; // should be 0.1706 for Salpeter
    double eta_SNII = run_globals.params.physics.eta_SNII; // total number of type II SN per solar mass of burst

    double eta_sn = const_phi * 1.0 / exponent * (pow(m_high, exponent) - pow(m_low, exponent));

    *snII_frac = eta_sn / eta_SNII;
    eta_sn *= 1.e10;

    if (*snII_frac > 1)
        *snII_frac = 1.0;

    assert((eta_sn >= 0) && (*snII_frac >= 0));
    return eta_sn;
}


static inline double calc_sn_reheat_eff(galaxy_t *gal, int snapshot)
{
    double Vmax = gal->Vmax;    // Vmax is in a unit of km/s
    double zplus1 = 1. + run_globals.ZZ[snapshot];
    physics_params_t *params = &run_globals.params.physics;
    int SnModel = params->SnModel;
    double SnReheatRedshiftDep = params->SnReheatRedshiftDep;
    double SnReheatEff = params->SnReheatEff;
    double SnReheatLimit = params->SnReheatLimit;
    if (SnModel == 1) {
        double SnReheatScaling = params->SnReheatScaling;
        double SnReheatNorm = params->SnReheatNorm;
        SnReheatEff *= pow(zplus1/4., SnReheatRedshiftDep) \
                       *(.5 + pow(Vmax/SnReheatNorm, -SnReheatScaling));
    }
    else {
        if (Vmax < 60.)
            SnReheatEff *= pow(zplus1/4., SnReheatRedshiftDep)*pow(Vmax/60., -3.2);
        else
            SnReheatEff *= pow(zplus1/4., SnReheatRedshiftDep)*pow(Vmax/60., -1);
    }
    if (SnReheatEff < SnReheatLimit)
        return SnReheatEff;
    else
        return SnReheatLimit;
}


static inline double calc_sn_ejection_eff(galaxy_t *gal, int snapshot)
{
    double Vmax = gal->Vmax;    // Vmax is in a unit of km/s
    double zplus1 = 1. + run_globals.ZZ[snapshot];
    physics_params_t *params = &run_globals.params.physics;
    int SnModel = params->SnModel;
    double SnEjectionRedshiftDep = params->SnEjectionRedshiftDep;
    double SnEjectionEff = params->SnEjectionEff;
    if (SnModel == 1) {
        double SnEjectionNorm = params->SnEjectionNorm;
        double SnEjectionScaling = params->SnEjectionScaling;
        SnEjectionEff *= pow(zplus1/4., SnEjectionRedshiftDep)\
                         *(.5 + pow(Vmax/SnEjectionNorm, -SnEjectionScaling));
    }
    else {
        if (Vmax < 60.)
            SnEjectionEff *= pow(zplus1/4., SnEjectionRedshiftDep)*pow(Vmax/60., -3.2);
        else
            SnEjectionEff *= pow(zplus1/4., SnEjectionRedshiftDep)*pow(Vmax/60., -1);
    }
    if (SnEjectionEff < 1.)
        return SnEjectionEff;
    else
        return 1.;
}


static inline double calc_sn_energy(double stars, double Vmax, double eta_sn, int snapshot)
{
    // work out the energy produced by the supernova and add it to our total at this snapshot
    double zplus1 = 1. + run_globals.ZZ[snapshot];
    double E_sn = run_globals.params.physics.EnergyPerSN / run_globals.units.UnitEnergy_in_cgs;
    double SnEjectionScaling = run_globals.params.physics.SnEjectionScaling;
    double SnEjectionNorm = run_globals.params.physics.SnEjectionNorm;
    double SnEjectionEff = run_globals.params.physics.SnEjectionEff;
    double sn_energy;

    if (SnEjectionScaling != 0) {
        //SnEjectionEff *= 0.5 + pow(Vmax / SnEjectionNorm, -SnEjectionScaling);
        if (Vmax < 60.)
            SnEjectionEff *= pow(zplus1, 1.3)*pow(Vmax/60., -3.2);
        else
            SnEjectionEff *= pow(zplus1, 1.3)*pow(Vmax/60., -1);
        if (SnEjectionEff > 1.0)
            SnEjectionEff = 1.0;
    }

    sn_energy = SnEjectionEff * stars * E_sn * eta_sn;
    assert(sn_energy >= 0);

    return sn_energy;
}

double calc_recycled_frac(double m_high, double m_low, double* burst_mass_frac)
{
    // calculate the mass ejected (from fraction of total SN-II that have gone off) from this burst
    double const_phi = run_globals.params.physics.IMFNormConst; // should be 0.1706 for Salpeter
    double exponent = run_globals.params.physics.IMFSlope + 2.0;

    double burst_recycled_frac = const_phi * 1.0 / exponent * (pow(m_high, exponent) - pow(m_low, exponent));
    double frac_mass_SSP_above_SNII = run_globals.params.physics.frac_mass_SSP_above_SNII; // Fraction of SSP with M>8Msol

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

    // This is a fit to the H+He core burning lifetimes of Z=0.004 stars of varying
    // masses from Table 14 of Portinari, L., Chiosi, C. & Bressan, A.
    // Galactic chemical enrichment with new metallicity dependent stellar
    // yields.  Astronomy and Astrophysics 334, 505--539 (1998).

    double const_a = 0.74729454;
    double const_b = -2.69790558;
    double const_c = -4.76591765;
    double const_d = 0.59339486;
    double m_high = 120.0; // highest mass star produced in stellar mass burst (Msol)
    double m_low;

    m_low = pow(10.0, (const_a / log_dt) + const_b * exp(const_c / log_dt) + const_d);

    // the fitting formula for m_low is only valid until t=const_d
    if (m_low > m_high)
        m_low = m_high;
    else if (m_low < 0.0)
        m_low = 0.0;

    return m_low;
}

void delayed_supernova_feedback(galaxy_t* gal, int snapshot)
{
    double sn_energy = 0.0;
    double m_reheat = 0.0;
    double m_eject = 0.0;
    double m_recycled = 0.0;
    double new_metals = 0.0;
    double fof_Vvir;
    // If we are at snapshot < N_HISTORY_SNAPS-1 then only try to look back to snapshot 0
    int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;

    // Loop through each of the last `N_HISTORY_SNAPS` recorded stellar mass
    // bursts and calculate the amount of energy and mass that they will release
    // in the current time step.
    for (int i_burst = 1; i_burst < n_bursts; i_burst++) {
        double m_stars = gal->NewStars[i_burst];

        // Only need to do this if any stars formed in this history bin
        if (m_stars > 1e-10) {
            double metallicity = calc_metallicity(m_stars, gal->NewMetals[i_burst]);
            // Calculate recycled mass and metals by yield tables
            m_recycled += m_stars * get_yield(i_burst, metallicity, Y_TOTAL);
            new_metals += m_stars * get_yield(i_burst, metallicity, Y_TOTAL_METAL);
            // Calculate SNII energy
            sn_energy += get_energy(i_burst, metallicity)*m_stars;
        }
    }

    m_reheat = calc_sn_reheat_eff(gal, snapshot)*sn_energy/get_total_energy();
    sn_energy *= calc_sn_ejection_eff(gal, snapshot);
    // We can only reheat as much gas as we have available.  Let's inforce this
    // now, to ensure that the maximal amount of available energy is used to
    // eject gas from the system.
    if (m_reheat > gal->ColdGas)
        m_reheat = gal->ColdGas;

    assert(m_reheat >= 0);
    assert(m_recycled >= 0);
    assert(new_metals >= 0);

    // how much mass is ejected due to this star formation episode?
    if (!gal->ghost_flag)
        fof_Vvir = gal->Halo->FOFGroup->Vvir;
    else
        fof_Vvir = -1;

    m_eject = calc_ejected_mass(&m_reheat, sn_energy, gal->Vvir, fof_Vvir);

    // Note that m_eject returned for ghosts by calc_ejected_mass() is
    // meaningless in the current physical prescriptions.  This fact is dealt
    // with in update_reservoirs_from_sn_feedback().

    assert(m_reheat >= 0);
    assert(m_eject >= 0);

    // update the baryonic reservoirs
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, m_recycled, new_metals);
}

static void backfill_ghost_NewStars(galaxy_t* gal, double m_stars, int snapshot)
{
    if ((snapshot - gal->LastIdentSnap) <= N_HISTORY_SNAPS) {
        double* LTTime = run_globals.LTTime;
        double burst_time = LTTime[gal->LastIdentSnap] - gal->dt * 0.5;

        for (int ii = 1; ii < N_HISTORY_SNAPS; ii++)
            if (LTTime[snapshot - ii] > burst_time) {
                gal->NewStars[ii - 1] += m_stars;
                break;
            }
    }
}

void contemporaneous_supernova_feedback(
    galaxy_t* gal,
    double* m_stars,
    int snapshot,
    double* m_reheat,
    double* m_eject,
    double* m_recycled,
    double* new_metals)
{
    // Here we approximate a constant SFR accross the timestep by a single burst
    // at t=0.5*dt. This is a pretty good approximation (to within ~15% of the
    // true number of SN that would have gone of by the end of the timestep for a
    // constant SFR). SN feedback due to merger driven starbursts adopts the same 
    // approximation.

    double sn_energy = 0.0;
    // init (just in case!)
    *m_reheat = *m_recycled = *new_metals = *m_eject = 0.0;

    // ** The IRA is broken in this version. **
    bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
    /* // N.B. If Flag_IRA is true then m_low and burst_recycled_frac will equal values in above declaration
    if (!Flag_IRA) {
        assert(snapshot > 0);
        double log_dt = log10(gal->dt * 0.5 *
                             run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h);
        double m_high = 120.0; // Msol
        double m_low = sn_m_low(log_dt); // Msol
        double burst_mass_frac;
        // calculate the mass reheated (from fraction of total SN-II that have gone off) from this burst
        burst_recycled_frac = calc_recycled_frac(m_high, m_low, &burst_mass_frac);
    }
    */

    // At this point, the baryonic reservoirs have not been updated. Thus, use the metallicity
    // of cold gas for new formed stars.
    double metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);
    // calculate recycled mass and metals by yield tables
    *m_recycled += *m_stars * get_yield(0, metallicity, Y_TOTAL);
    *new_metals = *m_stars * get_yield(0, metallicity, Y_TOTAL_METAL);

    // calculate the SNII energy and total reheated mass
    sn_energy = get_energy(0, metallicity) * *m_stars;
    *m_reheat = calc_sn_reheat_eff(gal, snapshot) * sn_energy/get_total_energy();
    sn_energy *= calc_sn_ejection_eff(gal, snapshot);

    // We can only reheat as much gas as we have available.  Let's inforce this
    // now, to ensure that the maximal amount of available energy is used to
    // eject gas from the system.
    if (*m_reheat > gal->ColdGas)
        *m_reheat = gal->ColdGas;

    // attenuate the star formation if necessary, so that we are being consistent
    // if (*m_reheat + *m_stars - *m_recycled > gal->ColdGas)
    if (*m_reheat + *m_stars > gal->ColdGas) {
        // double frac = gal->ColdGas / (*m_reheat + *m_stars - *m_recycled);
        double frac = gal->ColdGas / (*m_reheat + *m_stars);
        *m_reheat *= frac;
        *m_stars *= frac;
        // *m_recycled *= frac;
    }
    assert(*m_reheat >= 0);
    assert(*m_recycled >= 0);
    assert(*new_metals >= 0);

    // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
    *m_eject = calc_ejected_mass(m_reheat, sn_energy, gal->Vvir, gal->Halo->FOFGroup->Vvir);

    assert(*m_reheat >= 0);
    assert(*m_eject >= 0);

    // If this is a reidentified ghost, then back fill NewStars to reflect this
    // new SF burst.
    if (!Flag_IRA && (gal->LastIdentSnap < (snapshot - 1)))
        backfill_ghost_NewStars(gal, *m_stars, snapshot);
}
