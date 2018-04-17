#include "meraxes.h"
#include <assert.h>
#include <math.h>

/*
 * Mostly follows Visbal, Haiman & Bryan 2018
 * IMF: all PopIII are M_POPIII
 * M_POPIII = 40 means all PopIII reaches SN within the snapshot
 *   no needs for delayed supernovae
 *
 * SN remnants increases BH seeds (still have BlackHoleSeed)
 *
 * Fo any questions, ask Yuxiang Qin
 *
 */

#define E_LW 2e-11 //erg
#define NU_LW 5.8e14 //erg

double calculate_J_21_LW(int snapshot){
    // TODO consider J_21_LW from nearby objects
    // only background is considered at this moment

    physics_params_t* params = &(run_globals.params.physics);
    run_units_t*       units = &(run_globals.units);
    galaxy_t*            gal = run_globals.FirstGal;

    double z1cube = powf(1.0 + run_globals.ZZ[snapshot], 3);
    double volume = powf(run_globals.params.BoxSize, 3);
    double factor = 1e60 * (PROTONMASS / units->UnitMass_in_g) *\
                    (1.0 - 0.75 * params->Y_He);

    double Nion = 0; //1e60 photons
    while (gal != NULL) {
        Nion += (gal->GrossStellarMass * params->ReionNionPhotPerBary +\
                 gal->PopIIIMass * params->ReionNionPhotPerBaryPopIII) /\
                 factor + gal->BHemissivity;
        gal = gal->Next;
    }

    /*Nion *= (C / units->UnitVelocity_in_cm_per_s) / (4. * M_PI) *\
            (E_LW / units.UnitEnergy_in_cgs) /\
            (NU_LW * units->UnitTime_in_s) / volume * z1cube /\
            (1e-21 /  units.UnitEnergy_in_cgs * pow(units->UnitLength_in_cm, 2.0));*/
    Nion *= C / (4. * M_PI) * E_LW / NU_LW / volume * z1cube * 1e21/ \
            pow(units->UnitLength_in_cm, 3);

    run_globals.J_21_LW_bg[snapshot] = Nion;
}

double calculate_Mvir_crit_LW(int snapshot){
    double redshift      = run_globals.ZZ[snapshot];
    double term_J_21_LW  = pow(4 * M_PI * run_globals.J_21_LW_bg[snapshot], 0.47);
    double term_redshift = pow((1. + redshift)/26., -1.5);
    double Mvir_crit_LW  = 3.73e-5 * term_redshift * (1. + 6.96 * term_J_21_LW);
    double Matomic       = Tvir_to_Mvir(1e4, redshift);
	return Mvir_crit_LW >= Matomic ? Mvir_crit_LW : Matomic;
}

void popiii_supernova_feedback(galaxy_t* gal, double m_popiii){
    // all popiii reach SN

    double* m_reheat;
	double sn_energy, m_eject, new_metals;
    double SnReheatScaling = run_globals.params.physics.SnReheatScaling;
    double SnReheatNorm = run_globals.params.physics.SnReheatNorm;
    double SnReheatEff = run_globals.params.physics.SnReheatEff;
    double SnReheatLimit = run_globals.params.physics.SnReheatLimit;

    // scale the reheating efficiency
    if (SnReheatScaling != 0)
        SnReheatEff *= 0.5 + pow(gal->Vmax / SnReheatNorm, -SnReheatScaling);
    if (SnReheatEff > SnReheatLimit)
        SnReheatEff = SnReheatLimit;

    new_metals    = m_popiii * run_globals.params.physics.Yield;
	*m_reheat     = m_popiii * SnReheatEff;
    if (*m_reheat > gal->ColdGas)
        *m_reheat = gal->ColdGas;

    sn_energy = calc_sn_energy(m_popiii, gal->Vmax, 1.0);
    m_eject   = calc_ejected_mass(m_reheat, sn_energy, gal->Vvir, gal->Halo->FOFGroup->Vvir);

    // all popiii remnants added to blackhole
    // TODO: a fraction can be used for feedback
    // a fraction can be recycled to coldgas
    // the remnants can be used as BlackHoleSeed
    gal->BlackHoleMass += m_popiii;

    assert(*m_reheat >= 0);
    assert(m_eject >= 0);

    update_reservoirs_from_sn_feedback(gal, *m_reheat, m_eject, 0.0, new_metals);
}

void popiii_star_formation(galaxy_t* gal, int snapshot){

    int    n_popiii;
    double m_popiii;
    double M_POPIII = run_globals.params.physics.M_POPIII;

    m_popiii  = run_globals.params.physics.PopIIIEfficiency *\
                (gal->ColdGas + gal->HotGas);
    n_popiii  = (int)(m_popiii / M_POPIII);
    n_popiii += ((double)rand() / RAND_MAX) > (m_popiii - n_popiii * M_POPIII);
    m_popiii  = n_popiii * M_POPIII;

    // TODO: going to starve popiii SF when molecular cooling is implemented.
    // now galaxy form popiii as long as cold+hot is enough
    // no matter how small the cold gas is...

    if(m_popiii < gal->ColdGas)
        gal->ColdGas -= m_popiii;
    else{
        gal->ColdGas  = 0.0;
        gal->HotGas  -= (m_popiii - gal->ColdGas);
    }
    assert(gal->ColdGas >= 0);
    assert(gal->HotGas >= 0);
    gal->PopIIIMass += m_popiii;
    gal->Sfr        += m_popiii / gal->dt;

    popiii_supernova_feedback(gal, m_popiii);
}

void evolve_PopIII(galaxy_t* gal, int snapshot){
    double Mvir_crit_LW = calculate_Mvir_crit_LW(snapshot);
    if(gal->Mvir > Mvir_crit_LW)
        popiii_star_formation(gal, snapshot);
}
