#include "meraxes.h"
#include <assert.h>
#include <math.h>

#define ETA 0.06
//standard efficiency, 6% accreted mass is radiated

#define EMISSIVITY_CONVERTOR 8.40925088e-8
//=1e60 * PROTONMASS / 1e10 / SOLAR_MASS
// Unit_convertor is to convert emissivity from 1e60 photons to photons per atom

#define LUMINOSITY_CONVERTOR 32886.5934
// = (1 solar mass *(speed of light)^2/450e6year) /3.828e26watt
// where 450e6year is eddington time scale

#define LB2EMISSIVITY 67126822.0217
// this is a little bit complicated...
// firstly B band luminosity is:   LB = Lbol/kb 1e10Lsun
// then B band magnitude is:       MB = 4.74 - 2.5log10(1e10LB)
// then convert from Vega to AB:   MAB,B = MB-0.09
// then UV mag:                    M1450 = MAB,B+0.524
// then UV lum:                    LUV = 10**((M1450_1-51.594843)/-2.5) #erg/s/Hz
// then 912 lum:                   L912 = LUV*(1200./1450)**0.44*(912/1200.)**1.57  #erg/s/Hz
// then BHemissivity:              L912/Planck constant/1.57 #photons/s
// then total BHemissivity in this step: BHemissivity*accretion_time*SEC_PER_MEGAYEAR
// BHemissivity = fobs*Lbol/kb*2.1276330276278045e+54*SEC_PER_MEGAYEAR*accretion_time/1e60; // photon numbers/1e60
// BHemissivity = fobs*Lbol/kb*67126822.0217*accretion_time; // photon numbers/1e60

double calculate_BHemissivity(double BlackHoleMass, double accreted_mass)
{
    double accretion_time;
    double Lbol; //bolometric luminotisy
    double kb; //bolometric correction
    physics_params_t* physics = &(run_globals.params.physics);

    accretion_time = log1p(accreted_mass / BlackHoleMass) * 450.514890 * ETA / physics->EddingtonRatio; // Myr

    // Bolometric luminosity in 1e10 Lsun at the MIDDLE of accretion time
    // TODO: this introduce inconsistency compared to the calculation of luminosity.
    // Here we assume Lbol(t) = Lbol(t = accretion_time/2), which is only an approximation!
    Lbol = sqrt(1. + accreted_mass / BlackHoleMass) * physics->EddingtonRatio
        * BlackHoleMass / run_globals.params.Hubble_h * LUMINOSITY_CONVERTOR;
    kb = 6.25 * pow(Lbol, -0.37) + 9.0 * pow(Lbol, -0.012);

    return physics->quasar_fobs * Lbol / kb * LB2EMISSIVITY * accretion_time; //1e60 photon numbers
}

// quasar feedback suggested by Croton et al. 2016
void update_reservoirs_from_quasar_mode_bh_feedback(galaxy_t* gal, double m_reheat)
{
    double metallicity;
    galaxy_t* central;

    if (gal->ghost_flag)
        central = gal;
    else
        central = gal->Halo->FOFGroup->FirstOccupiedHalo->Galaxy;

    if (m_reheat < gal->ColdGas) {
        metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);
        gal->ColdGas -= m_reheat;
        gal->MetalsColdGas -= m_reheat * metallicity;
        central->MetalsHotGas += m_reheat * metallicity;
        central->HotGas += m_reheat;
    } else {
        metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);
        gal->ColdGas = 0.0;
        gal->MetalsColdGas = 0.0;
        central->HotGas -= m_reheat;
        central->MetalsHotGas -= m_reheat * metallicity;
        central->EjectedGas += m_reheat;
        central->MetalsEjectedGas += m_reheat * metallicity;
    }

    // Check the validity of the modified reservoir values (HotGas can be negtive for too strong quasar feedback)
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

double radio_mode_BH_heating(galaxy_t* gal, double cooling_mass, double x)
{
    double heated_mass = 0.0;

    // if there is any hot gas
    if (gal->HotGas > 0.0) {
        double Vvir;
        run_units_t* units = &(run_globals.units);

        if (gal->Type == 0)
            Vvir = gal->Halo->FOFGroup->Vvir;
        else
            Vvir = gal->Vvir;

        //bondi-hoyle accretion model
        double accreted_mass = run_globals.params.physics.RadioModeEff
            * run_globals.G * 3.4754 * x * gal->BlackHoleMass * gal->dt;
        // 15/8*pi*mu=3.4754, with mu=0.59; x=k*m_p*t/lambda

        // eddington rate
        double eddington_mass = exp(gal->dt * units->UnitTime_in_Megayears / run_globals.params.Hubble_h / 450.514890
                                     / ETA * run_globals.params.physics.EddingtonRatio) * gal->BlackHoleMass;

        // limit accretion by the eddington rate
        if (accreted_mass > eddington_mass)
            accreted_mass = eddington_mass;

        // limit accretion by amount of hot gas available
        if (accreted_mass > gal->HotGas)
            accreted_mass = gal->HotGas;

        // mass heated by AGN following Croton et al. 2006
        heated_mass = 2. * ETA * run_globals.Csquare / Vvir / Vvir * accreted_mass;

        // limit the amount of heating to the amount of cooling
        if (heated_mass > cooling_mass) {
            accreted_mass = cooling_mass / heated_mass * accreted_mass;
            heated_mass = cooling_mass;
        }

        gal->BlackHoleAccretedHotMass = accreted_mass;

        // add the accreted mass to the black hole
        double metallicity = calc_metallicity(gal->HotGas, gal->MetalsHotGas);

        //Assuming all energy from radio mode is going to heat the cooling flow
        //So no emissivity from radio mode!
        //TODO: we could add heating effienciency to split the energy into
        //heating and reionization.
        gal->BlackHoleMass += (1. - ETA) * accreted_mass;
        gal->HotGas -= accreted_mass;
        gal->MetalsHotGas -= metallicity * accreted_mass;
    }
    return heated_mass;
}

void merger_driven_BH_growth(galaxy_t* gal, double merger_ratio, int snapshot)
{
    if (gal->ColdGas > 0) {
        // If there is any cold gas to feed the black hole...
        double Vvir;

        // If this galaxy is the central of it's FOF group then use the FOF Halo properties
        // TODO: This needs closer thought as to if this is the best thing to do...
        if (gal->Type == 0)
            Vvir = gal->Halo->FOFGroup->Vvir;
        else
            Vvir = gal->Vvir;

        // Suggested by Bonoli et al. 2009 and Wyithe et al. 2003
        double z_scaling = pow((1 + run_globals.ZZ[snapshot]), run_globals.params.physics.quasar_mode_scaling);

        double accreting_mass = run_globals.params.physics.BlackHoleGrowthRate * merger_ratio / (1.0 + (280.0 * 280.0 / Vvir / Vvir)) * gal->ColdGas * z_scaling;

        // limit accretion to what is available
        if (accreting_mass > gal->ColdGas)
            accreting_mass = gal->ColdGas;

        // put the mass onto the accretion disk and let the black hole accrete it in the next snapshot
        // TODO: since the merger is put in the end of galaxy evolution, this is following the
        // inconsistence consistently
        gal->BlackHoleAccretingColdMass += accreting_mass;
    }
}

void previous_merger_driven_BH_growth(galaxy_t* gal)
{
    // If there is any cold gas to feed the black hole...
    double m_reheat;
    double accreted_mass;
    double BHemissivity;
    double Vvir;
    run_units_t* units = &(run_globals.units);

    // If this galaxy is the central of it's FOF group then use the FOF Halo properties
    // TODO: This needs closer thought as to if this is the best thing to do...
    if (gal->Type == 0)
        Vvir = gal->Halo->FOFGroup->Vvir;
    else
        Vvir = gal->Vvir;

    // Eddington rate
    accreted_mass = expm1(gal->dt * units->UnitTime_in_Megayears / run_globals.params.Hubble_h / 450.514890
                         / ETA * run_globals.params.physics.EddingtonRatio) * gal->BlackHoleMass;

    // limit accretion to what is need
    if (accreted_mass > gal->BlackHoleAccretingColdMass)
        accreted_mass = gal->BlackHoleAccretingColdMass;

    gal->BlackHoleAccretedColdMass += accreted_mass;
    gal->BlackHoleAccretingColdMass -= accreted_mass;

    //N_gamma,q * N_bh; later 1e60*BHemissivity * PROTONMASS/1e10/SOLAR_MASS will be N_gamma,q * M_bh
    BHemissivity = calculate_BHemissivity(gal->BlackHoleMass, accreted_mass);
    gal->BHemissivity += BHemissivity;
    gal->BlackHoleMass += (1. - ETA) * accreted_mass;
    gal->EffectiveBHM += BHemissivity * EMISSIVITY_CONVERTOR
        * run_globals.params.physics.ReionEscapeFracBH / run_globals.params.physics.ReionNionPhotPerBary;

    // quasar mode feedback
    m_reheat = run_globals.params.physics.QuasarModeEff * 2. * ETA * run_globals.Csquare * accreted_mass / Vvir / Vvir;
    update_reservoirs_from_quasar_mode_bh_feedback(gal, m_reheat);
}
