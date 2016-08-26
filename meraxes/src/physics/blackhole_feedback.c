#include "meraxes.h"
#include <math.h>
#include <assert.h> 

#define eta 0.06 //standard efficiency, 6% accreted mass is radiated 

double calculate_emissivity(run_globals_t *run_globals, double BlackHoleMass, double accreted_mass)
{
    double accretion_time;
    double Lbol;//bolometric luminotisy
    double kb;//bolometric correction
    run_units_t *units = &(run_globals->units);

    accretion_time = log10(1.0+accreted_mass/BlackHoleMass)/run_globals->params.physics.EddingtonRatio*eta / (1.402e37 / (units->UnitEnergy_in_cgs / units->UnitTime_in_Megayears)); //in second

    // this will be bolometric luminosity in 1e10Lsun, just google this
    // (1 solar mass *(speed of light)^2/450e6year) /3.828e26watt
    // where 450e6year is eddington time scale
    Lbol = run_globals->params.physics.EddingtonRatio *BlackHoleMass /run_globals->params.Hubble_h * 32886.5934;
    kb = 6.25*pow(Lbol, -0.37) + 9.0*pow(Lbol,-0.012);
    // this is very complicated...
    // firstly B band luminosity is:   LB = Lbol/kb 1e10Lsun
    // then B band magnitude is:       MB = 4.74 - 2.5log10(1e10LB)
    // then convert from Vega to AB:   MAB,B = MB-0.09
    // then UV mag:                    M1450 = MAB,B+0.524
    // then UV lum:                    LUV = 10**((M1450_1-51.594843)/-2.5) #erg/s/Hz
    // then 912 lum:                   L912 = LUV*(1200./1450)**0.44*(912/1200.)**1.57  #erg/s/Hz
    // then emissivity:                L912/Planck constant/1.57 #photons/s 
    // then total emissivity in this step: emissivity*accretion_time
    // emissivity = Lbol/kb*2.1276330276278045e+54*accretion_time; //photon numbers

    return Lbol/kb*2.1276330276278045e+54*accretion_time; //photon numbers 
}

// quasar feedback suggested by Croton et al. 2016
void update_reservoirs_from_quasar_mode_bh_feedback(run_globals_t *run_globals, galaxy_t *gal, double m_reheat)
{
    double metallicity; 
    galaxy_t *central;

    if (gal->ghost_flag)
        central = gal; 
    else
        central = gal->Halo->FOFGroup->FirstOccupiedHalo->Galaxy;   

    if (m_reheat < gal->ColdGas)
    {
        metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);
        gal->ColdGas          -= m_reheat;
        gal->MetalsColdGas    -= m_reheat * metallicity; 
        central->MetalsHotGas += m_reheat * metallicity;
        central->HotGas       += m_reheat;              
    }
    else
    {
        metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);
        gal->ColdGas               = 0.0;
        gal->MetalsColdGas         = 0.0;
        central->HotGas           -= m_reheat;
        central->MetalsHotGas     -= m_reheat * metallicity;
        central->EjectedGas       += m_reheat;              
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

double radio_mode_BH_heating(run_globals_t *run_globals, galaxy_t *gal, double cooling_mass, double x)
{
  double eddington_mass;
  double accreted_mass;
  double heated_mass;
  double metallicity;
  double emissivity;

  fof_group_t *fof_group = gal->Halo->FOFGroup;

  run_units_t *units = &(run_globals->units);

  // if there is any hot gas
  if (gal->HotGas > 0.0)
  {

    //Bondi-Hoyle accretion model
    accreted_mass = run_globals->params.physics.RadioModeEff
                    * run_globals->G * 1.7377 * x * gal->BlackHoleMass*gal->dt;
    // 15/16*pi*mu=1.7377, with mu=0.59; x=k*m_p*T/Lambda

    // Eddington rate
    eddington_mass = (exp(1.402e37 / (units->UnitEnergy_in_cgs / units->UnitTime_in_s)*gal->dt/eta*run_globals->params.physics.EddingtonRatio)-1.) * gal->BlackHoleMass;

    // limit accretion by the eddington rate
    if (accreted_mass > eddington_mass)
      accreted_mass = eddington_mass;

    // limit accretion by amount of hot gas available
    if (accreted_mass > gal->HotGas)
      accreted_mass = gal->HotGas;

    gal->BlackHoleAccretedHotMass = accreted_mass;

    // mass heated by AGN following Croton et al. 2006
    heated_mass = 2.*eta*8.98755e10 / fof_group->Vvir / fof_group->Vvir * accreted_mass;

    // limit the amount of heating to the amount of cooling
    if (heated_mass > cooling_mass)
    {
      accreted_mass = cooling_mass / heated_mass * accreted_mass;
      heated_mass   = cooling_mass;
    }

    // add the accreted mass to the black hole
    metallicity         = calc_metallicity(gal->HotGas, gal->MetalsHotGas);

    //N_gamma,q * N_bh; later emissivity * PROTONMASS/1e10/SOLAR_MASS will be  N_gamma,q * M_bh
    emissivity = calculate_emissivity(run_globals, gal->BlackHoleMass, accreted_mass);

    gal->emissivity = emissivity;
    gal->BlackHoleMass      += (1.-eta)*accreted_mass;
    gal->BlackHoleGrossMass += (1.-eta)*accreted_mass;
    gal->EffectiveBHM       += emissivity*PROTONMASS/1e10/SOLAR_MASS*run_globals->params.physics.ReionEscapeFracBH/run_globals->params.physics.ReionNionPhotPerBary/run_globals->params.physics.ReionEscapeFrac;
    gal->FescWeightedEBHM   += emissivity*PROTONMASS/1e10/SOLAR_MASS*run_globals->params.physics.ReionEscapeFracBH/run_globals->params.physics.ReionNionPhotPerBary;
    gal->HotGas        -= accreted_mass;
    gal->MetalsHotGas  -= metallicity * accreted_mass;
  }
  else
  {
    // if there is no hot gas
    gal->BlackHoleAccretedHotMass = 0.0;
    gal->emissivity = 0.0;
    heated_mass = 0.0;
  }

  return heated_mass;
}


void merger_driven_BH_growth(run_globals_t *run_globals, galaxy_t *gal, double merger_ratio, int snapshot)
{
  if (gal->ColdGas > 0)
  {
    // If there is any cold gas to feed the black hole...
    double m_reheat;
    double accreted_mass;
    double accreted_metals;
    double emissivity;
    double Vvir;
    double zplus1to1pt5;
    run_units_t *units = &(run_globals->units);

    // If this galaxy is the central of it's FOF group then use the FOF halo properties
    // TODO: This needs closer thought as to if this is the best thing to do...
    if (gal->Type == 0)
      Vvir = gal->Halo->FOFGroup->Vvir;
    else
      Vvir = gal->Vvir;

    // Suggested by Bonoli et al. 2009 and Wyithe et al. 2003
    zplus1to1pt5 = pow((1 + run_globals->ZZ[snapshot]), 1.5);

    assert(gal->BlackHoleAccretingColdMass >=0);
    gal->BlackHoleAccretingColdMass += run_globals->params.physics.BlackHoleGrowthRate * merger_ratio /
                    (1.0 + (280.0 * 280.0 / Vvir / Vvir)) * gal->ColdGas* zplus1to1pt5;

    // Eddington rate
    accreted_mass = (exp(1.402e37 / (units->UnitEnergy_in_cgs / units->UnitTime_in_s)*gal->dt/eta*run_globals->params.physics.EddingtonRatio)-1.) * gal->BlackHoleMass;

    // limit accretion to what is need
    if (accreted_mass > gal->BlackHoleAccretingColdMass)
      accreted_mass = gal->BlackHoleAccretingColdMass;

    // limit accretion to what is available
    if (accreted_mass > gal->ColdGas)
      accreted_mass = gal->ColdGas;

    gal->BlackHoleAccretedColdMass = accreted_mass;
    gal->BlackHoleAccretingColdMass -= accreted_mass;
    accreted_metals     = calc_metallicity(gal->ColdGas, gal->MetalsColdGas) * accreted_mass;

    //N_gamma,q * N_bh; later emissivity * PROTONMASS/1e10/SOLAR_MASS will be  N_gamma,q * M_bh
    emissivity = calculate_emissivity(run_globals, gal->BlackHoleMass, accreted_mass); 
    
    gal->emissivity +=emissivity;
    gal->BlackHoleMass      += (1.-eta)*accreted_mass;
    gal->BlackHoleGrossMass += (1.-eta)*accreted_mass; 
    gal->EffectiveBHM       += emissivity*PROTONMASS/1e10/SOLAR_MASS*run_globals->params.physics.ReionEscapeFracBH/run_globals->params.physics.ReionNionPhotPerBary/run_globals->params.physics.ReionEscapeFrac;
    gal->FescWeightedEBHM   += emissivity*PROTONMASS/1e10/SOLAR_MASS*run_globals->params.physics.ReionEscapeFracBH/run_globals->params.physics.ReionNionPhotPerBary;
    gal->ColdGas       -= accreted_mass;
    gal->MetalsColdGas -= accreted_metals;

    m_reheat = run_globals->params.physics.QuasarModeEff *eta*8.98755e10 * accreted_mass /Vvir /Vvir;
    update_reservoirs_from_quasar_mode_bh_feedback(run_globals, gal, m_reheat);
  }
  else
  {
    gal->BlackHoleAccretedColdMass = 0.0;
  }
}

void previous_merger_driven_BH_growth(run_globals_t *run_globals, galaxy_t *gal)
{
  if ((gal->ColdGas > 0) && (gal->BlackHoleAccretingColdMass >0))
  {
    // If there is any cold gas to feed the black hole...
    double m_reheat;
    double accreted_mass;
    double accreted_metals;
    double emissivity;
    double Vvir;
    run_units_t *units = &(run_globals->units);

    // If this galaxy is the central of it's FOF group then use the FOF halo properties
    // TODO: This needs closer thought as to if this is the best thing to do...
    if (gal->Type == 0)
      Vvir = gal->Halo->FOFGroup->Vvir;
    else
      Vvir = gal->Vvir;

    // Eddington rate
    accreted_mass = (exp(1.402e37 / (units->UnitEnergy_in_cgs / units->UnitTime_in_s)*gal->dt/eta*run_globals->params.physics.EddingtonRatio)-1.) * gal->BlackHoleMass;

    // limit accretion to what is need
    if (accreted_mass > gal->BlackHoleAccretingColdMass)
      accreted_mass = gal->BlackHoleAccretingColdMass;

    // limit accretion to what is available
    if (accreted_mass > gal->ColdGas)
      accreted_mass = gal->ColdGas;

    gal->BlackHoleAccretedColdMass = accreted_mass;
    gal->BlackHoleAccretingColdMass -= accreted_mass;
    accreted_metals     = calc_metallicity(gal->ColdGas, gal->MetalsColdGas) * accreted_mass;

    //N_gamma,q * N_bh; later emissivity * PROTONMASS/1e10/SOLAR_MASS will be  N_gamma,q * M_bh
    emissivity = calculate_emissivity(run_globals, gal->BlackHoleMass, accreted_mass);

    gal->emissivity +=emissivity;
    gal->BlackHoleMass      += (1.-eta)*accreted_mass;
    gal->BlackHoleGrossMass += (1.-eta)*accreted_mass;
    gal->EffectiveBHM       += emissivity*PROTONMASS/1e10/SOLAR_MASS*run_globals->params.physics.ReionEscapeFracBH/run_globals->params.physics.ReionNionPhotPerBary/run_globals->params.physics.ReionEscapeFrac;
    gal->FescWeightedEBHM   += emissivity*PROTONMASS/1e10/SOLAR_MASS*run_globals->params.physics.ReionEscapeFracBH/run_globals->params.physics.ReionNionPhotPerBary;
    gal->ColdGas       -= accreted_mass;
    gal->MetalsColdGas -= accreted_metals;

    m_reheat = run_globals->params.physics.QuasarModeEff *eta*8.98755e10 * accreted_mass /Vvir /Vvir;
    update_reservoirs_from_quasar_mode_bh_feedback(run_globals, gal, m_reheat);
  }
  else
  {
    gal->BlackHoleAccretedColdMass = 0.0;
  }
}
