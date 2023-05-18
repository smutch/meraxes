#include <assert.h>
#include <math.h>

#include "core/misc_tools.h"
#include "core/stellar_feedback.h"
#include "core/PopIII.h"
#include "core/virial_properties.h"
#include "meraxes.h"
#include "supernova_feedback.h"

void update_reservoirs_from_sn_feedback(galaxy_t* gal,
                                        double m_reheat,
                                        double m_eject,
                                        double m_recycled,
                                        double m_remnant,
                                        double new_metals)
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

  gal->StellarMass -= (m_recycled + m_remnant);
  if (gal->Galaxy_Population == 2)
    gal->StellarMass_II -= m_recycled;
  else if (gal->Galaxy_Population == 3) {
    gal->StellarMass_III -= (m_recycled + m_remnant);
    gal->Remnant_Mass += m_remnant;
    }
  // N.B. Stellar metallicity does not work properly. Metals are generated by
  // nuclear reaction in stars, and modelling this implicit evolution is tricky.
  // Stellar metallicity does not influence galaxy evolution, and shoud not use
  // properties.
  gal->MetalsStellarMass -= new_metals;
  gal->ColdGas += m_recycled;

  // assuming instantaneous recycling approximation and enrichment from SNII
  // only, work out the mass of metals returned to the ISM by this SF burst
  if (gal->ColdGas > 1e-10)
    gal->MetalsColdGas += new_metals;
  else
    central->MetalsHotGas += new_metals;

  // make sure we aren't trying to use more cold gas than is available...
  if (m_reheat > gal->ColdGas){
    m_reheat = gal->ColdGas;
    }
  metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

  gal->ColdGas -= m_reheat;
  gal->MetalsColdGas -= m_reheat * metallicity;
  central->MetalsHotGas += m_reheat * metallicity;
  central->HotGas += m_reheat;

  // If this is a ghost then we don't know what the real ejected mass is as we
  // don't know the properties of the halo!
  if (!gal->ghost_flag) {
    metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);

    if (m_eject > central->HotGas){
      m_eject = central->HotGas;
      //m_eject_III = central->HotGas - m_eject_II;
      }

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
  if (gal->StellarMass_II < 0)
    gal->StellarMass_II = 0.0;
  if (gal->StellarMass_III < 0)
    gal->StellarMass_III = 0.0;
  if (gal->Remnant_Mass < 0)
    gal->Remnant_Mass = 0.0;
  if (gal->MetalsStellarMass < 0)
    gal->MetalsStellarMass = 0.0;
  if (central->EjectedGas < 0)
    central->EjectedGas = 0.0;
  if (central->MetalsEjectedGas < 0)
    central->MetalsEjectedGas = 0.0;
}

void update_reservoirs_from_delayed_sn_feedback(galaxy_t* gal, //You don't need reheat_III/II and eject III/II but you might want to keep those to check conditions!
                                        double m_reheat,
                                        double m_eject,
                                        double m_recycled,
                                        double m_recycled_III,
                                        double m_recycled_II,
                                        double m_remnant,
                                        double new_metals)
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

  gal->StellarMass -= (m_recycled + m_remnant);
  gal->StellarMass_II -= m_recycled_II;
  gal->StellarMass_III -= (m_recycled_III + m_remnant);
  gal->Remnant_Mass += m_remnant;
  // N.B. Stellar metallicity does not work properly. Metals are generated by
  // nuclear reaction in stars, and modelling this implicit evolution is tricky.
  // Stellar metallicity does not influence galaxy evolution, and shoud not use
  // properties.
  gal->MetalsStellarMass -= new_metals;
  gal->ColdGas += m_recycled;

  // assuming instantaneous recycling approximation and enrichment from SNII
  // only, work out the mass of metals returned to the ISM by this SF burst
  if (gal->ColdGas > 1e-10)
    gal->MetalsColdGas += new_metals;
  else
    central->MetalsHotGas += new_metals;

  // make sure we aren't trying to use more cold gas than is available...
  if (m_reheat > gal->ColdGas){
    m_reheat = gal->ColdGas;
    }
  metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

  gal->ColdGas -= m_reheat;
  gal->MetalsColdGas -= m_reheat * metallicity;
  central->MetalsHotGas += m_reheat * metallicity;
  central->HotGas += m_reheat;

  // If this is a ghost then we don't know what the real ejected mass is as we
  // don't know the properties of the halo!
  if (!gal->ghost_flag) {
    metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);

    if (m_eject > central->HotGas){
      m_eject = central->HotGas;
      }

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
  if (gal->StellarMass_II < 0)
    gal->StellarMass_II = 0.0;
  if (gal->StellarMass_III < 0)
    gal->StellarMass_III = 0.0;
  if (gal->MetalsStellarMass < 0)
    gal->MetalsStellarMass = 0.0;
  if (central->EjectedGas < 0)
    central->EjectedGas = 0.0;
  if (central->MetalsEjectedGas < 0)
    central->MetalsEjectedGas = 0.0;
}

static inline double calc_ejected_mass(double* m_reheat, double sn_energy, double Vvir, double fof_Vvir)
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
    } else if (fof_Vvir > 0) {
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

static inline double calc_sn_reheat_eff(galaxy_t* gal, int snapshot, int flag_population)
{
  double Vmax = gal->Vmax; // Vmax is in a unit of km/s
  double zplus1 = 1. + run_globals.ZZ[snapshot];
  physics_params_t* params = &run_globals.params.physics;
  int SnModel = params->SnModel;
  double SnReheatRedshiftDep;
  double SnReheatEff;
  double SnReheatScaling;
  double SnReheatNorm;
  double SnReheatLimit;
  if (flag_population == 2) {
    SnReheatRedshiftDep = params->SnReheatRedshiftDep;
    SnReheatEff = params->SnReheatEff;
    SnReheatScaling = params->SnReheatScaling;
    SnReheatNorm = params->SnReheatNorm;
    SnReheatLimit = params->SnReheatLimit;
    }
  else if (flag_population == 3) {
    SnReheatRedshiftDep = params->SnReheatRedshiftDep_III;
    SnReheatEff = params->SnReheatEff_III;
    SnReheatScaling = params->SnReheatScaling_III;
    SnReheatNorm = params->SnReheatNorm_III;
    SnReheatLimit = params->SnReheatLimit_III;
    }
  switch (SnModel) {
    case 1: // Guo et al. 2011 with redshift dependence
      SnReheatEff *= pow(zplus1 / 4., SnReheatRedshiftDep) * (.5 + pow(Vmax / SnReheatNorm, -SnReheatScaling));
      break;
    case 2: // Muratov et al. 2015
      if (Vmax < SnReheatNorm)
        SnReheatScaling = params->SnReheatScaling2;
      SnReheatEff *= pow(zplus1 / 4., SnReheatRedshiftDep) * pow(Vmax / SnReheatNorm, -SnReheatScaling);
      break;
    default:
      mlog_error("Unknonw SnModel!");
      ABORT(EXIT_FAILURE);
      break;
  }
  if (SnReheatEff < SnReheatLimit)
    return SnReheatEff;
  else
    return SnReheatLimit;
}

static inline double calc_sn_ejection_eff(galaxy_t* gal, int snapshot, int flag_population)
{
  double Vmax = gal->Vmax;    // Vmax is in a unit of km/s
  double zplus1 = 1. + run_globals.ZZ[snapshot];
  physics_params_t *params = &run_globals.params.physics;
  int SnModel = params->SnModel;
  double SnEjectionRedshiftDep;
  double SnEjectionEff;
  double SnEjectionScaling;
  double SnEjectionNorm;
  if (flag_population == 2) {
    SnEjectionRedshiftDep = params->SnEjectionRedshiftDep;
    SnEjectionEff = params->SnEjectionEff;
    SnEjectionScaling = params->SnEjectionScaling;
    SnEjectionNorm = params->SnEjectionNorm;
    }
  else if (flag_population == 3) {
    SnEjectionRedshiftDep = params->SnEjectionRedshiftDep_III;
    SnEjectionEff = params->SnEjectionEff_III;
    SnEjectionScaling = params->SnEjectionScaling_III;
    SnEjectionNorm = params->SnEjectionNorm_III;
    }
  switch (SnModel) {
  case 1:    // Guo et al. 2011 with redshift dependence
      SnEjectionEff *= pow(zplus1/4., SnEjectionRedshiftDep) \
                       *(.5 + pow(Vmax/SnEjectionNorm, -SnEjectionScaling));
      break;
  case 2:
      // Use the same value with that is used for the mass loading
      if (Vmax < SnEjectionNorm)
        SnEjectionScaling = params->SnEjectionScaling2; 
        SnEjectionEff *= pow(zplus1/4., SnEjectionRedshiftDep) \
                         *pow(Vmax/SnEjectionNorm, -SnEjectionScaling);
        break;
  default:
    mlog_error("Unknonw SnModel!");
    ABORT(EXIT_FAILURE);
    break;
  }
  if (SnEjectionEff < 1.)
    return SnEjectionEff;
  else
    return 1.;
}

void delayed_supernova_feedback(galaxy_t* gal, int snapshot) 
{
  double sn_energy = 0.0;
  double sn_energy_II = 0.0;
  double sn_energy_III = 0.0;
  double m_reheat = 0.0;
  double m_reheat_II = 0.0;
  double m_reheat_III = 0.0;
  double m_eject = 0.0;
  double m_eject_II = 0.0;
  double m_eject_III = 0.0;
  double m_recycled = 0.0;
  double m_recycled_II = 0.0;
  double m_recycled_III = 0.0;
  double new_metals = 0.0;
  double m_remnant = 0.0;
  double fof_Vvir;
  
  double NumberPISN = run_globals.NumberPISN;
  double MassPISN = run_globals.MassPISN;
  double NumberSNII = run_globals.NumberSNII;
  double MassSNII = run_globals.MassSNII;
  double MassBHs = run_globals.MassBHs;
  
  double energy_unit = run_globals.units.UnitEnergy_in_cgs; 
  // If we are at snapshot < N_HISTORY_SNAPS-1 then only try to look back to snapshot 0
  int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;

  // Loop through each of the last `N_HISTORY_SNAPS` recorded stellar mass
  // bursts and calculate the amount of energy and mass that they will release
  // in the current time step.
  for (int i_burst = 1; i_burst < n_bursts; i_burst++) {
    double m_stars_II = gal->NewStars_II[i_burst];
    double m_stars_III = gal->NewStars_III[i_burst];
    double m_stars = m_stars_II + m_stars_III;

    // Only need to do this if any stars formed in this history bin
    if (m_stars > 1e-10) {
      double metallicity = calc_metallicity(m_stars, gal->NewMetals[i_burst]);
      // Calculate recycled mass and metals by yield tables
      
      m_recycled_II += m_stars_II * get_recycling_fraction(i_burst, metallicity);
      m_recycled_III += m_stars_III * CCSN_PopIII_Yield(i_burst, snapshot, 0); // Only CCSN have delayed feedback
      m_recycled += (m_recycled_II + m_recycled_III);
      new_metals += (m_stars_II * get_metal_yield(i_burst, metallicity) + m_stars_III * CCSN_PopIII_Yield(i_burst, snapshot, 1));
      m_remnant += m_stars_III * CCSN_PopIII_Yield(i_burst, snapshot, 2); // Pop. III will have remnants
      // Calculate SNII energy
      
      sn_energy_II += get_SN_energy(i_burst, metallicity) * m_stars_II; 
      sn_energy_III += get_SN_energy_PopIII(i_burst, snapshot, 0) * m_stars_III; // Tis is DeltaM reheat (eq.16 Mutch+16) * ENOVA
    }
  }

  m_reheat_II = calc_sn_reheat_eff(gal, snapshot, 2) * sn_energy_II / get_total_SN_energy();
  sn_energy_II *= calc_sn_ejection_eff(gal, snapshot, 2);
  m_reheat_III = calc_sn_reheat_eff(gal, snapshot, 3) * sn_energy_III / ENOVA_CC; 
  sn_energy_III *= (calc_sn_ejection_eff(gal, snapshot, 3) * NumberSNII * 1e10 / run_globals.params.Hubble_h); //Maybe for the SN ejection efficiency is more important to distinguish between PISN/CC rather than Pop.III/II 
  m_reheat = m_reheat_II + m_reheat_III;
  sn_energy = sn_energy_II + sn_energy_III / energy_unit; //Convert from erg to internal units! 10^10Msol/h * (km/s)^2
  
  // We can only reheat as much gas as we have available.  Let's inforce this
  // now, to ensure that the maximal amount of available energy is used to
  // eject gas from the system.
  if (m_reheat != m_reheat_III + m_reheat_II)
    m_reheat = m_reheat_III + m_reheat_II;
  if (m_reheat > gal->ColdGas)
    m_reheat = gal->ColdGas;

  assert(m_reheat >= 0);
  assert(m_recycled >= 0);
  assert(new_metals >= 0);
  assert(m_reheat_III >= 0);
  assert(m_recycled_III >= 0);
  assert(m_reheat_II >= 0);
  assert(m_recycled_II >= 0);
  assert(m_remnant >= 0);
  
  // how much mass is ejected due to this star formation episode?
  if (!gal->ghost_flag)
    fof_Vvir = gal->Halo->FOFGroup->Vvir;
  else
    fof_Vvir = -1;

  m_eject_III = calc_ejected_mass(&m_reheat_III, sn_energy_III / energy_unit, gal->Vvir, fof_Vvir);
  m_eject_II = calc_ejected_mass(&m_reheat_II, sn_energy_II, gal->Vvir, fof_Vvir);
  m_eject = m_eject_II + m_eject_III;

  // Note that m_eject returned for ghosts by calc_ejected_mass() is
  // meaningless in the current physical prescriptions.  This fact is dealt
  // with in update_reservoirs_from_sn_feedback().

  assert(m_reheat >= 0);
  assert(m_eject >= 0);
  assert(m_reheat_III >= 0);
  assert(m_eject_III >= 0);
  assert(m_reheat_II >= 0);
  assert(m_eject_II >= 0);
  
  if (m_recycled_II + m_recycled_III != m_recycled)
    m_recycled = m_recycled_II + m_recycled_III;
    
  if (m_eject != m_eject_II + m_eject_III)  
    m_eject = m_eject_II + m_eject_III;

  // update the baryonic reservoirs
  update_reservoirs_from_delayed_sn_feedback(gal, m_reheat, m_eject, m_recycled, m_recycled_III, m_recycled_II, m_remnant, new_metals);
}

void contemporaneous_supernova_feedback(galaxy_t* gal,
                                        double* m_stars,
                                        int snapshot,
                                        double* m_reheat,
                                        double* m_eject,
                                        double* m_recycled,
                                        double* m_remnant,
                                        double* new_metals)
{
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
  double sn_energy = 0.0;
  double energy_unit = run_globals.units.UnitEnergy_in_cgs;
  double NumberPISN = run_globals.NumberPISN;
  double MassPISN = run_globals.MassPISN;
  double NumberSNII = run_globals.NumberSNII;
  double MassSNII = run_globals.MassSNII;
  double MassBHs = run_globals.MassBHs;
  
  // init (just in case!)
  *m_reheat = *m_recycled = *new_metals = *m_eject = *m_remnant = 0.0;

  // Here we approximate a constant SFR accross the timestep by a single burst
  // at t=0.5*dt. This is a pretty good approximation (to within ~15% of the
  // true number of SN that would have gone of by the end of the timestep for a
  // constant SFR). SN feedback due to merger driven starbursts adopts the same
  // approximation.

  // At this point, the baryonic reservoirs have not been updated. Thus, use the metallicity
  // of cold gas for new formed stars.
  double metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);
  if (!Flag_IRA) {
    // Calculate recycled mass and metals by yield tables
    // Total yield includes H and He and all other elements
    // Total metal yield includes all elements except H and He
    if (gal->Galaxy_Population == 2){
      *m_recycled = *m_stars * get_recycling_fraction(0, metallicity);
      *new_metals = *m_stars * get_metal_yield(0, metallicity);
      }
    else if (gal->Galaxy_Population == 3){
      *m_recycled = *m_stars * (CCSN_PopIII_Yield(0, snapshot, 0) + PISN_PopIII_Yield(0));
      *m_remnant = *m_stars * (MassBHs + CCSN_PopIII_Yield(0, snapshot, 2));
      *new_metals = *m_stars * CCSN_PopIII_Yield(0, snapshot, 1);
      if (MassPISN > 0) // Account for PISN
        *new_metals += *m_stars * PISN_PopIII_Yield(1) - (20.0 / 1e10 * run_globals.params.Hubble_h);
      }
  } else {
    // Recycling fraction and metals yield are input parameters when using IRA
    if (gal->Galaxy_Population == 2){
      *m_recycled = *m_stars * run_globals.params.physics.SfRecycleFraction;
      *new_metals = *m_stars * run_globals.params.physics.Yield;
      }
    else if (gal->Galaxy_Population == 3){
      *m_recycled = *m_stars * run_globals.params.physics.SfRecycleFraction_III;
      *new_metals = *m_stars * run_globals.params.physics.Yield_III;
    }
  }
  
  if (gal->Galaxy_Population == 2){
  // calculate the SNII energy and total reheated mass
    sn_energy = *m_stars * get_SN_energy(0, metallicity);
    *m_reheat = calc_sn_reheat_eff(gal, snapshot, 2) * sn_energy / get_total_SN_energy();
    sn_energy *= calc_sn_ejection_eff(gal, snapshot, 2);
    }
  else if (gal->Galaxy_Population == 3){ 
    sn_energy = (*m_stars) * (get_SN_energy_PopIII(0, snapshot, 0) + get_SN_energy_PopIII(0, snapshot, 1)); //erg
    sn_energy /= energy_unit; //Convert this because you need in internal units it for m_ejected
    *m_reheat = calc_sn_reheat_eff(gal, snapshot, 3) * (*m_stars) * (get_SN_energy_PopIII(0, snapshot, 1) / ENOVA_PISN + get_SN_energy_PopIII(0, snapshot, 0) / ENOVA_CC );
    sn_energy *= calc_sn_ejection_eff(gal, snapshot, 3);  
    }

  // We can only reheat as much gas as we have available.  Let's inforce this
  // now, to ensure that the maximal amount of available energy is used to
  // eject gas from the system.
  if (*m_reheat > gal->ColdGas)
    *m_reheat = gal->ColdGas;
    
  // You might add here a condition for m_remnant!
    
  // attenuate the star formation if necessary, so that we are being consistent
  // if (*m_reheat + *m_stars - *m_recycled > gal->ColdGas)
  if (*m_reheat + *m_stars > gal->ColdGas) {
    double frac = gal->ColdGas / (*m_reheat + *m_stars);
    
    *m_reheat *= frac;
    *m_stars *= frac;
    *m_recycled *= frac;
    *m_remnant *= frac;
  }
  if (*new_metals < 0) // Just to be sure
    *new_metals = 0.0;
  assert(*m_recycled >= 0);
  assert(*m_reheat >= 0);
  assert(*m_remnant >= 0);

  // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
  *m_eject = calc_ejected_mass(m_reheat, sn_energy, gal->Vvir, gal->Halo->FOFGroup->Vvir);
  
  assert(*m_reheat >= 0);
  assert(*m_eject >= 0);
  
}

/*void calc_metal_bubble(galaxy_t* gal, int snapshot) // Result in internal units (Mpc/h) You need to update this function for Pop III/Pop II, but before check it with Yuxiang
{
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
  double mm_stars = gal->NewStars[0]; //The last episode of SF
  
  double UnitMass_in_g = run_globals.units.UnitMass_in_g;
  double UnitLength_in_cm = run_globals.units.UnitLength_in_cm;
  double time_unit = run_globals.units.UnitTime_in_s;
  
  int A = gal->count_SF;
  
  if (Flag_IRA == false) {
  
    if (mm_stars > 1e-10) { 
      if (gal->Galaxy_Population == 3) //Crucial to update the galaxy index! 
        gal->Galaxy_Population = 2;

      gal->count_SF += 1;
      double gas_density;
      
      if (gal->count_SF > 70)
        mlog_error("Too many SF episodes"); 
      gas_density = (gal->HotGas + gal->ColdGas) * UnitMass_in_g / PROTONMASS / (4.0 * M_PI / 3.0 * pow(gal->Rvir * UnitLength_in_cm, 3.)); // cm^-3
    
      gal->Prefactor[A] = pow(EnergySN * N_SN_Pop2 * mm_stars * UnitMass_in_g / SOLAR_MASS / (PROTONMASS * gas_density), 0.2) / UnitLength_in_cm; //Mpc s^-0.4
      gal->Times[A] = run_globals.LTTime[snapshot] * time_unit; // s 
    }
    if (gal->count_SF > 0) {
      for (int i_SF = 0; i_SF < gal->count_SF; i_SF++)
        gal->Radii[i_SF] = gal->Prefactor[i_SF] * pow((gal->Times[i_SF] - run_globals.LTTime[snapshot] * time_unit), 0.4); 
    }
  }
  
  else {
    int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;
    mlog_error("So far, you can't relax the IRA");
  }
      
  double max = gal->Radii[0];    
    
  for (int i = 0; i < 70; i++) {       
     if(gal->Radii[i] > max)    
         max = gal->Radii[i];    
    }
    
  gal->RmetalBubble = max; 
  
  if (gal->RmetalBubble < 0.0)
    gal->RmetalBubble = 0.0;
}*/
void calc_metal_bubble(galaxy_t* gal, int snapshot) // result in internal units (Mpc/h) This new function assumes that a bubble will overtake a previous one in no more than 17 snapshots!
{
  int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;
  
  double UnitMass_in_g = run_globals.units.UnitMass_in_g;
  double UnitLength_in_cm = run_globals.units.UnitLength_in_cm;
  double energy_unit = run_globals.units.UnitEnergy_in_cgs;
  double time_unit = run_globals.units.UnitTime_in_s;
  
  double mm_stars = gal->NewStars[0]; //The last episode of SF
  double sn_energy = 0.0;
  
  // First evolve the existing RmetalBubble
  
  if (gal->RmetalBubble > 0.){
    gal->RmetalBubble = gal->PrefactorBubble * pow((gal->TimeBubble - run_globals.LTTime[snapshot] * time_unit), 0.4);
    mlog("Current Bubble = %f", MLOG_MESG, gal->RmetalBubble);
    if (gal->RmetalBubble > 10.)
      mlog("StrangeBubble = %f, Prefactor = %f, TimeBubble = %f", MLOG_MESG, gal->RmetalBubble, gal->PrefactorBubble, gal->TimeBubble);
    }
  
  // Now compute the last N_HISTORY_SNAPS bubble to see if any of those gets bigger than the existing one.
  
  double gas_density;  
  gas_density = (gal->HotGas + gal->ColdGas) * UnitMass_in_g / PROTONMASS / (4.0 * M_PI / 3.0 * pow(gal->Rvir * UnitLength_in_cm, 3.)); // cm^-3
  
  if (mm_stars > 1e-10) { 
    if (gal->Galaxy_Population == 3) //Crucial to update the galaxy index! 
      gal->Galaxy_Population = 2;
      }

  for (int i_burst = 0; i_burst < n_bursts; i_burst++) {
    double m_stars_II = gal->NewStars_II[i_burst];
    double m_stars_III = gal->NewStars_III[i_burst];
    double m_stars = m_stars_II + m_stars_III;
    
    // Compute the SN energy that drives the metal bubble. Should you include the ejection efficiency?  
    if (m_stars_II > 1e-10) {
      double metallicity = calc_metallicity(m_stars, gal->NewMetals[i_burst]);
      sn_energy += m_stars_II * get_SN_energy(0, metallicity) * energy_unit * calc_sn_ejection_eff(gal, snapshot, 2); // Are you sure of the unit? (You want this in erg)
      }
    else if (m_stars_III > 1e-10) {
      if (i_burst == 0) // You have both CC and PISN
        sn_energy += (get_SN_energy_PopIII(i_burst, snapshot, 0) + get_SN_energy_PopIII(i_burst, snapshot, 1))  * m_stars_III * calc_sn_ejection_eff(gal, snapshot, 3);
      else
        sn_energy += get_SN_energy_PopIII(i_burst, snapshot, 0) * m_stars_III * calc_sn_ejection_eff(gal, snapshot, 3);
      }
    if (i_burst != 0) {
      gal->Prefactor[n_bursts - i_burst] = gal->Prefactor[n_bursts - i_burst - 1];
      gal->Times[n_bursts - i_burst] = gal->Times[n_bursts - i_burst - 1];
      if (gal->Prefactor[n_bursts - i_burst] > 0.0)
        gal->Radii[n_bursts - i_burst] = gal->Prefactor[n_bursts - i_burst] * pow((gal->Times[n_bursts - i_burst] - run_globals.LTTime[snapshot] * time_unit), 0.4);
      else
        gal->Radii[n_bursts - i_burst] = 0.0;
      if (gal->Radii[n_bursts - i_burst] > gal->RmetalBubble) { //Look if one of the new bubbles is bigger than RmetalBubble
        //mlog("New Bubble = %f, Old Bubble = %f, SNenergy = %f", MLOG_MESG, gal->Radii[n_bursts - i_burst], gal->RmetalBubble, log10(sn_energy));
        gal->RmetalBubble = gal->Radii[n_bursts - i_burst];
        gal->PrefactorBubble = gal->Prefactor[n_bursts - i_burst];
        gal->TimeBubble = gal->Times[n_bursts - i_burst];
        }
      }
    }
  gal->Prefactor[0] = pow(sn_energy / (PROTONMASS * gas_density), 0.2) / UnitLength_in_cm; //Mpc s^-0.4
  if (gal->Prefactor[0] > 1e20)
    mlog("StrangePrefactor = %f, sn_energy = %f, gas_density = %f", MLOG_MESG, gal->Prefactor[0], log10(sn_energy), gas_density);
  gal->Times[0] = run_globals.LTTime[snapshot] * time_unit; // s 
  //gal->Radii[0] = gal->Prefactor[0] * pow((gal->Times[0] - run_globals.LTTime[snapshot] * time_unit), 0.4); //This is 0, so I could just put it as a 0.
  gal->Radii[0] = 0.0;
}
