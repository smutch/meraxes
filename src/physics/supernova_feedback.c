#include <assert.h>
#include <math.h>

#include "core/misc_tools.h"
#include "core/stellar_feedback.h"
#include "core/virial_properties.h"
#include "meraxes.h"
#include "supernova_feedback.h"

void update_reservoirs_from_sn_feedback(galaxy_t* gal,
                                        double m_reheat,
                                        double m_reheat_III,
                                        double m_reheat_II,
                                        double m_eject,
                                        double m_eject_III,
                                        double m_eject_II,
                                        double m_recycled,
                                        double m_recycled_III,
                                        double m_recycled_II,
                                        double new_metals)
{
  double metallicity;
  galaxy_t* central;
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);

  // If this is a ghost then it doesn't have an identified halo at this
  // snapshot.  We will therefore dump all of the reheated gas into the ghost's
  // hot halo, to be recollected and distributed when the ghost is reidentified
  // at a later time.
  if (gal->ghost_flag)
    central = gal;
  else
    central = gal->Halo->FOFGroup->FirstOccupiedHalo->Galaxy;

  gal->StellarMass -= m_recycled;
  gal->StellarMass_II -= m_recycled_II;
  gal->StellarMass_III -= m_recycled_III;
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
    m_reheat_III = gal->ColdGas - m_reheat_II;
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
      m_eject_III = central->HotGas - m_eject_II;
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

static inline double calc_sn_reheat_eff(galaxy_t* gal, int snapshot)
{
  double Vmax = gal->Vmax; // Vmax is in a unit of km/s
  double zplus1 = 1. + run_globals.ZZ[snapshot];
  physics_params_t* params = &run_globals.params.physics;
  int SnModel = params->SnModel;
  double SnReheatRedshiftDep = params->SnReheatRedshiftDep;
  double SnReheatEff = params->SnReheatEff;
  double SnReheatScaling = params->SnReheatScaling;
  double SnReheatNorm = params->SnReheatNorm;
  double SnReheatLimit = params->SnReheatLimit;
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

static inline double calc_sn_ejection_eff(galaxy_t* gal, int snapshot)
{
    double Vmax = gal->Vmax;    // Vmax is in a unit of km/s
    double zplus1 = 1. + run_globals.ZZ[snapshot];
    physics_params_t *params = &run_globals.params.physics;
    int SnModel = params->SnModel;
    double SnEjectionRedshiftDep = params->SnEjectionRedshiftDep;
    double SnEjectionEff = params->SnEjectionEff;
    double SnEjectionScaling = params->SnEjectionScaling;
    double SnEjectionNorm = params->SnEjectionNorm;
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

void delayed_supernova_feedback(galaxy_t* gal, int snapshot) // THIS FUNCTION IS STILL TO BE UPDATED!!!! 
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
  double fof_Vvir;
  // If we are at snapshot < N_HISTORY_SNAPS-1 then only try to look back to snapshot 0
  int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;

  // Loop through each of the last `N_HISTORY_SNAPS` recorded stellar mass
  // bursts and calculate the amount of energy and mass that they will release
  // in the current time step.
  for (int i_burst = 1; i_burst < n_bursts; i_burst++) {
    double m_stars = gal->NewStars[i_burst]; // IF YOU WANT TO DO THIS FOR DIFFERENT STELLAR POPULATIONS YOU NEED TO SAVE MORE OUTPUTS (NewStars_III, NewStars_II!!)
    double m_stars_II = gal->NewStars[i_burst];
    //double m_stars_III = gal->NewStars[i_burst];

    // Only need to do this if any stars formed in this history bin
    if (m_stars > 1e-10) {
      double metallicity = calc_metallicity(m_stars, gal->NewMetals[i_burst]);
      // Calculate recycled mass and metals by yield tables
      m_recycled += m_stars * get_recycling_fraction(i_burst, metallicity);
      m_recycled_II += m_stars_II * get_recycling_fraction(i_burst, metallicity);
      //m_recycled_III += m_stars_III * get_recycling_fraction(i_burst, metallicity);
      new_metals += m_stars * get_metal_yield(i_burst, metallicity);
      // Calculate SNII energy
      sn_energy += get_SN_energy(i_burst, metallicity) * m_stars;
      sn_energy_II += get_SN_energy(i_burst, metallicity) * m_stars_II;
      //sn_energy_III += get_SN_energy(i_burst, metallicity) * m_stars_III;
    }
  }

  m_reheat = calc_sn_reheat_eff(gal, snapshot) * sn_energy / get_total_SN_energy();
  sn_energy *= calc_sn_ejection_eff(gal, snapshot);
  m_reheat_II = calc_sn_reheat_eff(gal, snapshot) * sn_energy_II / get_total_SN_energy();
  sn_energy_II *= calc_sn_ejection_eff(gal, snapshot);
  //m_reheat_III = calc_sn_reheat_eff(gal, snapshot) * sn_energy_III / get_total_SN_energy();
  //sn_energy_III *= calc_sn_ejection_eff(gal, snapshot);
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
  //m_eject_III = calc_ejected_mass(&m_reheat_III, sn_energy_III, gal->Vvir, fof_Vvir);
  m_eject_II = calc_ejected_mass(&m_reheat_II, sn_energy_II, gal->Vvir, fof_Vvir);

  // Note that m_eject returned for ghosts by calc_ejected_mass() is
  // meaningless in the current physical prescriptions.  This fact is dealt
  // with in update_reservoirs_from_sn_feedback().

  assert(m_reheat >= 0);
  assert(m_eject >= 0);

  // update the baryonic reservoirs
  //update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, m_recycled, new_metals);
  update_reservoirs_from_sn_feedback(gal, m_reheat, m_reheat_III, m_reheat_II, m_eject, m_eject_III, m_eject_II, m_recycled, m_recycled_III, m_recycled_II, new_metals);
}

void contemporaneous_supernova_feedback(galaxy_t* gal,
                                        double* m_stars,
                                        double* m_stars_III,
                                        double* m_stars_II,
                                        int snapshot,
                                        double* m_reheat,
                                        double* m_reheat_III,
                                        double* m_reheat_II,
                                        double* m_eject,
                                        double* m_eject_III,
                                        double* m_eject_II,
                                        double* m_recycled,
                                        double* m_recycled_III,
                                        double* m_recycled_II,
                                        double* new_metals)
{
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
  double sn_energy = 0.0;
  double sn_energy_II = 0.0;
  double sn_energy_III = 0.0;
  // init (just in case!)
  *m_reheat = *m_recycled = *new_metals = *m_eject = 0.0;
  *m_reheat_III = *m_recycled_III = *m_eject_III = 0.0;
  *m_reheat_II = *m_recycled_II = *m_eject_II = 0.0;

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
    *m_recycled = *m_stars * get_recycling_fraction(0, metallicity);
    *m_recycled_II = *m_stars_II * get_recycling_fraction(0, metallicity);
    *m_recycled_III = *m_stars_III * get_recycling_fraction(0, metallicity);
    //*new_metals = *m_stars * get_metal_yield(0, metallicity);
    *new_metals = (*m_stars_II * get_metal_yield(0, metallicity) + *m_stars_III * get_metal_yield(0, metallicity)); // Not sure about this
  } else {
    // Recycling fraction and metals yield are input parameters when using IRA
    *m_recycled = *m_stars * run_globals.params.physics.SfRecycleFraction;
    *m_recycled_II = *m_stars_II * run_globals.params.physics.SfRecycleFraction; //Here you might want add different parameters
    *m_recycled_III = *m_stars_III * run_globals.params.physics.SfRecycleFraction;
    *new_metals = *m_stars * run_globals.params.physics.Yield;
    //*new_metals_II = *m_stars_II * run_globals.params.physics.Yield;
    //*new_metals_III = *m_stars_III * run_globals.params.physics.Yield;
  }
  // calculate the SNII energy and total reheated mass
  sn_energy = *m_stars * get_SN_energy(0, metallicity);
  sn_energy_II = *m_stars_II * get_SN_energy(0, metallicity);
  sn_energy_III = *m_stars_III * get_SN_energy(0, metallicity);
  *m_reheat = calc_sn_reheat_eff(gal, snapshot) * sn_energy / get_total_SN_energy();
  *m_reheat_II = calc_sn_reheat_eff(gal, snapshot) * sn_energy_II / get_total_SN_energy();
  *m_reheat_III = calc_sn_reheat_eff(gal, snapshot) * sn_energy_III / get_total_SN_energy();
  //*m_reheat = *m_reheat_II + *m_reheat_III;
  sn_energy *= calc_sn_ejection_eff(gal, snapshot);
  sn_energy_II *= calc_sn_ejection_eff(gal, snapshot);
  sn_energy_III *= calc_sn_ejection_eff(gal, snapshot);

  // We can only reheat as much gas as we have available.  Let's inforce this
  // now, to ensure that the maximal amount of available energy is used to
  // eject gas from the system.
  //if (*m_reheat > gal->ColdGas)
  //  *m_reheat = gal->ColdGas;
  if ((*m_reheat_III + *m_reheat_II) > gal->ColdGas){
    *m_reheat = gal->ColdGas;
    *m_reheat_III = gal->ColdGas - *m_reheat_II; // Why should you change III and not II ?
    }
    
  // attenuate the star formation if necessary, so that we are being consistent
  // if (*m_reheat + *m_stars - *m_recycled > gal->ColdGas)
  /*if (*m_reheat + *m_stars > gal->ColdGas) {
    // double frac = gal->ColdGas / (*m_reheat + *m_stars - *m_recycled);
    double frac = gal->ColdGas / (*m_reheat + *m_stars);
    *m_reheat *= frac;
    *m_stars *= frac;
    *m_recycled *= frac;
  }*/
  if ((*m_reheat_III + *m_reheat_II) + (*m_stars_III + *m_stars_II) > gal->ColdGas) {
    double frac_II = gal->ColdGas / (*m_reheat_II + *m_stars_II);
    double frac_III = gal->ColdGas / (*m_reheat_III + *m_stars_III);
    double frac = gal->ColdGas / *m_reheat;
    //double frac = frac_II + frac_III;
    *m_reheat_III *= frac_III;
    *m_reheat_II *= frac_II;
    *m_stars *= frac;
    *m_stars_III *= frac_III;
    *m_stars_II *= frac_II;
    *m_recycled *= frac;
    *m_recycled_III *= frac_III;
    *m_recycled_II *= frac_II;
  }
  //assert(*m_reheat_III >= 0);
  //assert(*m_recycled_III >= 0);
  assert(*new_metals >= 0);
  assert(*m_reheat_II >= 0);
  assert(*m_recycled_II >= 0);
  assert(*m_recycled >= 0);

  // how much mass is ejected due to this star formation episode? (ala Croton+ 2006)
  //*m_eject = calc_ejected_mass(m_reheat, sn_energy, gal->Vvir, gal->Halo->FOFGroup->Vvir);
  *m_eject = calc_ejected_mass(m_reheat, sn_energy, gal->Vvir, gal->Halo->FOFGroup->Vvir);
  *m_eject_III = calc_ejected_mass(m_reheat_III, sn_energy_III, gal->Vvir, gal->Halo->FOFGroup->Vvir);
  *m_eject_II = calc_ejected_mass(m_reheat_II, sn_energy_II, gal->Vvir, gal->Halo->FOFGroup->Vvir);

  assert(*m_reheat_III + *m_reheat_II >= 0);
  assert(*m_eject >= 0);
  assert(*m_reheat_III >= 0);
  assert(*m_eject_III >= 0);
  assert(*m_reheat_II >= 0);
  assert(*m_eject_II >= 0);
  
  mlog("Eject %f, EjectIII %f, Recyled %f, RecyledIII %f, reheat %f, reheatIII %f", MLOG_MESG, *m_eject, *m_eject_III, *m_recycled, *m_recycled_III, *m_reheat, * m_reheat_III);
}

void calc_metal_bubble(galaxy_t* gal, int snapshot) // Done! Result in internal units (Mpc/h)
{
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
  double mm_stars = gal->NewStars[0]; //The last episode of SF
  
  double UnitMass_in_g = run_globals.units.UnitMass_in_g;
  double UnitLength_in_cm = run_globals.units.UnitLength_in_cm;
  double time_unit = run_globals.units.UnitTime_in_s;
  
  int A = gal->count_SF;
  
  if (Flag_IRA == false) {
  
    if (mm_stars > 1e-10) { 

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
}
