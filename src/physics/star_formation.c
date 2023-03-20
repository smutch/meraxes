#include <assert.h>
#include <gsl/gsl_integration.h>

#include "core/magnitudes.h"
#include "core/misc_tools.h"
#include "core/reionization.h"
#include "meraxes.h"
#include "star_formation.h"
#include "supernova_feedback.h"

static void backfill_ghost_star_formation(galaxy_t* gal, double m_stars, double m_stars_III, double m_stars_II, double sfr, double metallicity, int snapshot)
{
  double* LTTime = run_globals.LTTime;
  double burst_time = LTTime[gal->LastIdentSnap] - gal->dt * 0.5;

  for (int snap = snapshot - 1, ii = 1; snap >= gal->LastIdentSnap; --snap, ++ii) {

    if (LTTime[snap] > burst_time) {
#ifdef CALC_MAGS
      if (sfr > 0.)
        add_luminosities(&run_globals.mag_params, gal, snap, metallicity, sfr);
#endif
      if (ii < N_HISTORY_SNAPS) {
        gal->NewStars[ii] += m_stars;
        gal->NewStars_II[ii] += m_stars_II;
        gal->NewStars_III[ii] += m_stars_III;
        gal->NewMetals[0] += m_stars * metallicity;
      }
      update_galaxy_fesc_vals(gal, m_stars, snap);
      break;
    }
  }
}

void update_reservoirs_from_sf(galaxy_t* gal, double new_stars, int snapshot, SFtype type) //If you want to save Sfr PopIII / PopII change this
{
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);
  if (new_stars > 0) {
    double metallicity;
    bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

    // instantaneous recycling approximation of stellar mass
    metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

    // update the galaxy's SFR value
    double sfr = new_stars / gal->dt;
    gal->Sfr += sfr;
    assert(gal->Sfr >= 0);

    gal->ColdGas -= new_stars;
    gal->MetalsColdGas -= new_stars * metallicity;
    gal->StellarMass += new_stars;
    gal->GrossStellarMass += new_stars;
    gal->MetalsStellarMass += new_stars * metallicity;
    
    if (Flag_Metals == true) {
      if (gal->Galaxy_Population == 2)
        gal->StellarMass_II += new_stars;
      else
        gal->StellarMass_III += new_stars;
    }

    if ((type == INSITU) && !Flag_IRA && (gal->LastIdentSnap < (snapshot - 1))) {
      // If this is a reidentified ghost, then back fill NewStars and
      // escape fraction dependent properties to reflect this new insitu
      // SF burst.
      if (gal->Galaxy_Population == 2)
        backfill_ghost_star_formation(gal, new_stars, 0, new_stars, sfr, metallicity, snapshot);
      else
        backfill_ghost_star_formation(gal, new_stars, new_stars, 0, sfr, metallicity, snapshot);
      } else {
      // update the stellar mass history assuming the burst is happening in this snapshot
#ifdef CALC_MAGS
      if (sfr > 0.)
        add_luminosities(&run_globals.mag_params, gal, snapshot, metallicity, sfr);
#endif
      gal->NewStars[0] += new_stars;
      if (gal->Galaxy_Population == 2)
        gal->NewStars_II[0] += new_stars;
      else
        gal->NewStars_III[0] += new_stars;
      gal->NewMetals[0] += new_stars * metallicity;

      update_galaxy_fesc_vals(gal, new_stars, snapshot);
      
    }

    // Check the validity of the modified reservoir values.
    // Note that the ColdGas reservers *can* be negative at this point.  This
    // is because some fraction of the stars in this burst will go nova and
    // return mass to the ISM.  This will be accounted for when we update the
    // reservoirs due to supernova feedback.
    if (gal->StellarMass < 0)
      gal->StellarMass = 0.0;
    if (gal->MetalsStellarMass < 0)
      gal->MetalsStellarMass = 0.0;
  }
  // The same halo cannot form Pop.III twice, so I read if it experienced previously a SF episode.
  if (Flag_Metals == true) {
    if ((gal->NewStars[1] > 1e-10) && (gal->Galaxy_Population == 3))
      gal->Galaxy_Population = 2;
  }
}

void insitu_star_formation(galaxy_t* gal, int snapshot, int flag_population) // Added Flag Population
{
  // there is no point doing anything if there is no cold gas!
  if (gal->ColdGas > 1e-10) {
    double r_disk;
    double v_disk;
    double m_crit;
    double m_stars; // Do you want to save sum? Probably not here
    double m_stars_II; // Pop.II
    double m_stars_III; // Pop.III

    double m_reheat;
    double m_reheat_II;
    double m_reheat_III;
    double m_eject;
    double m_eject_II;
    double m_eject_III;
    double m_recycled;
    double m_recycled_II;
    double m_recycled_III;
    double new_metals;

    double zplus1;
    double zplus1_n;

    zplus1 = 1.0 + run_globals.ZZ[snapshot];
    zplus1_n = pow(zplus1, run_globals.params.physics.SfEfficiencyScaling); //Later added separate parameters for Pop III and Pop II

    double SfEfficiency = run_globals.params.physics.SfEfficiency; //Later added separate parameters for Pop III and Pop II
    double SfCriticalSDNorm = run_globals.params.physics.SfCriticalSDNorm;
    int SfDiskVelOpt = run_globals.params.physics.SfDiskVelOpt;
    int SfPrescription = run_globals.params.physics.SfPrescription;

    // What velocity are we going to use as a proxy for the disk rotation velocity?
    switch (SfDiskVelOpt) {
      case 1:
        v_disk = gal->Vmax;
        break;
      case 2:
        v_disk = gal->Vvir;
        break;
      default:
        mlog_error("Unrecognised value for SfVelocityOpt parameter! Defaulting to v_disk=Vmax.");
        v_disk = gal->Vmax;
        break;
    }

    // calculate disk scalelength using Mo, Mao & White (1998) eqn. 12 and
    // multiply it by 3 to approximate the star forming region size (ala
    // Croton+ 2006).
    r_disk = gal->DiskScaleLength * 3.0;

    switch (SfPrescription) {
      case 1:
        // what is the critical mass within r_crit?
        // from Kauffmann (1996) eq7 x piR^2, (Vvir in km/s, reff in Mpc/h) in units of 10^10Msun/h
        m_crit = SfCriticalSDNorm * v_disk * r_disk;
        if (gal->ColdGas > m_crit){
          if (flag_population == 2){
            m_stars_II = zplus1_n * SfEfficiency * (gal->ColdGas - m_crit) / r_disk * v_disk * gal->dt;
            m_stars = m_stars_II;
            }
          else if (flag_population == 3){
            m_stars_III = zplus1_n * SfEfficiency * (gal->ColdGas - m_crit) / r_disk * v_disk * gal->dt;
            m_stars = m_stars_III;
            }
          }
        else
          // no star formation
          return;
        break;

      case 2:
        // f_h2 from Blitz & Rosolowski 2006 abd Bigiel+11 SF law
        if (flag_population == 2){
          m_stars_II = pressure_dependent_star_formation(gal, snapshot) * gal->dt;
          m_stars = m_stars_II;
          }
        else{
          m_stars_III = pressure_dependent_star_formation(gal, snapshot) * gal->dt;
          m_stars = m_stars_III;
          }
        break;

      case 3:
        // GALFORM
        if (flag_population == 2){
          m_stars_II = gal->ColdGas / (r_disk / v_disk / 0.029 * pow(200. / v_disk, 1.5)) * gal->dt;
          m_stars = m_stars_II;
          }
        else{
          m_stars_III = gal->ColdGas / (r_disk / v_disk / 0.029 * pow(200. / v_disk, 1.5)) * gal->dt;
          m_stars = m_stars_III;
          }
        break;

      default:
        m_stars = m_stars_III = m_stars_II = 0;
        mlog_error("Unknown SfPrescription!");
        ABORT(EXIT_FAILURE);
        break;
    }
    //m_stars = m_stars_II + m_stars_III;
    if (m_stars > gal->ColdGas)
      m_stars = gal->ColdGas;
      //m_stars_III = gal->ColdGas - m_stars_II;
    if (m_stars_III > m_stars)
      m_stars_III = m_stars;
    if (m_stars_II > m_stars)
      m_stars_II = m_stars;

    // calculate the total supernova feedback which would occur if this star
    // formation happened continuously and evenly throughout the snapshot
    contemporaneous_supernova_feedback(gal, &m_stars, &m_stars_III, &m_stars_II, snapshot, &m_reheat, &m_reheat_III, &m_reheat_II, 
                                       &m_eject, &m_eject_III, &m_eject_II, &m_recycled, &m_recycled_III, &m_recycled_II, &new_metals);
    // update the baryonic reservoirs (note that the order we do this in will change the result!)
    update_reservoirs_from_sf(gal, m_stars, snapshot, INSITU);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_reheat_III, m_reheat_II, m_eject, m_eject_III, m_eject_II, m_recycled, m_recycled_III, m_recycled_II, new_metals);
    //} ? Sure
  }
}

static double integrand_p_dependent_SFR(double q, void* gal)
{
  struct FR_parameters* params = (struct FR_parameters*)gal;

  double sigma_gas0 = (params->a);
  double sigma_stars0 = (params->b);
  double v_ratio = (params->c);
  double reff = (params->d);

  double G_SI = GRAVITY * 1.e-3;
  double sf_effH = 1.;

  double p_ext = M_PI / 2.0 * G_SI * sigma_gas0 * exp(-q / reff) *
                 (sigma_gas0 * exp(-q / reff) + v_ratio * sqrt(sigma_stars0 * exp(-q / reff)));
  double fmol = 1.0 / (1.0 + pow(p_ext / 4.79e-13, -0.92));

  // double Surface0 = 200. * 1.989e30 / 3.086e16 / 3.086e16;  // 200M_sun/pc^(-2)
  // sf_effH=sf_effH*(1.+pow(sigma_gas0*exp(-q/reff)/Surface0,0.4));  //Lagos et al. 2011 paper

  double my_integrandR = q * fmol * sigma_gas0 * exp(-q / reff) * sf_effH;

  return my_integrandR;
}

static double p_dependent_SFR(double lower_limit,
                              double upper_limit,
                              double sigma_gas0,
                              double sigma_stars0,
                              double v_ratio,
                              double reff)
{
  static gsl_function FR;
  static gsl_integration_workspace* workspace;
  double result, abserr;
  size_t worksize = 512;

  struct FR_parameters parameters = { sigma_gas0, sigma_stars0, v_ratio, reff };

  workspace = gsl_integration_workspace_alloc(worksize);

  FR.function = &integrand_p_dependent_SFR;
  FR.params = &parameters;

  gsl_integration_qag(
    &FR, lower_limit, upper_limit, 1.0e-8, 1.0e-8, worksize, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
  gsl_integration_workspace_free(workspace);

  return result;
}

double pressure_dependent_star_formation(galaxy_t* gal, int snapshot)
{
  /*
   * Credit: Hansik Kim
   * Based on the SF prescription of Blitz & Rosolowski (2006).
   */

  double SfEfficiency = run_globals.params.physics.SfEfficiency;
  double Y_He = run_globals.params.physics.Y_He;
  double zplus1_n = pow(1.0 + run_globals.ZZ[snapshot], run_globals.params.physics.SfEfficiencyScaling);
  run_units_t* units = &(run_globals.units);
  double G_SI = GRAVITY * 1.e-3;

  // SF timescale:
  double sf_eff = 1.0 / 3.0e8 * SfEfficiency * zplus1_n;
  double MSFRR = 0.0;

  if (gal->DiskScaleLength > 0.0) {
    double reff = 1. * gal->DiskScaleLength;
    double sigma_gas0 = 0.76 * gal->ColdGas / (2.0 * M_PI * reff * reff);
    sigma_gas0 = sigma_gas0 * units->UnitMass_in_g / pow(units->UnitLength_in_cm, 2); // in g.cm^-2
    sigma_gas0 = sigma_gas0 * 1.0e-3 * 1.0e4; // in kg.m^-2 //gas surface density

    double sigma_stars0 = gal->StellarMass / (2.0 * M_PI * reff * reff);                  // DRAGONS units
    sigma_stars0 = sigma_stars0 * units->UnitMass_in_g / pow(units->UnitLength_in_cm, 2); // in g.cm^-2
    sigma_stars0 = sigma_stars0 * 1.0e-3 * 1.0e4; // in kg.m^-2 // stellar surface density

    reff = reff * units->UnitLength_in_cm / 100.0; // in m

    double vdisp_gas = 10.0e3;                                    // m.s^-1 , following BR06
    double v_ratio = vdisp_gas / sqrt(M_PI * G_SI * 0.14 * reff); // relation between stellar and gas dispersion

    if (sigma_gas0 > 0.0) {
      double p_ext = M_PI / 2.0 * G_SI * sigma_gas0 * (sigma_gas0 + v_ratio * sqrt(sigma_stars0));
      double MSFR = 1.0 / (1.0 + pow(p_ext / 4.79e-13, -0.92));
      gal->H2Frac = MSFR; // Molecular hydrogen fraction, f_(H2Mass)

      if ((MSFR < 0.0) || (MSFR > 1.0)) {
        mlog_error("BR06 - H2Frac = %e\n", MSFR);
        ABORT(66);
      }

      // Bigiel+11 SF law
      // TODO: PUT THIS BACK!
      MSFRR = p_dependent_SFR(0, 5 * reff, sigma_gas0, sigma_stars0, v_ratio, reff);
      gal->H2Mass = 2. * M_PI * MSFRR * 1.0e3 / units->UnitMass_in_g; // Molecular hydrogen mass
      if (gal->H2Mass > (1. - Y_He) * gal->ColdGas)
        gal->H2Mass = (1. - Y_He) * gal->ColdGas;
      gal->HIMass = (1. - Y_He) * gal->ColdGas - gal->H2Mass; // hydrogen mass
      MSFRR = MSFRR * 2.0 * M_PI * sf_eff / SEC_PER_YEAR;
      MSFRR = MSFRR * 1.0e3; // in g/s
    } else {
      MSFRR = 0.0;
      gal->H2Frac = 0.0;
      gal->H2Mass = 0.0;
      gal->HIMass = 0.0;
    }
  } else {
    MSFRR = 0.0;
    gal->H2Frac = 0.0;
    gal->H2Mass = 0.0;
    gal->HIMass = 0.0;
  }

  MSFRR = MSFRR / units->UnitMass_in_g * units->UnitTime_in_s; // SFR in DRAGONS units

  return MSFRR;
}
