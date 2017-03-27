#include <math.h>
#include "meraxes.h"
#include <assert.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

void update_reservoirs_from_sf(galaxy_t *gal, double new_stars)
{
  if (new_stars > 0)
  {
    double metallicity;
    double current_time;

    // update the galaxy's SFR value
    gal->Sfr += new_stars / gal->dt;
    assert(gal->Sfr >= 0);

    // update the stellar mass history
    gal->NewStars[0]       += new_stars;

    // instantaneous recycling approximation of stellar mass
    metallicity             = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

    gal->ColdGas           -= new_stars;
    gal->MetalsColdGas     -= new_stars * metallicity;
    gal->StellarMass       += new_stars;
    gal->GrossStellarMass  += new_stars;
    gal->FescWeightedGSM   += new_stars * run_globals.params.physics.ReionEscapeFrac;
    gal->MetalsStellarMass += new_stars * metallicity;

    // update the luminosities
    current_time            = run_globals.LTTime[gal->LastIdentSnap] - 0.5 * gal->dt;
    add_to_luminosities(gal, new_stars, metallicity, current_time);

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
}


void insitu_star_formation(galaxy_t *gal, int snapshot)
{
  // there is no point doing anything if there is no cold gas!
  if (gal->ColdGas > 1e-10)
  {
    double r_disk;
    double v_disk;
    double m_crit;
    double m_stars;

    double m_reheat;
    double m_eject;
    double m_recycled;
    double new_metals;

    double zplus1;
    double zplus1_n;

    zplus1   = 1.0 + run_globals.ZZ[snapshot];
    zplus1_n = pow(zplus1,run_globals.params.physics.SfEfficiencyScaling);

    double SfEfficiency     = run_globals.params.physics.SfEfficiency;
    double SfCriticalSDNorm = run_globals.params.physics.SfCriticalSDNorm;
    int    SfDiskVelOpt     = run_globals.params.physics.SfDiskVelOpt;
    int    SfPrescription   = run_globals.params.physics.SfPrescription;

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

    switch (SfPrescription)
    {
      case 1:
        // what is the critical mass within r_crit?
        // from Kauffmann (1996) eq7 x piR^2, (Vvir in km/s, reff in Mpc/h) in units of 10^10Msun/h
        m_crit = SfCriticalSDNorm * v_disk * r_disk;
        if (gal->ColdGas > m_crit)
          m_stars = zplus1_n * SfEfficiency * (gal->ColdGas - m_crit) / r_disk * v_disk * gal->dt;
        else
          // no star formation
          return;
        break;

      case 2:
        // f_h2 from Blitz & Rosolowski 2006 abd Bigiel+11 SF law
        m_stars = pressure_dependent_star_formation(gal, snapshot) * gal->dt;
        break;

      case 3:
        // GALFORM
        m_stars = gal->ColdGas / (r_disk / v_disk / 0.029 * pow(200. / v_disk,1.5)) * gal->dt;
        break;

      default:
        mlog_error("Unknown SfPrescription!");
        ABORT(EXIT_FAILURE);
        break;
    }

    if (m_stars > gal->ColdGas)
      m_stars = gal->ColdGas;

    // calculate the total supernova feedback which would occur if this star
    // formation happened continuously and evenly throughout the snapshot
    contemporaneous_supernova_feedback(gal, &m_stars, snapshot, &m_reheat, &m_eject, &m_recycled, &new_metals);

    // update the baryonic reservoirs (note that the order we do this in will change the result!)
    update_reservoirs_from_sf(gal, m_stars);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, m_recycled, new_metals);
  }
}


struct FR_parameters { double a; double b; double c; double d; };

static double integrand_p_dependent_SFR(double q, void *gal)
{
  struct FR_parameters * params       = (struct FR_parameters *)gal;

  double                 sigma_gas0   = (params->a);
  double                 sigma_stars0 = (params->b);
  double                 v_ratio      = (params->c);
  double                 reff         = (params->d);

  double                 G_SI         = GRAVITY * 1.e-3;
  double                 sf_effH      = 1.;

  double                 p_ext        = M_PI / 2.0 * G_SI * sigma_gas0 * exp(-q / reff) * (sigma_gas0 * exp(-q / reff)  + v_ratio * sqrt(sigma_stars0 * exp(-q / reff)) );
  double                 fmol         = 1.0 / (1.0 + pow(p_ext / 4.79e-13, -0.92));

  // double Surface0 = 200. * 1.989e30 / 3.086e16 / 3.086e16;  // 200M_sun/pc^(-2)
  // sf_effH=sf_effH*(1.+pow(sigma_gas0*exp(-q/reff)/Surface0,0.4));  //Lagos et al. 2011 paper

  double my_integrandR = q * fmol * sigma_gas0 * exp(-q / reff) * sf_effH;

  return my_integrandR;
}


static double p_dependent_SFR(double lower_limit, double upper_limit, double sigma_gas0, double sigma_stars0, double v_ratio, double reff)
{
  static gsl_function               FR;
  static gsl_integration_workspace *workspace;
  double                            result, abserr;
  size_t                            worksize   = 512;

  struct FR_parameters              parameters = {sigma_gas0, sigma_stars0, v_ratio, reff};

  workspace   = gsl_integration_workspace_alloc(worksize);

  FR.function = &integrand_p_dependent_SFR;
  FR.params   = &parameters;

  gsl_integration_qag(&FR, lower_limit, upper_limit, 1.0e-8, 1.0e-8, worksize, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
  gsl_integration_workspace_free(workspace);

  return result;
}


double pressure_dependent_star_formation(galaxy_t *gal, int snapshot)
{
  /*
   * Credit: Hansik Kim
   * Based on the SF prescription of Blitz & Rosolowski (2006).
   */

  double       SfEfficiency = run_globals.params.physics.SfEfficiency;
  double       Y_He         = run_globals.params.physics.Y_He;
  double       zplus1_n     = pow(1.0 + run_globals.ZZ[snapshot], run_globals.params.physics.SfEfficiencyScaling);
  run_units_t *units        = &(run_globals.units);
  double       G_SI         = GRAVITY * 1.e-3;

  // SF timescale:
  double sf_eff = 1.0 / 3.0e8 * SfEfficiency * zplus1_n;
  double MSFRR  = 0.0;

  if(gal->DiskScaleLength > 0.0)
  {
    double reff       = 1. * gal->DiskScaleLength;
    double sigma_gas0 = 0.76 * gal->ColdGas / (2.0 * M_PI * reff * reff);
    sigma_gas0   = sigma_gas0 * units->UnitMass_in_g / pow(units->UnitLength_in_cm,2);   // in g.cm^-2
    sigma_gas0   = sigma_gas0 * 1.0e-3 * 1.0e4;                                          // in kg.m^-2 //gas surface density

    double sigma_stars0 = gal->StellarMass / (2.0 * M_PI * reff * reff);                 // DRAGONS units
    sigma_stars0 = sigma_stars0 * units->UnitMass_in_g / pow(units->UnitLength_in_cm,2); // in g.cm^-2
    sigma_stars0 = sigma_stars0 * 1.0e-3 * 1.0e4;                                        // in kg.m^-2 // stellar surface density

    reff         = reff * units->UnitLength_in_cm / 100.0;                               // in m

    double vdisp_gas = 10.0e3;                                                           // m.s^-1 , following BR06
    double v_ratio   = vdisp_gas / sqrt(M_PI * G_SI * 0.14 * reff);                      //relation between stellar and gas dispersion

    if(sigma_gas0 > 0.0)
    {
      double p_ext = M_PI / 2.0 * G_SI * sigma_gas0 * (sigma_gas0  + v_ratio * sqrt(sigma_stars0) );
      double MSFR  = 1.0 / (1.0 + pow(p_ext / 4.79e-13, -0.92));
      gal->H2Frac = MSFR; // Molecular hydrogen fraction, f_(H2Mass)

      if ((MSFR < 0.0) || (MSFR > 1.0))
      {
        mlog_error("BR06 - H2Frac = %e\n", MSFR);
        ABORT(66);
      }

      // Bigiel+11 SF law
      // TODO: PUT THIS BACK!
      MSFRR       = p_dependent_SFR(0, 5 * reff, sigma_gas0, sigma_stars0, v_ratio,reff);
      gal->H2Mass = 2. * M_PI * MSFRR * 1.0e3 / units->UnitMass_in_g; // Molecular hydrogen mass
      if (gal->H2Mass > (1.-Y_He) * gal->ColdGas)
        gal->H2Mass = (1.-Y_He) * gal->ColdGas;
      gal->HIMass = (1.-Y_He) * gal->ColdGas - gal->H2Mass;                //hydrogen mass
      MSFRR       = MSFRR * 2.0 * M_PI * sf_eff / SEC_PER_YEAR;
      MSFRR       = MSFRR * 1.0e3; // in g/s
    }
    else
    {
      MSFRR       = 0.0;
      gal->H2Frac = 0.0;
      gal->H2Mass = 0.0;
      gal->HIMass = 0.0;
    }
  }
  else
  {
    MSFRR       = 0.0;
    gal->H2Frac = 0.0;
    gal->H2Mass = 0.0;
    gal->HIMass = 0.0;
  }

  MSFRR = MSFRR / units->UnitMass_in_g * units->UnitTime_in_s; // SFR in DRAGONS units

  return MSFRR;
}
