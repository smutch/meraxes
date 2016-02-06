#include "meraxes.h"
#include <math.h>
#include <assert.h>

static double inline M0(double z)
{
  return Tvir_to_Mvir(run_globals.params.physics.ReionSobacchi_T0, z);
}


static double inline Mcool(double z)
{
  return Tvir_to_Mvir(run_globals.params.physics.ReionTcool, z);
}


static double sobacchi_Mvir_min(double z)
{
  // Calculate the minimum halo mass capable of hosting star forming galaxies
  // following Sobacchi & Mesinger 2013b.

  double current_Mcool = Mcool(z);
  double current_M0    = M0(z);
  double Mvir_min;
  double g_term;
  physics_params_t *params = &(run_globals.params.physics);

  g_term = 1. / (1. + exp((z - (params->ReionSobacchi_Zre - params->ReionSobacchi_DeltaZsc)) / params->ReionSobacchi_DeltaZre));
  Mvir_min = current_Mcool * pow(current_M0 / current_Mcool, g_term);

  return Mvir_min;
}


double sobacchi2013_modifier(double Mvir, double redshift)
{

  // This asserion is a check for validity of using Mvir_min here.  Really
  // Mvir_crit should be used.  Mvir_min = Mvir_crit as long as the halo mass
  // is a bit larger than Mcool.
  // double current_Mcool = Mcool(redshift);
  // assert((redshift < 6.0) || (Mvir > 1.05*current_Mcool));

  double Mvir_min = sobacchi_Mvir_min(redshift);

  return pow(2.0, -Mvir_min / Mvir);
}


static double precomputed_Mcrit_modifier(galaxy_t *gal, double Mvir, int snapshot)
{
  // Use precomputed Mcrit (Mfilt) mass values as a function of redshif in
  // order to calculate the baryon fraction modifier using the Sobacchi
  // formalism.

  double Mvir_crit = run_globals.params.MvirCrit[snapshot];
  if (gal != NULL)
    gal->MvirCrit = Mvir_crit;
  return pow(2.0, -Mvir_crit / Mvir);
}


double gnedin2000_modifer(double Mvir, double redshift)
{
  // NOTE THAT PART OF THIS CODE IS COPIED VERBATIM FROM THE CROTON ET AL. 2006 SEMI-ANALYTIC MODEL.
  // WITH A COUPLE OF BUGFIXES SO THAT EQUATIONS MATCH KRAVTSOV+ 2004

  double a0;
  double ar;
  double alpha;
  double a;
  double f_of_a;
  double a_on_a0;
  double a_on_ar;
  double Mfiltering;
  double Mjeans;
  double Mchar;
  double mass_to_use;
  double modifier;

  a0 = 1.0 / (1.0 + run_globals.params.physics.ReionGnedin_z0);
  ar = 1.0 / (1.0 + run_globals.params.physics.ReionGnedin_zr);

  // we employ the reionization recipie described in Gnedin (2000), however use the fitting
  // formulas given by Kravtsov et al (2004) Appendix B

  // alpha=6 gives the best fit to the Gnedin data
  alpha = 6.0;

  // calculate the filtering mass
  a       = 1.0 / (1.0 + redshift);
  a_on_a0 = a / a0;
  a_on_ar = a / ar;

  if (a <= a0)
    f_of_a = 3.0 * a / ((2.0 + alpha) * (5.0 + 2.0 * alpha)) * pow(a_on_a0, alpha);
  else if ((a > a0) && (a < ar))
    f_of_a =
      (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
      a * a / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5)));
  else
    f_of_a =
      (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) /
                              (5.0 + 2.0 * alpha)) + (ar * ar / 10.0) * (5.0 - 4.0 * pow(a_on_ar,
                                                                                         -0.5)) - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5)) + a
                   * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * pow(a_on_ar, -0.5)));

  // this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21
  Mjeans     = 25.0 * pow(run_globals.params.OmegaM, -0.5) * 2.21;
  Mfiltering = Mjeans * pow(f_of_a, 1.5);

  // calculate the characteristic atomic cooling mass coresponding to a halo temperature of 10^4K
  Mchar = Mcool(redshift);

  // we use the maximum of Mfiltering and Mchar
  mass_to_use = (Mfiltering > Mchar) ? Mfiltering : Mchar;
  modifier    = 1.0 / pow(1.0 + 0.26 * (mass_to_use / Mvir), 3.0);

  return modifier;
}


double reionization_modifier(galaxy_t *gal, double Mvir, int snapshot)
{
  double redshift;
  double modifier;

  redshift = run_globals.ZZ[snapshot];

  if ((tocf_params.uvb_feedback) && (run_globals.params.TOCF_Flag))
  {
    modifier = tocf_modifier(gal, Mvir);
    return modifier;
  }

  switch (run_globals.params.physics.Flag_ReionizationModifier)
  {
  case 1:
    // Sobacchi & Mesinger 2013 global reionization scheme
    modifier = sobacchi2013_modifier(Mvir, redshift);
    break;

  case 2:
    // Gnedin 2000 global reionization modifier
    modifier = gnedin2000_modifer(Mvir, redshift);
    break;

  case 3:
    // Precomputed mean Mcrit values
    modifier = precomputed_Mcrit_modifier(gal, Mvir, snapshot);
    break;

  default:
    modifier = 1.0;
    break;
  }

  return modifier;
}


