#include "meraxes.h"
#include <math.h>

static inline double E_z(double z, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double one_plus_z;
  double one_plus_z_sq;
  double one_plus_z_cu;
  double result;

  one_plus_z    = 1. + z;
  one_plus_z_sq = one_plus_z * one_plus_z;
  one_plus_z_cu = one_plus_z_sq * one_plus_z;
  result        = sqrt(OmegaM * one_plus_z_cu + OmegaK * one_plus_z_sq + OmegaLambda);

  return result;
}

static inline double Omega_z(double redshift, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double Ez;
  double one_plus_z_cube;

  Ez              = E_z(redshift, OmegaM, OmegaK, OmegaLambda);
  one_plus_z_cube = (1. + redshift) * (1. + redshift) * (1. + redshift);

  return OmegaM * one_plus_z_cube / (Ez * Ez);
}

static inline double Delta_vir(double redshift)
{
  // Function stolen and adapted from gbpCosmo
  double x;
  double Omega;
  double OmegaM      = run_globals.params.OmegaM;
  double OmegaK      = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;

  Omega = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);
  x     = Omega - 1.;

  return (18. * M_PI * M_PI + 82 * x - 39 * x * x) / Omega;
}


//! Calculates Mvir in internal units (1.e10 h^{-1}Msol), given Tvir (in K) and a redshift (z)
double Tvir_to_Mvir(double T, double z)
{
  double OmegaM      = run_globals.params.OmegaM;
  double OmegaK      = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;

  double mu;         //!< Mean molecular weight (ionized gas)
  if (T < 9.99999e3) // Neutral IGM
    mu = 1.22;
  else               // Ionised IGM
    mu = 0.59;

  double z_term     = pow((1. + z) / 10., -1.5);
  double T_term     = pow(T / 1.98e4, 1.5);
  double cosmo_term = pow(OmegaM / Omega_z(z, OmegaM, OmegaK, OmegaLambda) *
                          Delta_vir(z) / 18. / (M_PI * M_PI), -0.5);
  double mol_term = pow(mu / 0.6, -1.5);

  return 0.01 * mol_term * cosmo_term * T_term * z_term;
}


double calculate_Mvir(double Mvir, int len)
{
  if ((len < 0) && (Mvir > 0))
    return Mvir;
  else
    return (double)len * run_globals.params.PartMass;
}


double hubble_at_snapshot(int snapshot)
{
  double Hubble      = run_globals.Hubble;
  double OmegaM      = run_globals.params.OmegaM;
  double OmegaK      = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1      = run_globals.ZZ[snapshot] + 1;

  return Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda);
}


double hubble_time(int snapshot)
{
  return 1.0 / hubble_at_snapshot(snapshot);
}


double calculate_Rvir(double Mvir, int snapshot)
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;

  hubble_of_z_sq = pow(hubble_at_snapshot(snapshot), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);

  Delta = Delta_vir(run_globals.ZZ[snapshot]);

  fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);

  return cbrt(Mvir * fac);
}

double calculate_Vvir(double Mvir, double Rvir)
{
  return sqrt((run_globals.G) * Mvir / Rvir);
}

double calculate_spin_param(halo_t *halo)
{
  float spin;

  spin = sqrt(halo->AngMom[0] * halo->AngMom[0] +
              halo->AngMom[1] * halo->AngMom[1] +
              halo->AngMom[2] * halo->AngMom[2]);

  spin = spin / (1.414213562 * halo->Vvir * halo->Rvir);

  return (double)spin;
}
