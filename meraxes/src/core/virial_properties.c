#include "meraxes.h"
#include <math.h>

static inline double E_z(double z, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double one_plus_z;
  double one_plus_z_sq;
  double one_plus_z_cu;
  double result;

  one_plus_z    = 1.+z;
  one_plus_z_sq = one_plus_z*one_plus_z;
  one_plus_z_cu = one_plus_z_sq*one_plus_z;
  result        = sqrt(OmegaM*one_plus_z_cu+OmegaK*one_plus_z_sq+OmegaLambda);

  return result;
}

inline double Omega_z(double redshift, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double Ez;
  double one_plus_z_cube;

  Ez              = E_z(redshift, OmegaM, OmegaK, OmegaLambda);
  one_plus_z_cube = (1.+redshift)*(1.+redshift)*(1.+redshift);

  return OmegaM*one_plus_z_cube/(Ez*Ez);
}

inline double Delta_vir(double redshift, run_globals_t *run_globals)
{
  // Function stolen and adapted from gbpCosmo
  double x;
  double Omega;
  double OmegaM      = run_globals->params.OmegaM;
  double OmegaK      = run_globals->params.OmegaK;
  double OmegaLambda = run_globals->params.OmegaLambda;

  Omega = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);
  x     = Omega-1.;

  return (18.*M_PI*M_PI+82*x-39*x*x)/Omega;
}

double calculate_Mvir(run_globals_t *run_globals, halo_t *halo)
{
  if(halo->Type==0 && halo->Mvir)
    return halo->Mvir;
  else
    return (double)halo->Len * run_globals->params.PartMass;
}

float calculate_Rvir(run_globals_t *run_globals, halo_t *halo, double Mvir, int snapshot)
{

  double zplus1;
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble      = run_globals->Hubble;
  double OmegaM      = run_globals->params.OmegaM;
  double OmegaK      = run_globals->params.OmegaK;
  double OmegaLambda = run_globals->params.OmegaLambda;

  zplus1 = 1 + run_globals->ZZ[snapshot];
  hubble_of_z_sq = Hubble * Hubble *(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals->G);

  Delta = Delta_vir(run_globals->ZZ[snapshot], run_globals);

  fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);

  return cbrt((float)Mvir * fac);

}

float calculate_Vvir(run_globals_t *run_globals, double Mvir, float Rvir)
{
  return sqrt((float)(run_globals->G) * (float)Mvir / Rvir);
}
