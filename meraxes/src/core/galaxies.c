#include "meraxes.h"
#include <math.h>

galaxy_struct* new_galaxy(run_globals_struct *run_globals, int *unique_ID)
{
 
  galaxy_struct *gal = NULL;

  gal = SID_malloc(sizeof(galaxy_struct));

  // Initialise the properties
  gal->id_MBP            = 0;
  gal->ID                = (*unique_ID)++;
  gal->Type              = -1;
  gal->OldType           = -1;
  gal->SnapSkipCounter   = 0;
  gal->HaloDescIndex     = -1;
  gal->TreeFlags         = -1;
  gal->Halo              = NULL;
  gal->FirstGalInHalo    = NULL;
  gal->NextGalInHalo     = NULL;
  gal->Next              = NULL;
  gal->MergerTarget      = NULL;
  gal->Len               = 0;
  gal->dt                = 0.0;
  gal->LTTime            = 0.0;
  gal->Mvir              = 0.0;
  gal->dM                = 0.0;
  gal->Rvir              = 0.0;
  gal->Vvir              = 0.0;
  gal->Vmax              = 0.0;
  gal->StellarMass       = 0.0;
  gal->Cos_Inc           = gsl_rng_uniform(run_globals->random_generator);
  gal->MergTime          = 99999.9;
  gal->CellIonization    = 0.0;

  for(int ii=0; ii<3; ii++)
  {
    gal->Pos[ii] = -99999.9;
    gal->Vel[ii] = -99999.9;
  }
  for(int ii=0; ii<NOUT; ii++)
    gal->Sfr[ii] = 0;

  gal->output_index = -1;
  gal->ghost_flag = false;

  init_luminosities(run_globals, gal);
  return gal;
}

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

static inline double Omega_z(double redshift, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double Ez;
  double one_plus_z_cube;

  Ez              = E_z(redshift, OmegaM, OmegaK, OmegaLambda);
  one_plus_z_cube = (1.+redshift)*(1.+redshift)*(1.+redshift);

  return OmegaM*one_plus_z_cube/(Ez*Ez);
}

static inline double Delta_vir(double redshift, run_globals_struct *run_globals)
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

static double calculate_Mvir(run_globals_struct *run_globals, halo_struct *halo)
{
  if(halo->Type==0 && halo->Mvir)
    return halo->Mvir;
  else
    return (double)halo->Len * run_globals->params.PartMass;
}

static double calculate_Rvir(run_globals_struct *run_globals, halo_struct *halo, double Mvir, int snapshot)
{

  if(halo->Type==0 && halo->Rvir)
    return halo->Rvir;
  else
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
    // Delta = 200.;

    fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);

    return cbrt(Mvir * fac);
  }

}

static double calculate_Vvir(run_globals_struct *run_globals, double Mvir, double Rvir)
{
  return sqrt(run_globals->G * Mvir / Rvir);
}

void copy_halo_to_galaxy(run_globals_struct *run_globals, halo_struct *halo, galaxy_struct *gal, int snapshot)
{
  gal->id_MBP          = halo->id_MBP;
  gal->Type            = halo->Type;
  gal->Len             = halo->Len;
  gal->SnapSkipCounter = halo->SnapOffset;
  gal->HaloDescIndex   = halo->DescIndex;
  gal->Mvir            = calculate_Mvir(run_globals, halo);
  gal->Rvir            = calculate_Rvir(run_globals, halo, gal->Mvir, snapshot);
  gal->Vvir            = calculate_Vvir(run_globals, gal->Mvir, gal->Rvir);
  gal->Vmax            = halo->Vmax;
  gal->TreeFlags       = halo->TreeFlags;
  for (int ii=0; ii<3; ii++)
  {
    gal->Pos[ii]       = halo->Pos[ii];
    gal->Vel[ii]       = halo->Vel[ii];
  }
}
