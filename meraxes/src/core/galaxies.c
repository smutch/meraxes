#include "meraxes.h"
#include <math.h>

galaxy_struct* new_galaxy(int *unique_ID)
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
  gal->Cos_Inc           = 0.0;
  gal->MergTime          = 99999.9;

  for(int ii=0; ii<3; ii++)
  {
    gal->Pos[ii] = -99999.9;
    gal->Vel[ii] = -99999.9;
  }
  for(int ii=0; ii<NOUT; ii++)
    gal->Sfr[ii] = 0;

  gal->output_index = -1;
  gal->ghost_flag = false;

  return gal;
}

static double calculate_Mvir(run_globals_struct *run_globals, halo_struct *halo)
{
  if(halo->Type==0 && halo->Mvir)
    return halo->Mvir;
  else
    return (double)halo->Len * run_globals->params.PartMass;
}

static double calculate_Rvir(run_globals_struct *run_globals, double Mvir, int snapshot)
{

	double zplus1, hubble_of_z_sq, rhocrit, fac;
  double Hubble      = run_globals->Hubble;
  double Omega       = run_globals->params.Omega;
  double OmegaLambda = run_globals->params.OmegaLambda;
	
	zplus1 = 1 + run_globals->ZZ[snapshot];
	hubble_of_z_sq =
	  Hubble * Hubble *(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
	  OmegaLambda);
	
	rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals->G);
	fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
	
	return cbrt(Mvir * fac);
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
  gal->Rvir            = calculate_Rvir(run_globals, gal->Mvir, snapshot);
  gal->Vvir            = calculate_Vvir(run_globals, gal->Mvir, gal->Rvir);
  gal->Vmax            = halo->Vmax;
  gal->TreeFlags       = halo->TreeFlags;
  for (int ii=0; ii<3; ii++)
  {
    gal->Pos[ii]       = halo->Pos[ii];
    gal->Vel[ii]       = halo->Vel[ii];
  }
}
