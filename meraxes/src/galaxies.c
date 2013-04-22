#include "meraxes.h"
#include <math.h>

void init_galaxy(galaxy_struct *gal)
{

  // Initialise the properties
  gal->Type              = -1;
  gal->Halo              = NULL;
  gal->NextGalInHalo     = NULL;
  gal->Next              = NULL;
  gal->MergerTarget      = NULL;
  gal->HaloDesc          = -1;
  gal->HaloNGal          = 0;
  gal->CentralGal        = -1;
  gal->Len               = 0;
  gal->CentralMvir       = 0.0;
  gal->Mvir              = 0.0;
  gal->dM                = 0.0;
  gal->dMdt              = 0.0;
  gal->Rvir              = 0.0;
  gal->Vvir              = 0.0;
  gal->Vmax              = 0.0;
  gal->StellarMass       = 0.0;
  gal->BlackHoleMass     = 0.0;
  gal->Cos_Inc           = 0.0;
  gal->MergTime          = 99999.9;

  for(int ii=0; ii<3; ii++)
  {
    gal->Pos[ii] = -99999.9;
    gal->Vel[ii] = -99999.9;
  }
  for(int ii=0; ii<NOUT; ii++)
    gal->Sfr[ii] = -99999.9;

}

static double calculate_Vvir(run_globals_struct *run_globals, halo_struct *halo)
{
  return sqrt(run_globals->G * halo->Mvir/1.0e10 / halo->Rvir);
}

void copy_halo_to_galaxy(run_globals_struct *run_globals, halo_struct *halo, galaxy_struct *gal)
{
  gal->Halo            = halo;
  gal->Type            = halo->Type;
  gal->Len             = halo->Len;
  gal->HaloDescIndex   = halo->DescIndex;
  gal->Mvir            = halo->Mvir/1.0e10;
  gal->Rvir            = halo->Rvir;
  gal->Vvir            = calculate_Vvir(run_globals, halo);
  gal->Vmax            = halo->Vmax;
  gal->TreeFlags       = halo->TreeFlags;
}
