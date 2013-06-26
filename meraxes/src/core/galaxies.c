#include "meraxes.h"
#include <math.h>

galaxy_struct* new_galaxy(int *unique_ID)
{
 
  galaxy_struct *gal = NULL;

  gal = SID_malloc(sizeof(galaxy_struct));

  // Initialise the properties
  gal->ID                = (*unique_ID)++;
  gal->Type              = -1;
  gal->SnapSkipCounter   = 0;
  gal->HaloDescIndex     = -1;
  gal->TreeFlags         = -1;
  gal->Halo              = NULL;
  gal->FirstGalInHalo    = NULL;
  gal->NextGalInHalo     = NULL;
  gal->Next              = NULL;
  gal->MergerTarget      = NULL;
  gal->Len               = 0;
  gal->Mvir              = 0.0;
  gal->dM                = 0.0;
  gal->dMdt              = 0.0;
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

  return gal;
}

static double calculate_Vvir(run_globals_struct *run_globals, halo_struct *halo)
{
  return sqrt(run_globals->G * halo->Mvir / halo->Rvir);
}

void copy_halo_to_galaxy(run_globals_struct *run_globals, halo_struct *halo, galaxy_struct *gal)
{
  gal->Halo            = halo;
  gal->Type            = halo->Type;
  gal->Len             = halo->Len;
  gal->SnapSkipCounter = halo->SnapOffset;
  gal->HaloDescIndex   = halo->DescIndex;
  gal->Mvir            = halo->Mvir;
  gal->Rvir            = halo->Rvir;
  gal->Vvir            = calculate_Vvir(run_globals, halo);
  gal->Vmax            = halo->Vmax;
  gal->TreeFlags       = halo->TreeFlags;
  for (int ii=0; ii<3; ii++)
  {
    gal->Pos[ii]       = halo->Pos[ii];
    gal->Vel[ii]       = halo->Vel[ii];
  }
}
