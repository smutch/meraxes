#include "meraxes.h"
#include <math.h>

galaxy_t* new_galaxy(run_globals_t *run_globals, int *unique_ID)
{
 
  galaxy_t *gal = NULL;

  gal = SID_malloc(sizeof(galaxy_t));

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


void copy_halo_to_galaxy(halo_t *halo, galaxy_t *gal, int snapshot)
{
  gal->id_MBP          = halo->id_MBP;
  gal->Type            = halo->Type;
  gal->Len             = halo->Len;
  gal->SnapSkipCounter = halo->SnapOffset;
  gal->HaloDescIndex   = halo->DescIndex;
  gal->Mvir            = halo->Mvir;
  gal->Rvir            = (double)(halo->Rvir);
  gal->Vvir            = (double)(halo->Vvir);
  gal->Vmax            = halo->Vmax;
  gal->TreeFlags       = halo->TreeFlags;
  for (int ii=0; ii<3; ii++)
  {
    gal->Pos[ii]       = halo->Pos[ii];
    gal->Vel[ii]       = halo->Vel[ii];
  }
}
