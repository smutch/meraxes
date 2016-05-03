#include "meraxes.h"
#include <assert.h>
#include <math.h>

galaxy_t* new_galaxy(run_globals_t *run_globals, int snapshot, int halo_ID)
{
  galaxy_t *gal = NULL;

  gal = SID_malloc(sizeof(galaxy_t));

  // Initialise the properties
  gal->id_MBP             = 0;
  gal->ID                 = snapshot * 10000000 + halo_ID;
  gal->Type               = -1;
  gal->OldType            = -1;
  gal->SnapSkipCounter    = 0;
  gal->HaloDescIndex      = -1;
  gal->TreeFlags          = -1;
  gal->LastIdentSnap      = -1;
  gal->Halo               = NULL;
  gal->FirstGalInHalo     = NULL;
  gal->NextGalInHalo      = NULL;
  gal->Next               = NULL;
  gal->MergerTarget       = NULL;
  gal->Len                = 0;
  gal->MaxLen             = 0;
  gal->dt                 = 0.0;
  gal->Mvir               = 0.0;
  gal->Rvir               = 0.0;
  gal->Vvir               = 0.0;
  gal->Vmax               = 0.0;
  gal->Spin               = 0.0;
  gal->DiskScaleLength    = 0.0;
  gal->HotGas             = 0.0;
  gal->MetalsHotGas       = 0.0;
  gal->ColdGas            = 0.0;
  gal->MetalsColdGas      = 0.0;
  gal->EjectedGas         = 0.0;
  gal->MetalsEjectedGas   = 0.0;
  gal->Mcool              = 0.0;
  gal->StellarMass        = 0.0;
  gal->GrossStellarMass   = 0.0;
  gal->MetalsStellarMass  = 0.0;
  gal->mwmsa_num          = 0.0;
  gal->mwmsa_denom        = 0.0;
  gal->BlackHoleMass      = 0.0;
  gal->MaxReheatFrac      = 0.0;
  gal->MaxEjectFrac       = 0.0;
  gal->Sfr                = 0.0;
  gal->Cos_Inc            = gsl_rng_uniform(run_globals->random_generator);
  gal->MergTime           = 99999.9;
  gal->BaryonFracModifier = 1.0;
  gal->MvirCrit           = 0.0;
  gal->PhysicsFlags       = 0;
  gal->MergerStartRadius  = 0.0;
  gal->MergerStartRadius  = 0.0;
  gal->MergSnap           = 0;

  for (int ii = 0; ii < 3; ii++)
  {
    gal->Pos[ii] = -99999.9;
    gal->Vel[ii] = -99999.9;
  }

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
    gal->NewStars[ii] = 0.0;

  gal->output_index = -1;
  gal->ghost_flag   = false;

  init_luminosities(run_globals, gal);

  return gal;
}

void copy_halo_to_galaxy(run_globals_t *run_globals, halo_t *halo, galaxy_t *gal, int snapshot)
{
  gal->id_MBP          = halo->id_MBP;
  gal->Type            = halo->Type;
  gal->Len             = halo->Len;
  gal->SnapSkipCounter = halo->SnapOffset;
  gal->HaloDescIndex   = halo->DescIndex;
  gal->Mvir            = halo->Mvir;
  gal->Rvir            = (double)(halo->Rvir);
  gal->Vvir            = (double)(halo->Vvir);
  gal->TreeFlags       = halo->TreeFlags;
  gal->Spin            = calculate_spin_param(halo);

  double sqrt_2 = 1.414213562;
  if (gal->Type == 0)
  {
    gal->Vmax            = halo->Vmax;
    gal->DiskScaleLength = gal->Spin * gal->Rvir / sqrt_2;
  }
  else 
  {
    if (!run_globals->params.physics.Flag_FixVmaxOnInfall)
      gal->Vmax = halo->Vmax;
    if (!run_globals->params.physics.Flag_FixDiskRadiusOnInfall)
      gal->DiskScaleLength = gal->Spin * gal->Rvir / sqrt_2;
  }

  for (int ii = 0; ii < 3; ii++)
  {
    gal->Pos[ii] = halo->Pos[ii];
    gal->Vel[ii] = halo->Vel[ii];
  }

  // record the maximum Len value if necessary
  if (halo->Len > gal->MaxLen)
    gal->MaxLen = halo->Len;
}

void reset_galaxy_properties(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{
  // Here we reset any galaxy properties which are calculated on a snapshot by
  // snapshot basis.
  gal->Mcool              = 0.0;
  gal->Sfr                = 0.0;
  gal->BaryonFracModifier = 1.0;
  gal->MvirCrit           = 0.0;
  gal->PhysicsFlags       = 0;

  // update the stellar mass weighted mean age values
  assert(snapshot > 0);
  gal->mwmsa_denom += gal->NewStars[N_HISTORY_SNAPS - 1];
  gal->mwmsa_num   += gal->NewStars[N_HISTORY_SNAPS - 1] * run_globals->LTTime[snapshot - N_HISTORY_SNAPS];

  // roll over the baryonic history arrays
  for (int ii = N_HISTORY_SNAPS - 1; ii > 0; ii--)
    gal->NewStars[ii] = gal->NewStars[ii - 1];

  gal->NewStars[0] = 0.0;

  gal->MaxEjectFrac = 0.0;
  gal->MaxReheatFrac = 0.0;
}


void assign_galaxy_to_halo(galaxy_t *gal, halo_t *halo)
{
  if (halo->Galaxy == NULL)
    halo->Galaxy = gal;
  else
  {
    SID_log_error("Trying to assign first galaxy to a halo which already has a first galaxy!");
#ifdef DEBUG
    mpi_debug_here();
#endif
    ABORT(EXIT_FAILURE);
  }
}


void create_new_galaxy(
  run_globals_t *run_globals,
  int            snapshot,
  halo_t        *halo,
  int           *NGal,
  int           *new_gal_counter)
{
  galaxy_t *gal;

  gal       = new_galaxy(run_globals, snapshot, halo->ID);
  gal->Halo = halo;

  if (snapshot > 0)
    gal->LastIdentSnap = snapshot - 1;
  else
    gal->LastIdentSnap = snapshot;

  assign_galaxy_to_halo(gal, halo);

  if (run_globals->LastGal != NULL)
    run_globals->LastGal->Next = gal;
  else
    run_globals->FirstGal = gal;

  run_globals->LastGal = gal;
  gal->FirstGalInHalo  = gal;
  gal->dt              = run_globals->LTTime[gal->LastIdentSnap] - run_globals->LTTime[snapshot];
  *NGal                = *NGal + 1;
  *new_gal_counter     = *new_gal_counter + 1;
}


void kill_galaxy(
  run_globals_t *run_globals,
  galaxy_t      *gal,
  galaxy_t      *prev_gal,
  int           *NGal,
  int           *kill_counter)
{
  galaxy_t *cur_gal;

  // Remove it from the global linked list
  if (prev_gal != NULL)
    prev_gal->Next = gal->Next;
  else
    run_globals->FirstGal = gal->Next;

  cur_gal = gal->FirstGalInHalo;

  if (cur_gal != gal)
  {
    // If it is a type 2 then also remove it from the linked list of galaxies in its halo
    while ((cur_gal->NextGalInHalo != NULL) && (cur_gal->NextGalInHalo != gal))
      cur_gal = cur_gal->NextGalInHalo;
    cur_gal->NextGalInHalo = gal->NextGalInHalo;
  }
  else
  {
    // If it is a type 0 or 1 (i.e. first galaxy in it's halo) and there are
    // other galaxies in this halo, reset the FirstGalInHalo pointer so that
    // the satellites can be killed later
    cur_gal = gal->NextGalInHalo;
    while (cur_gal != NULL)
    {
      cur_gal->FirstGalInHalo = gal->NextGalInHalo;
      cur_gal                 = cur_gal->NextGalInHalo;
    }
  }

  // Finally deallocated the galaxy and decrement any necessary counters
  SID_free(SID_FARG gal);
  *NGal         = *NGal - 1;
  *kill_counter = *kill_counter + 1;
}

