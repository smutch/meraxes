#include <assert.h>

#include "galaxies.h"
#include "magnitudes.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "tree_flags.h"
#include "virial_properties.h"

galaxy_t* new_galaxy(int snapshot, unsigned long halo_ID)
{
  galaxy_t* gal = malloc(sizeof(galaxy_t));

  // Initialise the properties
  gal->ID = (unsigned long)(snapshot * 1e10 + halo_ID);
  gal->Type = -1;
  gal->OldType = -1;
  gal->SnapSkipCounter = 0;
  gal->HaloDescIndex = -1;
  gal->TreeFlags = 0;
  gal->LastIdentSnap = -1;
  gal->Halo = NULL;
  gal->FirstGalInHalo = NULL;
  gal->NextGalInHalo = NULL;
  gal->Next = NULL;
  gal->MergerTarget = NULL;
  gal->Len = 0;
  gal->MaxLen = 0;
  gal->dt = 0.0;
  gal->Mvir = 0.0;
  gal->Rvir = 0.0;
  gal->Vvir = 0.0;
  gal->Vmax = 0.0;
  gal->Spin = 0.0;
  gal->DiskScaleLength = 0.0;
  gal->HotGas = 0.0;
  gal->MetalsHotGas = 0.0;
  gal->ColdGas = 0.0;
  gal->MetalsColdGas = 0.0;
  gal->H2Frac = 0.0;
  gal->H2Mass = 0.0;
  gal->HIMass = 0.0;
  gal->EjectedGas = 0.0;
  gal->MetalsEjectedGas = 0.0;
  gal->Mcool = 0.0;
  gal->Rcool = 0.0;
  gal->StellarMass = 0.0;
  gal->GrossStellarMass = 0.0;
  gal->Fesc = 1.0;
  gal->FescWeightedGSM = 0.0;
  gal->MetalsStellarMass = 0.0;
  gal->mwmsa_num = 0.0;
  gal->mwmsa_denom = 0.0;
  gal->BlackHoleMass = run_globals.params.physics.BlackHoleSeed;
  gal->FescBH = 1.0;
  gal->BHemissivity = 0.0;
  gal->EffectiveBHM = 0.0;
  gal->BlackHoleAccretedHotMass = 0.0;
  gal->BlackHoleAccretedColdMass = 0.0;
  gal->BlackHoleAccretingColdMass = 0.0;
  gal->Sfr = 0.0;
  gal->Cos_Inc = gsl_rng_uniform(run_globals.random_generator);
  gal->MergTime = 99999.9;
  gal->BaryonFracModifier = 1.0;
  gal->FOFMvirModifier = 1.0;
  gal->MvirCrit = 0.0;
  gal->MvirCrit_MC = 0.0;
  gal->MergerBurstMass = 0.0;
  gal->MergerStartRadius = 0.0;

#if USE_MINI_HALOS
  gal->StellarMass_II = 0.;
  gal->StellarMass_III = 0.;
  gal->GrossStellarMassIII = 0.0;
  gal->FescIII = 1.0;
  gal->FescIIIWeightedGSM = 0.0;
  gal->Remnant_Mass = 0.;
  gal->Metal_Probability = 0.0;
  gal->Metals_IGM = 0.0;
  gal->Gas_IGM = 0.0;
  gal->Metallicity_IGM = -50.0;
  gal->RmetalBubble = 0.0;
  gal->PrefactorBubble = 0.0;
  gal->TimeBubble = 0.0;
  gal->AveBubble = 0.;
  gal->MaxBubble = 0.;
  gal->Flag_ExtMetEnr = 0;
  gal->GalMetal_Probability = gsl_rng_uniform(run_globals.random_generator);

  if (run_globals.params.Flag_IncludeMetalEvo ==
      false) // If you don't have the external metal enrichment all galaxies will start as pristine (Pop.III forming)
    gal->Galaxy_Population = 3;
#else // If you are not computing PopIII , all the galaxies are PopII. Again you need to initialize the variable
      // otherwise star_formation.c will fail!
  gal->Galaxy_Population = 2;
#endif

  for (int ii = 0; ii < 3; ii++) {
    gal->Pos[ii] = (float)-99999.9;
    gal->Vel[ii] = (float)-99999.9;
  }

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++) {
    gal->NewStars[ii] = 0.0;
#if USE_MINI_HALOS
    gal->NewStars_II[ii] = 0.0;
    gal->NewStars_III[ii] = 0.0;
    if (run_globals.params.Flag_IncludeMetalEvo) {
      gal->Prefactor[ii] = 0.0;
      gal->Times[ii] = 0.0;
      gal->Radii[ii] = 0.0;
    }
#endif
  }

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
    gal->NewMetals[ii] = 0.0;

  gal->output_index = -1;
  gal->ghost_flag = false;

#ifdef CALC_MAGS
  init_luminosities(gal);
#endif

  return gal;
}

void copy_halo_props_to_galaxy(halo_t* halo, galaxy_t* gal)
{
  gal->Type = halo->Type;
  gal->Len = halo->Len;
  gal->SnapSkipCounter = halo->SnapOffset;
  gal->HaloDescIndex = halo->DescIndex;
  gal->Mvir = halo->Mvir;
  gal->Rvir = halo->Rvir;
  gal->Vvir = halo->Vvir;
  gal->TreeFlags = halo->TreeFlags;
  gal->Spin = calculate_spin_param(halo);
  gal->FOFMvirModifier = halo->FOFGroup->FOFMvirModifier;

  double sqrt_2 = 1.414213562;
  if (gal->Type == 0) {
    gal->Vmax = halo->Vmax;
    gal->DiskScaleLength = gal->Spin * gal->Rvir / sqrt_2;
  } else {
    if (!run_globals.params.physics.Flag_FixVmaxOnInfall)
      gal->Vmax = halo->Vmax;
    if (!run_globals.params.physics.Flag_FixDiskRadiusOnInfall)
      gal->DiskScaleLength = gal->Spin * gal->Rvir / sqrt_2;
  }

  for (int ii = 0; ii < 3; ii++) {
    gal->Pos[ii] = halo->Pos[ii];
    gal->Vel[ii] = halo->Vel[ii];
  }

  // record the maximum Len value if necessary
  if (halo->Len > gal->MaxLen)
    gal->MaxLen = halo->Len;
}

void reset_galaxy_properties(galaxy_t* gal, int snapshot)
{
  // Here we reset any galaxy properties which are calculated on a snapshot by
  // snapshot basis.
  gal->Sfr = 0.0;
  gal->Mcool = 0.0;
  gal->Rcool = 0.0;
  gal->MvirCrit = 0.0;
  gal->MvirCrit_MC = 0.0;
  gal->BHemissivity = 0.0;
  gal->BaryonFracModifier = 1.0;
  gal->FOFMvirModifier = 1.0;
  gal->BlackHoleAccretedHotMass = 0.0;
  gal->BlackHoleAccretedColdMass = 0.0;

  // Update the stellar mass weighted mean age values.  This only needs to be
  // done for snapshots shich are passing out of what we are able to track
  // with N_HISTORY_SNAPS.
  assert(snapshot > 0);
  if (snapshot >= N_HISTORY_SNAPS) {
    gal->mwmsa_denom += gal->NewStars[N_HISTORY_SNAPS - 1];
    gal->mwmsa_num += gal->NewStars[N_HISTORY_SNAPS - 1] * run_globals.LTTime[snapshot - N_HISTORY_SNAPS];
  }

  // roll over the baryonic history arrays
  for (int ii = N_HISTORY_SNAPS - 1; ii > 0; ii--) {
    gal->NewStars[ii] = gal->NewStars[ii - 1];
#if USE_MINI_HALOS
    gal->NewStars_II[ii] = gal->NewStars_II[ii - 1];
    gal->NewStars_III[ii] = gal->NewStars_III[ii - 1];
#endif
  }

  for (int ii = N_HISTORY_SNAPS - 1; ii > 0; ii--)
    gal->NewMetals[ii] = gal->NewMetals[ii - 1];

  gal->NewStars[0] = 0.0;
#if USE_MINI_HALOS
  gal->NewStars_II[0] = 0.0;
  gal->NewStars_III[0] = 0.0;
#endif
  gal->NewMetals[0] = 0.0;
}

static void push_galaxy_to_halo(galaxy_t* gal, halo_t* halo)
{
  if (halo->Galaxy == NULL)
    halo->Galaxy = gal;
  else {
    // Walk the galaxy list for the halo to find the end and then link the
    // new galaxy
    galaxy_t* cur_gal = halo->Galaxy;
    galaxy_t* prev_gal = cur_gal;
    while (cur_gal != NULL) {
      prev_gal = cur_gal;
      cur_gal->Halo = halo;
      cur_gal = cur_gal->NextGalInHalo;
    }
    prev_gal->NextGalInHalo = gal;
  }

  // Loop through the new galaxy, and all galaxies attached to it, and set
  // the first galaxy in halo and halo pointer
  gal = halo->Galaxy;
  while (gal != NULL) {
    gal->Halo = halo;
    gal->FirstGalInHalo = halo->Galaxy;
    gal = gal->NextGalInHalo;
  }
}

void connect_galaxy_and_halo(galaxy_t* gal, halo_t* halo, int* merger_counter)
{

  if (halo->Galaxy == NULL)
    push_galaxy_to_halo(gal, halo);
  else {
    // There is already a galaxy been assigned to this halo.  That means we
    // have a merger. Now we need to work out which galaxy is merging into
    // which.

    assert(merger_counter != NULL);
    (*merger_counter)++;

    galaxy_t* parent = NULL;
    galaxy_t* infaller = NULL;

    switch (run_globals.params.TreesID) {
      case GBPTREES_TREES:
        // For gbpTrees, we have the merger flags to give us guidance.  Let's use them...
        // TODO: Make sure I don't need to turn off the merger flag...
        parent = check_for_flag(TREE_CASE_MERGER, gal->TreeFlags) ? halo->Galaxy : gal;
        infaller = halo->Galaxy == parent ? gal : halo->Galaxy;
        assert(check_for_flag(TREE_CASE_MERGER, infaller->TreeFlags));
        break;

      case VELOCIRAPTOR_TREES:
        // For VELOCIraptor we have some guidance in the form of the progenitor indices.

        if (check_for_flag(TREE_CASE_NO_PROGENITORS, halo->TreeFlags)) {
          // The host halo has been marked as having no progenitors.
          // Since we have a merger though, their are clearly halos which
          // think this is their descendant.  In this case, none of the
          // halo progenitors are deemed to be good enough matches and so
          // we can't use the pointers to select the main progenitor.
          // Let's use the mass instead in this case.
          //
          // N.B. We haven't yet copied any halo properties into the
          // galaxies and so the galaxy.Mvir values still correspond to
          // the previous snapshot.
          parent = gal->Mvir > halo->Galaxy->Mvir ? gal : halo->Galaxy;
        } else {
          // Here we can use the pointers.
          parent = gal->HaloDescIndex == halo->ProgIndex ? gal : halo->Galaxy;
        }

        infaller = halo->Galaxy == parent ? gal : halo->Galaxy;
        break;

      default:
        mlog_error("Unrecognised input trees identifier (TreesID).");
        break;
    }

    infaller->Type = 2;
    // Make sure the halo is pointing to the right galaxy
    if (parent != halo->Galaxy)
      halo->Galaxy = parent;

    // Add the incoming galaxy to the end of the halo's linked list
    push_galaxy_to_halo(infaller, halo);
  }
}

void create_new_galaxy(int snapshot, halo_t* halo, int* NGal, int* new_gal_counter, int* merger_counter)
{
  galaxy_t* gal;

  gal = new_galaxy(snapshot, halo->ID);
  gal->Halo = halo;

  if (snapshot > 0)
    gal->LastIdentSnap = snapshot - 1;
  else
    gal->LastIdentSnap = snapshot;

  connect_galaxy_and_halo(gal, halo, merger_counter);

  if (run_globals.LastGal != NULL)
    run_globals.LastGal->Next = gal;
  else
    run_globals.FirstGal = gal;

  run_globals.LastGal = gal;
  gal->dt = run_globals.LTTime[gal->LastIdentSnap] - run_globals.LTTime[snapshot];
  *NGal = *NGal + 1;
  *new_gal_counter = *new_gal_counter + 1;
}

void kill_galaxy(galaxy_t* gal, galaxy_t* prev_gal, int* NGal, int* kill_counter)
{
  galaxy_t* cur_gal;

  // Remove it from the global linked list
  if (prev_gal != NULL)
    prev_gal->Next = gal->Next;
  else
    run_globals.FirstGal = gal->Next;

  cur_gal = gal->FirstGalInHalo;

  if (cur_gal != gal) {
    // If it is a type 2 then also remove it from the linked list of galaxies in its halo
    while ((cur_gal->NextGalInHalo != NULL) && (cur_gal->NextGalInHalo != gal))
      cur_gal = cur_gal->NextGalInHalo;
    cur_gal->NextGalInHalo = gal->NextGalInHalo;
  } else {
    // If it is a type 0 or 1 (i.e. first galaxy in it's halo) and there are
    // other galaxies in this halo, reset the FirstGalInHalo pointer so that
    // the satellites can be killed later
    cur_gal = gal->NextGalInHalo;
    while (cur_gal != NULL) {
      cur_gal->FirstGalInHalo = gal->NextGalInHalo;
      cur_gal = cur_gal->NextGalInHalo;
    }
  }

  // Finally deallocated the galaxy and decrement any necessary counters
  free(gal);
  *NGal = *NGal - 1;
  *kill_counter = *kill_counter + 1;
}