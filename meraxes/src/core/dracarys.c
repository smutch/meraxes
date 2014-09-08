#include <math.h>
#include "meraxes.h"
#include "tree_flags.h"


static void inline turn_off_merger_flag(galaxy_t *gal)
{
  gal->TreeFlags = gal->TreeFlags & (~TREE_CASE_MERGER);
}

static inline bool check_for_flag(int flag, int tree_flags)
{
  if ((tree_flags & flag) == flag)
    return true;
  else
    return false;
}

static inline bool check_for_merger(galaxy_t *gal, halo_t *new_halo)
{
  // if this is marked as a merger in the trees OR the galaxy stellar mass is
  // now greater than the host halo mass, then mark this as a merger
  // if( check_for_flag(TREE_CASE_MERGER, gal->TreeFlags) || (gal->StellarMass > new_halo->Mvir) )
  if (check_for_flag(TREE_CASE_MERGER, gal->TreeFlags))
    return true;
  else
    return false;
}

static inline bool check_if_valid_host(run_globals_t *run_globals, halo_t *halo)
{
  // int invalid_flags = (TREE_CASE_FRAGMENTED_RETURNED
  //     | TREE_CASE_FRAGMENTED_EXCHANGED
  //     | TREE_CASE_FRAGMENTED_STRAYED);
  // | TREE_CASE_STRAYED);
  // | TREE_CASE_SPUTTERED);

  if ((halo->Type == 0)
      && (halo->Galaxy == NULL))
    // && (check_for_flag(invalid_flags, halo->TreeFlags)==false))
    return true;
  else
    return false;
}


static int find_original_index(int index, int *lookup, int n_mappings)
{
  int *pointer  = NULL;
  int new_index = -1;

  pointer = bsearch(&index, lookup, (size_t)n_mappings, sizeof(int), compare_ints);
  if (pointer)
    new_index = (int)(pointer - lookup);

  return new_index;
}


//! Actually run the model
void dracarys(run_globals_t *run_globals)
{
  trees_info_t trees_info;
  halo_t *halo           = NULL;
  fof_group_t *fof_group = NULL;
  galaxy_t *gal          = NULL;
  galaxy_t *prev_gal     = NULL;
  galaxy_t *next_gal     = NULL;
  galaxy_t *cur_gal      = NULL;
  int *index_lookup      = NULL;
  int i_newhalo;
  int NGal            = 0;
  int nout_gals       = 0;
  int last_nout_gals  = 0;
  int last_snap       = 0;
  int kill_counter    = 0;
  int merger_counter  = 0;
  int new_gal_counter = 0;
  int ghost_counter   = 0;
  int i_snap;
  int NSteps                        = run_globals->params.NSteps;
  int n_store_snapshots             = run_globals->NStoreSnapshots;
  halo_t **snapshot_halo            = run_globals->SnapshotHalo;
  fof_group_t **snapshot_fof_group  = run_globals->SnapshotFOFGroup;
  int **snapshot_index_lookup       = run_globals->SnapshotIndexLookup;
  trees_info_t *snapshot_trees_info = run_globals->SnapshotTreesInfo;
  double *LTTime = run_globals->LTTime;

  // Find what the last requested output snapshot is
  for (int ii = 0; ii < NOUT; ii++)
    if (run_globals->ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals->ListOutputSnaps[ii];

  // Prep the output file
  sprintf(run_globals->FNameOut, "%s/%s_%d.hdf5", run_globals->params.OutputDir, run_globals->params.FileNameGalaxies, SID.My_rank);
  prep_hdf5_file(run_globals);

  // Loop through each snapshot
  for (int snapshot = 0; snapshot <= last_snap; snapshot++)
  {
    // Reset book keeping counters
    kill_counter    = 0;
    merger_counter  = 0;
    new_gal_counter = 0;
    ghost_counter   = 0;

    // Read in the halos for this snapshot
    if (run_globals->params.FlagInteractive)
      i_snap = snapshot;
    else
      i_snap = 0;

    trees_info = read_halos(run_globals, snapshot, &(snapshot_halo[i_snap]), &(snapshot_fof_group[i_snap]), &(snapshot_index_lookup[i_snap]), snapshot_trees_info);

    // Set the relevant pointers to this snapshot
    halo         = snapshot_halo[i_snap];
    fof_group    = snapshot_fof_group[i_snap];
    index_lookup = snapshot_index_lookup[i_snap];

    SID_log("Processing snapshot %d (z=%.2f)...", SID_LOG_OPEN | SID_LOG_TIMER, snapshot, run_globals->ZZ[snapshot]);

#ifdef USE_TOCF
    // Calculate the critical halo mass for cooling
    if ((run_globals->params.TOCF_Flag) && (tocf_params.uvb_feedback))
      calculate_Mvir_crit(run_globals, run_globals->ZZ[snapshot]);
#endif

    // Reset the halo pointers and ghost flags for all galaxies and decrement
    // the snapskip counter
    gal = run_globals->FirstGal;
    while (gal != NULL)
    {
      gal->Halo       = NULL;
      gal->ghost_flag = false;
      gal->SnapSkipCounter--;
      reset_galaxy_properties(run_globals, gal, snapshot);
      gal = gal->Next;
    }

    gal      = run_globals->FirstGal;
    prev_gal = NULL;
    while (gal != NULL)
    {
      i_newhalo = gal->HaloDescIndex;

      if (gal->SnapSkipCounter <= 0)
      {
        if ((index_lookup) && (i_newhalo > -1) && !(gal->ghost_flag) && (gal->Type < 2))
          i_newhalo = find_original_index(gal->HaloDescIndex, index_lookup, trees_info.n_halos);

        if (i_newhalo > -1)
        {
          gal->OldType = gal->Type;
          gal->dt      = LTTime[gal->LastIdentSnap] - LTTime[snapshot];
          if (gal->Type < 2)
          {
            if (check_for_merger(gal, &(halo[i_newhalo])))
            {
              // Here we have a new merger...  Mark it and deal with it below.
              gal->Type = 999;
              merger_counter++;

              // if this was marked as a halo merger in the input trees then
              // update the halo pointer correspondingly and turn off the
              // merger flag
              if (check_for_flag(TREE_CASE_MERGER, gal->TreeFlags))
              {
                gal->Halo = &(halo[i_newhalo]);
                turn_off_merger_flag(gal);
              }
              // otherwise, set the new halo to be the central halo of the corresponding FOF group
              else
                gal->Halo = halo[i_newhalo].FOFGroup->FirstHalo;
            }
            else
            {
              // Here we have the simplest case where a galaxy continues along in it's halo...
              gal->Halo = &(halo[i_newhalo]);
              assign_galaxy_to_halo(gal, &(halo[i_newhalo]));

              // Loop through all of the other galaxies in this halo and set their Halo pointer
              cur_gal = gal->NextGalInHalo;
              while (cur_gal != NULL)
              {
                cur_gal->Halo = &(halo[i_newhalo]);
                cur_gal       = cur_gal->NextGalInHalo;
              }
            }
          }
        }
        else    // this galaxy has been marked for death
        {
          if (gal->FirstGalInHalo == gal)
          {
            // We have marked the first galaxy in the halo for death. If there are any
            // other type 2 galaxies in this halo then we must kill them as well...
            // Unfortunately, we don't know if we have alrady processed these
            // remaining galaxies, so we have to just mark them for the moment
            // and then do another check below...
            cur_gal = gal->NextGalInHalo;
            while (cur_gal != NULL)
            {
              cur_gal->HaloDescIndex = -1;
              cur_gal                = cur_gal->NextGalInHalo;
            }
          }
          kill_galaxy(run_globals, gal, prev_gal, &NGal, &kill_counter);
          gal = prev_gal;
        }
      }
      else    // this galaxy's halo has skipped this snapshot
      {
        // This is a ghost galaxy for this snapshot.
        // We need to count all the other galaxies in this halo as ghosts as
        // well since they won't be reachable by traversing the FOF groups
        cur_gal = gal;
        while (cur_gal != NULL)
        {
          if (cur_gal->HaloDescIndex > -1)
          {
            ghost_counter++;
            cur_gal->ghost_flag = true;
          }
          cur_gal = cur_gal->NextGalInHalo;
        }
      }

      // gal may be NULL if we just killed the first galaxy
      if (gal != NULL)
      {
        prev_gal = gal;
        gal      = gal->Next;
      }
      else
        gal = run_globals->FirstGal;
    }

    // Do one more pass to make sure that we have killed all galaxies which we
    // should have (i.e. satellites in strayed halos etc.)
    prev_gal = NULL;
    gal      = run_globals->FirstGal;
    while (gal != NULL)
    {
      if (gal->HaloDescIndex < 0)
      {
        kill_galaxy(run_globals, gal, prev_gal, &NGal, &kill_counter);
        gal = prev_gal;
      }
      prev_gal = gal;
      gal      = gal->Next;
    }

    // Store the number of ghost galaxies present at this snapshot
    run_globals->NGhosts = ghost_counter;

    // Incase we ended up removing the last galaxy, update the LastGal pointer
    run_globals->LastGal = prev_gal;

    // Find empty (valid) type 0 halos and place new galaxies in them.
    // Also update the fof_group pointers.
    for (int i_fof = 0; i_fof < trees_info.n_fof_groups; i_fof++)
    {
      halo_t *first_occupied_halo = NULL;
      halo_t *cur_halo = fof_group[i_fof].FirstHalo;
      int total_subhalo_len = 0;

      while (cur_halo != NULL)
      {
        if (check_if_valid_host(run_globals, cur_halo))
          create_new_galaxy(run_globals, snapshot, cur_halo, &NGal, &new_gal_counter);

        if ((first_occupied_halo == NULL) && (cur_halo->Galaxy != NULL))
          first_occupied_halo = cur_halo;

        total_subhalo_len += cur_halo->Len;

        cur_halo = cur_halo->NextHaloInFOFGroup;
      }

      fof_group[i_fof].FirstOccupiedHalo = first_occupied_halo;
      fof_group[i_fof].TotalSubhaloLen = total_subhalo_len;
    }

    // Loop through each galaxy and deal with HALO mergers now that all other
    // galaxies have been processed and their halo pointers updated...
    gal = run_globals->FirstGal;
    while (gal != NULL)
    {
      if (gal->Type == 999)
      {
        if (gal->Halo->Galaxy == NULL)
        {
          // Here we have a halo with a galaxy that has just merged into an
          // empty halo.  From the point of view of the model, this isn't
          // actually a merger and so we need to catch these cases...
          gal->Halo->Galaxy = gal;
          gal->Type         = gal->Halo->Type;
          cur_gal           = gal->NextGalInHalo;
          while (cur_gal != NULL)
          {
            cur_gal->Halo = gal->Halo;
            cur_gal       = cur_gal->NextGalInHalo;
          }
        }
        else    // there are galaxies in the halo being merged into
        {
          // Remember that we have already set the galaxy's halo pointer so we
          // don't need to do that here.

          gal->Type = 2;

          // Add the incoming galaxy to the end of the halo's linked list
          cur_gal = gal->Halo->Galaxy;
          while (cur_gal != NULL)
          {
            prev_gal = cur_gal;
            cur_gal  = cur_gal->NextGalInHalo;
          }
          prev_gal->NextGalInHalo = gal;

          // Update the FirstGalInHalo pointer.
          gal->FirstGalInHalo = gal->Halo->Galaxy;

          // Loop through and update the FirstGalInHalo and Halo pointers of any other
          // galaxies that are attached to the incoming galaxy
          cur_gal = gal->NextGalInHalo;
          while (cur_gal != NULL)
          {
            cur_gal->FirstGalInHalo = gal->FirstGalInHalo;
            cur_gal->Halo           = gal->Halo;
            cur_gal                 = cur_gal->NextGalInHalo;
          }

          if (gal->FirstGalInHalo == NULL)
            SID_log_warning("Just set gal->FirstGalInHalo = NULL!", SID_LOG_COMMENT);

          // Set the merger target of the incoming galaxy and initialise the
          // merger clock.  Note that we *increment* the clock imemdiately
          // after calculating it. This is because we will decrement the clock
          // (by the same amount) when checking for mergers in evolve.c
          gal->MergerTarget = gal->FirstGalInHalo;
          gal->MergTime     = calculate_merging_time(run_globals, gal, snapshot);
          gal->MergTime    += gal->dt;
        }
      }

      gal = gal->Next;
    }

    // We finish by copying the halo properties into the galaxy structure of
    // all galaxies with type<2, passively evolving ghosts, and updating the dt
    // values for non-ghosts.
    gal = run_globals->FirstGal;
    while (gal != NULL)
    {
      if ((gal->Halo == NULL) && (!gal->ghost_flag))
      {
        SID_log_error("We missed a galaxy during processing!");
#ifdef DEBUG
        mpi_debug_here();
#endif
        ABORT(EXIT_FAILURE);
      }

      if (!gal->ghost_flag)
        gal->dt /= (double)NSteps;
      else
        passively_evolve_ghost(run_globals, gal, snapshot);

      if ((gal->Type < 2) && (!gal->ghost_flag))
        copy_halo_to_galaxy(gal->Halo, gal, snapshot);

      gal = gal->Next;
    }

#ifdef DEBUG
    check_counts(run_globals, fof_group, NGal, trees_info.n_fof_groups);
#endif

    // Do the physics
    if (NGal > 0)
      nout_gals = evolve_galaxies(run_globals, fof_group, snapshot, NGal, trees_info.n_fof_groups);
    else
      nout_gals = 0;

    // Add the ghost galaxies into the nout_gals count
    nout_gals += ghost_counter;

#ifdef USE_TOCF
    if (run_globals->params.TOCF_Flag)
    {
      if (!tocf_params.uvb_feedback)
      {
        // We are decoupled, so no need to run 21cmFAST unless we are ouputing this snapshot
        for (int i_out = 0; i_out < NOUT; i_out++)
          if (snapshot == run_globals->ListOutputSnaps[i_out])
            call_find_HII_bubbles(run_globals, snapshot, nout_gals);
      }
      else
        call_find_HII_bubbles(run_globals, snapshot, nout_gals);
    }
#endif

#ifdef DEBUG
    // print some statistics for this snapshot
    SID_Allreduce(SID_IN_PLACE, &merger_counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
    SID_Allreduce(SID_IN_PLACE, &kill_counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
    SID_Allreduce(SID_IN_PLACE, &new_gal_counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
    SID_Allreduce(SID_IN_PLACE, &ghost_counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);

    SID_log("Newly identified merger events    :: %d", SID_LOG_COMMENT, merger_counter);
    SID_log("Killed galaxies                   :: %d", SID_LOG_COMMENT, kill_counter);
    SID_log("Newly created galaxies            :: %d", SID_LOG_COMMENT, new_gal_counter);
    SID_log("Galaxies in ghost halos           :: %d", SID_LOG_COMMENT, ghost_counter);
#endif

    // Write the results if this is a requested snapshot
    for (int i_out = 0; i_out < NOUT; i_out++)
      if (snapshot == run_globals->ListOutputSnaps[i_out])
        write_snapshot(run_globals, nout_gals, i_out, &last_nout_gals, &trees_info);

#ifdef USE_TOCF
    if (run_globals->params.TOCF_Flag)
      check_if_reionization_complete(run_globals);
#endif

    // Update the LastIdentSnap values for non-ghosts
    gal = run_globals->FirstGal;
    while (gal != NULL)
    {
      if (!gal->ghost_flag)
        gal->LastIdentSnap = snapshot;
      gal = gal->Next;
    }

#ifdef DEBUG
    check_pointers(run_globals, halo, fof_group, &trees_info);
#endif

    SID_log("...done", SID_LOG_CLOSE);
  }

  if (run_globals->params.FlagInteractive)
  {
    // Tidy up counters and galaxies from this iteration
    NGal           = 0;
    nout_gals      = 0;
    last_nout_gals = 0;

    SID_log("Resetting halo->galaxy pointers", SID_LOG_COMMENT);
    for (int ii = 0; ii < n_store_snapshots; ii++)
      for (int jj = 0; jj < snapshot_trees_info[ii].n_halos; jj++)
        snapshot_halo[ii][jj].Galaxy = NULL;
  }

  SID_log("Freeing galaxies...", SID_LOG_OPEN);
  gal = run_globals->FirstGal;
  while (gal != NULL)
  {
    next_gal = gal->Next;
    SID_free(SID_FARG gal);
    gal = next_gal;
  }
  run_globals->FirstGal = NULL;
  SID_log(" ...done", SID_LOG_CLOSE);

  // Create the master file
  SID_Barrier(SID.COMM_WORLD);
  if (SID.My_rank == 0)
    create_master_file(run_globals);
}


