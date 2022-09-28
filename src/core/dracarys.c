#include "BrightnessTemperature.h"
#include "ComputePowerSpectrum.h"
#include "ConstructLightcone.h"
#include "debug.h"
#include "galaxies.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "physics/evolve.h"
#include "physics/mergers.h"
#include "physics/reionization.h"
#include "read_halos.h"
#include "reionization.h"
#include "save.h"
#include "tree_flags.h"

static inline bool check_if_valid_host(halo_t* halo)
{
  // We don't want to place new galaxies in any halos with the following flags set...
  int invalid_flags =
    (TREE_CASE_FRAGMENTED_NORMAL | TREE_CASE_FRAGMENTED_NEW | TREE_CASE_FRAGMENTED_EJECTED // TODO: Now marked as other
     | TREE_CASE_FRAGMENTED_STRAYED | TREE_CASE_MERGER); // TODO: Try off and think about closely

  if ((halo->Type == 0) && ((halo->Galaxy == NULL) || check_for_flag(TREE_CASE_MERGER, halo->Galaxy->TreeFlags)) &&
      !(invalid_flags & halo->TreeFlags))
    return true;
  else
    return false;
}

//! Actually run the model
void dracarys()
{
  trees_info_t trees_info;
  halo_t* halo = NULL;
  fof_group_t* fof_group = NULL;
  galaxy_t* gal = NULL;
  galaxy_t* prev_gal = NULL;
  galaxy_t* next_gal = NULL;
  galaxy_t* cur_gal = NULL;
  int i_newhalo;
  int NGal = 0;
  int nout_gals = 0;
  int last_nout_gals = 0;
  int last_snap = 0;
  int kill_counter = 0;
  int i_snap;
  int NSteps = run_globals.params.NSteps;
  int n_store_snapshots = run_globals.NStoreSnapshots;
  halo_t** snapshot_halo = run_globals.SnapshotHalo;
  fof_group_t** snapshot_fof_group = run_globals.SnapshotFOFGroup;
  int** snapshot_index_lookup = run_globals.SnapshotIndexLookup;
  trees_info_t* snapshot_trees_info = run_globals.SnapshotTreesInfo;
  double* LTTime = run_globals.LTTime;
  int NOutputSnaps = run_globals.NOutputSnaps;

  // Find what the last requested output snapshot is
  for (int ii = 0; ii < NOutputSnaps; ii++)
    if (run_globals.ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals.ListOutputSnaps[ii];

  // Prep the output file
  if (!run_globals.params.FlagMCMC) {
    sprintf(run_globals.FNameOut,
            "%s/%s_%d.hdf5",
            run_globals.params.OutputDir,
            run_globals.params.FileNameGalaxies,
            run_globals.mpi_rank);
    prep_hdf5_file();
  }

  // Initialize timer
  timer_info timer;
  timer_start(&timer);

  // Loop through each snapshot
  for (int snapshot = 0; snapshot <= last_snap; snapshot++) {
    int* index_lookup = NULL;
    int merger_counter = 0;
    int new_gal_counter = 0;
    int ghost_counter = 0;

    mlog("", MLOG_MESG);
    mlog("===============================================================", MLOG_MESG);
    mlog("Snapshot %d  (z = %.3f)", MLOG_MESG, snapshot, run_globals.ZZ[snapshot]);
    mlog("===============================================================", MLOG_MESG);

    // Reset book keeping counters
    kill_counter = 0;

    // Read in the halos for this snapshot
    if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC)
      i_snap = snapshot;
    else
      i_snap = 0;

    trees_info = read_halos(snapshot,
                            &(snapshot_halo[i_snap]),
                            &(snapshot_fof_group[i_snap]),
                            &(snapshot_index_lookup[i_snap]),
                            snapshot_trees_info);

    // Set the relevant pointers to this snapshot
    halo = snapshot_halo[i_snap];
    fof_group = snapshot_fof_group[i_snap];
    index_lookup = snapshot_index_lookup[i_snap];

    mlog("Processing snapshot %d (z = %.2f)...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, run_globals.ZZ[snapshot]);

    // Calculate the critical halo mass for cooling
    if ((run_globals.params.Flag_PatchyReion) && (run_globals.params.ReionUVBFlag)) {
      calculate_Mvir_crit(run_globals.ZZ[snapshot]);

      if (run_globals.params.Flag_IncludeLymanWerner)
        calculate_Mvir_crit_MC(run_globals.ZZ[snapshot]);
    }

    // Reset the halo pointers and ghost flags for all galaxies and decrement
    // the snapskip counter
    gal = run_globals.FirstGal;
    while (gal != NULL) {
      gal->Halo = NULL;
      gal->ghost_flag = false;
      gal->SnapSkipCounter--;
      reset_galaxy_properties(gal, snapshot);
      gal = gal->Next;
    }

    // Loop through each galaxy we already have
    gal = run_globals.FirstGal;
    prev_gal = NULL;
    while (gal != NULL) {
      // Get the index of this galaxies descendent halo (which will be the one
      // which exists at this snapshot unless the halo has skipped a snap).
      i_newhalo = gal->HaloDescIndex;

      // If the halo of this galaxy should exist at this snapshot.
      if (gal->SnapSkipCounter <= 0) {
        // If we are subsampling the trees, or for some other reason need to
        // find the corrected halo index, then do so.
        if ((index_lookup) && (i_newhalo > -1) && !(gal->ghost_flag) && (gal->Type < 2))
          i_newhalo = find_original_index(gal->HaloDescIndex, index_lookup, trees_info.n_halos);

        // If this galaxy hasn't been marked for death
        if (i_newhalo > -1) {
          gal->OldType = gal->Type;
          gal->dt = LTTime[gal->LastIdentSnap] - LTTime[snapshot];

          // If this is a central or a satellite
          if (gal->Type < 2)
            connect_galaxy_and_halo(gal, &halo[i_newhalo], &merger_counter);
        } else { // this galaxy has been marked for death
          if (gal->FirstGalInHalo == gal) {
            // We have marked the first galaxy in the halo for death. If there are any
            // other type 2 galaxies in this halo then we must kill them as well...
            // Unfortunately, we don't know if we have alrady processed these
            // remaining galaxies, so we have to just mark them for the moment
            // and then do another check below...
            cur_gal = gal->NextGalInHalo;
            while (cur_gal != NULL) {
              cur_gal->HaloDescIndex = -1;
              cur_gal = cur_gal->NextGalInHalo;
            }
          }
          kill_galaxy(gal, prev_gal, &NGal, &kill_counter);
          gal = prev_gal;
        }
      } else // this galaxy's halo has skipped this snapshot
      {
        // This is a ghost galaxy for this snapshot.
        // We need to count all the other galaxies in this halo as ghosts as
        // well since they won't be reachable by traversing the FOF groups
        cur_gal = gal;
        while (cur_gal != NULL) {
          if (cur_gal->HaloDescIndex > -1) {
            ghost_counter++;
            cur_gal->ghost_flag = true;
          }
          cur_gal = cur_gal->NextGalInHalo;
        }
      }

      // gal may be NULL if we just killed the first galaxy
      if (gal != NULL) {
        prev_gal = gal;
        gal = gal->Next;
      } else
        gal = run_globals.FirstGal;
    }

    // Do one more pass to make sure that we have killed all galaxies which we
    // should have (i.e. satellites in strayed halos etc.)
    prev_gal = NULL;
    next_gal = NULL;
    gal = run_globals.FirstGal;
    while (gal != NULL) {
      if (gal->HaloDescIndex < 0) {
        next_gal = gal->Next;
        kill_galaxy(gal, prev_gal, &NGal, &kill_counter);
        gal = prev_gal;
      }
      prev_gal = gal;
      if (gal != NULL)
        gal = gal->Next;
      else
        gal = next_gal;
    }

    // Store the number of ghost galaxies present at this snapshot
    run_globals.NGhosts = ghost_counter;

    // Incase we ended up removing the last galaxy, update the LastGal pointer
    run_globals.LastGal = prev_gal;

    // Find empty (valid) type 0 halos and place new galaxies in them.
    // Also update the fof_group pointers.
    // Note that we can (and sometimes do) have cases where halos with
    // galaxies have merged into haloes that don't have galaxies.  What do
    // do in this situation is debatable.  If we want to assume that these
    // empty halos could have formed galaxies before the merger event, then
    // this for loop must appear before the follwing while loop.  If we
    // want to assume that these halos wouldn't have formed galaxies then
    // it should come after the while loop...
    for (int i_fof = 0; i_fof < trees_info.n_fof_groups; i_fof++) {
      halo_t* cur_halo = fof_group[i_fof].FirstHalo;
      int total_subhalo_len = 0;

      while (cur_halo != NULL) {
        if (check_if_valid_host(cur_halo))
          create_new_galaxy(snapshot, cur_halo, &NGal, &new_gal_counter, &merger_counter);

        total_subhalo_len += cur_halo->Len;

        cur_halo = cur_halo->NextHaloInFOFGroup;
      }

      fof_group[i_fof].TotalSubhaloLen = total_subhalo_len;
    }

    // Loop through each galaxy and set the merger clocks for new infallers
    // now that all other galaxies have been processed and their halo
    // pointers updated...
    gal = run_globals.FirstGal;

    while (gal != NULL) {
      if ((gal->Type == 2) && (gal->MergerTarget == NULL)) {
        // Set the merger target of the incoming galaxy and initialise the
        // merger clock.  Note that we *increment* the clock immediately
        // after calculating it. This is because we will decrement the clock
        // (by the same amount) when checking for mergers in evolve.c
        gal->MergerTarget = gal->FirstGalInHalo;
        gal->MergTime = calculate_merging_time(gal, snapshot);
        gal->MergTime += gal->dt;
      }
      gal = gal->Next;
    }

    // Calculate the first occupied halo
    for (int i_fof = 0; i_fof < trees_info.n_fof_groups; i_fof++) {
      fof_group[i_fof].FirstOccupiedHalo = NULL;
      halo_t* cur_halo = fof_group[i_fof].FirstHalo;
      while (cur_halo != NULL) {
        if (cur_halo->Galaxy != NULL) {
          fof_group[i_fof].FirstOccupiedHalo = cur_halo;
          break;
        }
        cur_halo = cur_halo->NextHaloInFOFGroup;
      }
    }

    // We finish by copying the halo properties into the galaxy structure of
    // all galaxies with type<2, passively evolving ghosts, and updating the dt
    // values for non-ghosts.
    gal = run_globals.FirstGal;
    while (gal != NULL) {
      if ((gal->Halo == NULL) && (!gal->ghost_flag)) {
        mlog_error("We missed a galaxy during processing!");
#ifdef DEBUG
        mpi_debug_here();
#endif
        ABORT(EXIT_FAILURE);
      }

      if (!gal->ghost_flag)
        gal->dt /= (double)NSteps;
      else
        passively_evolve_ghost(gal, snapshot);

      if ((gal->Type < 2) && (!gal->ghost_flag))
        copy_halo_props_to_galaxy(gal->Halo, gal);

      gal = gal->Next;
    }

#ifdef DEBUG
    check_counts(fof_group, NGal, trees_info.n_fof_groups);
#endif

    if (run_globals.params.Flag_PatchyReion) {
      int ngals_in_slabs = map_galaxies_to_slabs(NGal);
      if (run_globals.params.ReionUVBFlag) {
        assign_Mvir_crit_to_galaxies(ngals_in_slabs, 1);
        if (run_globals.params.Flag_IncludeLymanWerner)
          assign_Mvir_crit_to_galaxies(ngals_in_slabs, 2);
      }
    }

    // Do the physics
    if (NGal > 0)
      nout_gals = evolve_galaxies(fof_group, snapshot, NGal, trees_info.n_fof_groups);
    else
      nout_gals = 0;

    // Add the ghost galaxies into the nout_gals count
    nout_gals += ghost_counter;

    if (run_globals.params.Flag_PatchyReion) {

      if (check_if_reionization_ongoing(snapshot)) {
        if (!run_globals.params.ReionUVBFlag) {
          // We are decoupled, so no need to run 21cmFAST unless we are ouputing this snapshot
          for (int i_out = 0; i_out < NOutputSnaps; i_out++) {
            if (snapshot == run_globals.ListOutputSnaps[i_out]) {
              call_find_HII_bubbles(snapshot, nout_gals, &timer);

              if (run_globals.params.Flag_Compute21cmBrightTemp) {
                ComputeBrightnessTemperatureBox(snapshot);
              }

              if (run_globals.params.Flag_ComputePS) {
                Compute_PS(snapshot);
              }
            }
          }
        } else {

          if (run_globals.params.Flag_IncludeSpinTemp) {
            call_ComputeTs(snapshot, nout_gals, &timer);
          }

          call_find_HII_bubbles(snapshot, nout_gals, &timer);

          if (run_globals.params.Flag_Compute21cmBrightTemp) {
            ComputeBrightnessTemperatureBox(snapshot);
          }

          if (run_globals.params.Flag_ComputePS) {
            Compute_PS(snapshot);
          }

          if (run_globals.params.Flag_ConstructLightcone) {
            ConstructLightcone(snapshot);
          }
        }
      }

      // if we have already created a mapping of galaxies to MPI slabs then we no
      // longer need them as they will need to be re-created for the new halo
      // positions in the next time step
      free(run_globals.reion_grids.galaxy_to_slab_map);
    }

#ifdef DEBUG
    // print some statistics for this snapshot
    MPI_Allreduce(MPI_IN_PLACE, &merger_counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &kill_counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &new_gal_counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &ghost_counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);

    mlog("Newly identified merger events    :: %d", MLOG_MESG, merger_counter);
    mlog("Killed galaxies                   :: %d", MLOG_MESG, kill_counter);
    mlog("Newly created galaxies            :: %d", MLOG_MESG, new_gal_counter);
    mlog("Galaxies in ghost halos           :: %d", MLOG_MESG, ghost_counter);
#endif

    // Write the results if this is a requested snapshot
    if (!run_globals.params.FlagMCMC)
      for (int i_out = 0; i_out < NOutputSnaps; i_out++)
        if (snapshot == run_globals.ListOutputSnaps[i_out])
          write_snapshot(nout_gals, i_out, &last_nout_gals);

    // Update the LastIdentSnap values for non-ghosts
    gal = run_globals.FirstGal;
    while (gal != NULL) {
      if (!gal->ghost_flag)
        gal->LastIdentSnap = snapshot;
      gal = gal->Next;
    }

#ifdef DEBUG
    check_pointers(halo, fof_group, &trees_info);
#endif

    if (run_globals.params.FlagMCMC)
      meraxes_mhysa_hook(run_globals.mhysa_self, snapshot, nout_gals);

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
  }

  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC) {
    // Tidy up counters and galaxies from this iteration
    NGal = 0;
    nout_gals = 0;
    last_nout_gals = 0;

    mlog("Resetting halo->galaxy pointers", MLOG_MESG);
    for (int ii = 0; ii < n_store_snapshots; ii++)
      for (int jj = 0; jj < snapshot_trees_info[ii].n_halos; jj++)
        snapshot_halo[ii][jj].Galaxy = NULL;

    // reset started and finished flags for reionization and reinialize the grids if needed
    if (run_globals.params.Flag_PatchyReion) {
      run_globals.reion_grids.started = 0;
      run_globals.reion_grids.finished = 0;

      init_reion_grids();
    }
  }

  mlog("Freeing galaxies...", MLOG_OPEN);
  gal = run_globals.FirstGal;
  while (gal != NULL) {
    next_gal = gal->Next;
    free(gal);
    gal = next_gal;
  }
  run_globals.FirstGal = NULL;
  mlog("...done", MLOG_CLOSE);

  // Create the master file
  MPI_Barrier(run_globals.mpi_comm);
  if (!run_globals.params.FlagMCMC)
    if (run_globals.mpi_rank == 0)
      create_master_file();
}
