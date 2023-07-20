#include <assert.h>
#include <gsl/gsl_sort_int.h>
#include <hdf5_hl.h>

#include "meraxes.h"
#include "misc_tools.h"
#include "modifiers.h"
#include "read_halos.h"

inline static void update_pointers_from_offsets(int n_fof_groups_kept,
                                                fof_group_t* fof_group,
                                                const size_t* fof_FirstHalo_os,
                                                int n_halos_kept,
                                                halo_t* halo,
                                                const size_t* halo_FOFGroup_os,
                                                const size_t* halo_NextHaloInFOFGroup_os)
{
  for (int ii = 0; ii < n_halos_kept; ii++) {
    halo[ii].FOFGroup = &(fof_group[halo_FOFGroup_os[ii]]);
    if (halo_NextHaloInFOFGroup_os[ii] != (size_t)-1)
      halo[ii].NextHaloInFOFGroup = &(halo[halo_NextHaloInFOFGroup_os[ii]]);
    else
      halo[ii].NextHaloInFOFGroup = NULL;
  }
  for (int ii = 0; ii < n_fof_groups_kept; ii++)
    fof_group[ii].FirstHalo = &(halo[fof_FirstHalo_os[ii]]);
}

static fof_group_t* init_fof_groups()
{
  mlog("Allocating fof_group array with %d elements...", MLOG_MESG, run_globals.NFOFGroupsMax);
  fof_group_t* fof_groups = malloc(sizeof(fof_group_t) * run_globals.NFOFGroupsMax);

  for (int ii = 0; ii < run_globals.NFOFGroupsMax; ii++) {
    fof_groups[ii].FirstHalo = NULL;
    fof_groups[ii].FirstOccupiedHalo = NULL;
    fof_groups[ii].Mvir = 0.0;
  }

  return fof_groups;
}

/**
 * Dump the RequestedForestId lists for all ranks.
 */
// static void dump_forest_list()
// {
//     for (int i_rank=0; i_rank < run_globals.mpi_size; ++i_rank) {
//         if (i_rank == run_globals.mpi_rank) {
//             hid_t fd;
//             const char fname[] = "requested_forest_ids.h5";

//             if (run_globals.mpi_rank == 0)
//                 fd = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//             else
//                 fd = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

//             char name[50];
//             sprintf(name, "core%03d", run_globals.mpi_rank);
//             H5LTmake_dataset_long(fd, name, 1, (hsize_t []){run_globals.NRequestedForests},
//             run_globals.RequestedForestId);

//             H5Fflush(fd, H5P_DEFAULT);
//             H5Fclose(fd);
//         }
//         MPI_Barrier(run_globals.mpi_comm);
//     }
// }

void initialize_halo_storage()
{
  int* n_store_snapshots = &(run_globals.NStoreSnapshots);
  halo_t*** snapshot_halo = &(run_globals.SnapshotHalo);
  fof_group_t*** snapshot_fof_group = &(run_globals.SnapshotFOFGroup);
  int*** snapshot_index_lookup = &(run_globals.SnapshotIndexLookup);
  trees_info_t** snapshot_trees_info = &(run_globals.SnapshotTreesInfo);

  int last_snap = 0;

  mlog("Initializing halo storage arrays...", MLOG_OPEN);

  // Find what the last requested output snapshot is
  for (int ii = 0; ii < run_globals.NOutputSnaps; ii++)
    if (run_globals.ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals.ListOutputSnaps[ii];

  // Allocate an array of last_snap halo array pointers
  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC)
    *n_store_snapshots = last_snap + 1;
  else
    *n_store_snapshots = 1;

  *snapshot_halo = (halo_t**)calloc((size_t)*n_store_snapshots, sizeof(halo_t*));
  *snapshot_fof_group = (fof_group_t**)calloc((size_t)*n_store_snapshots, sizeof(fof_group_t*));
  *snapshot_index_lookup = (int**)calloc((size_t)*n_store_snapshots, sizeof(int*));
  *snapshot_trees_info = (trees_info_t*)calloc((size_t)*n_store_snapshots, sizeof(trees_info_t));

  for (int ii = 0; ii < *n_store_snapshots; ii++) {
    (*snapshot_trees_info)[ii].n_halos = -1;
    (*snapshot_index_lookup)[ii] = NULL;
  }

  // loop through and read all snapshots
  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC) {
    mlog("Preloading input trees and halos...", MLOG_OPEN);
    for (int i_snap = 0; i_snap <= last_snap; i_snap++)
      read_halos(i_snap,
                 &((*snapshot_halo)[i_snap]),
                 &((*snapshot_fof_group)[i_snap]),
                 &((*snapshot_index_lookup)[i_snap]),
                 *snapshot_trees_info);
    mlog("...done", MLOG_CLOSE);
  }

  mlog("...done", MLOG_CLOSE);
}

static void select_forests()
{
  // search the input tree files for all unique forest ids, store them, sort
  // them, and then potentially split them amongst cores
  mlog("Calling select_forests()...", MLOG_MESG | MLOG_TIMERSTART);

  int* rank_n_assigned = NULL;
  long* assigned_ids = NULL;
  int* rank_max_contemp_halo = 0;
  int* rank_max_contemp_fof = 0;
  int n_forests = 0;
  int max_rank_n_forests = 0;

  if (run_globals.mpi_rank == 0) {
    char fname[STRLEN + 34];

    switch (run_globals.params.TreesID) {
      case VELOCIRAPTOR_TREES:
        sprintf(fname, "%s/trees/meraxes_augmented_stats.h5", run_globals.params.SimulationDir);
        break;
      case GBPTREES_TREES:
        sprintf(fname, "%s/trees/forests_info.hdf5", run_globals.params.SimulationDir);
        break;
      default:
        mlog_error("Unrecognised input trees identifier (TreesID).");
        break;
    }

    hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(EXIT_FAILURE);
    }

    // find out how many forests there are
    switch (run_globals.params.TreesID) {
      case GBPTREES_TREES:
        H5LTget_attribute_int(fd, "info", "n_forests", &n_forests);
        break;
      case VELOCIRAPTOR_TREES:
        H5LTget_attribute_int(fd, "forests", "n_forests", &n_forests);
        break;
      default:
        mlog_error("Unrecognised TreesID parameter.");
        break;
    }

    // read in the total halo counts for each forest at the last snapshot
    int last_snap = 0;
    for (int ii = 0; ii < run_globals.NOutputSnaps; ii++) {
      if (run_globals.ListOutputSnaps[ii] > last_snap) {
        last_snap = run_globals.ListOutputSnaps[ii];
      }
    }

    long* forest_ids = (long*)malloc(sizeof(long) * n_forests);
    int* final_counts = (int*)malloc(sizeof(int) * n_forests);
    int* max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
    int* max_contemp_fof = (int*)malloc(sizeof(int) * n_forests);

    {
      int* temp_ids;
      char dset_name[128] = { '\0' };

      switch (run_globals.params.TreesID) {
        case GBPTREES_TREES:
          temp_ids = (int*)malloc(sizeof(int) * n_forests);
          H5LTread_dataset_int(fd, "info/forest_id", temp_ids);
          for (int ii = 0; ii < n_forests; ++ii) {
            forest_ids[ii] = (long)temp_ids[ii];
          }
          free(temp_ids);
          H5LTread_dataset_int(fd, "info/max_contemporaneous_halos", max_contemp_halo);
          H5LTread_dataset_int(fd, "info/max_contemporaneous_fof_groups", max_contemp_fof);
          sprintf(dset_name, "snapshots/snap_%03d", last_snap);
          H5LTread_dataset_int(fd, dset_name, final_counts);
          break;
        case VELOCIRAPTOR_TREES:
          H5LTread_dataset_long(fd, "forests/forest_ids", forest_ids);
          H5LTread_dataset_int(fd, "forests/max_contemporaneous_halos", max_contemp_halo);
          H5LTread_dataset_int(fd, "forests/max_contemporaneous_fof_groups", max_contemp_fof);
          sprintf(dset_name, "snapshots/Snap%03d", last_snap);
          H5LTread_dataset_int(fd, dset_name, final_counts);
          break;
        default:
          mlog_error("Unrecognised TreesID parameter.");
          break;
      }
    }

    // If we have requested forest IDs already (ie. read in from a file) then
    // we will set the final_counts of all forest IDs not in this list to zero.
    // WARNING: This block is not well (if at all!) tested!!!
    // NOTE: This method assumes that forest_id is sorted by value.
    if (run_globals.RequestedForestId != NULL) {
      qsort(run_globals.RequestedForestId, (size_t)run_globals.NRequestedForests, sizeof(long), compare_longs);
      int jj = 0;
      for (int ii = 0; ii < run_globals.NRequestedForests; ++ii) {
        while (forest_ids[jj] < run_globals.RequestedForestId[ii]) {
          final_counts[jj++] = 0;
        }
        assert(forest_ids[jj] == run_globals.RequestedForestId[ii]);
      }
    }

    // We'll use this for a check below
    unsigned long true_total = 0;
    for (int ii = 0; ii < n_forests; ++ii) {
      true_total += final_counts[ii];
    }

    // indirectly sort the final counts and store the sort indices (descending order)
    size_t* sort_ind = calloc(n_forests, sizeof(size_t));
    gsl_sort_int_index(sort_ind, final_counts, 1, n_forests);
    {
      int ii = 0;
      int jj = 0;
      while (ii < jj) {
        int tmp = sort_ind[ii];
        sort_ind[ii] = sort_ind[jj];
        sort_ind[jj] = tmp;
        ii++;
        jj--;
      }
    }

    max_rank_n_forests = (int)((float)n_forests * 0.75);
    rank_n_assigned = calloc(run_globals.mpi_size, sizeof(int));
    if (rank_n_assigned == NULL) {
      mlog_error("Failed to allocate rank_n_assigned array.");
      ABORT(EXIT_FAILURE);
    }
    assigned_ids = calloc(max_rank_n_forests * run_globals.mpi_size, sizeof(long));
    if (assigned_ids == NULL) {
      mlog_error("Failed to allocate assigned_ids array.");
      ABORT(EXIT_FAILURE);
    }
    rank_max_contemp_halo = calloc(run_globals.mpi_size, sizeof(int));
    if (rank_max_contemp_halo == NULL) {
      mlog_error("Failed to allocate rank_max_contemp_halo array.");
      ABORT(EXIT_FAILURE);
    }
    rank_max_contemp_fof = calloc(run_globals.mpi_size, sizeof(int));
    if (rank_max_contemp_fof == NULL) {
      mlog_error("Failed to allocate rank_max_contemp_fof array.");
      ABORT(EXIT_FAILURE);
    }

    int* rank_counts = calloc(run_globals.mpi_size, sizeof(int));
    int* rank_argsort_ind = malloc(run_globals.mpi_size * sizeof(int));
    for (int ii = 0; ii < run_globals.mpi_size; ++ii) {
      rank_argsort_ind[ii] = ii;
    }

    int* snap_counts = NULL;
    snap_counts = calloc(n_forests, sizeof(int));
    assert(snap_counts != NULL);

    // Save old hdf5 error handler and turn off error handling
    herr_t (*old_func)(long, void*) = NULL;
    void* old_client_data = NULL;
    hid_t estack_id = 0;
    H5Eget_auto(estack_id, &old_func, &old_client_data);
    H5Eset_auto(estack_id, NULL, NULL);

    for (int snap = 0; snap < last_snap + 1; ++snap) {
      char dset_name[128] = { '\0' };
      switch (run_globals.params.TreesID) {
        case GBPTREES_TREES:
          sprintf(dset_name, "snapshots/snap_%03d", snap);
          break;
        case VELOCIRAPTOR_TREES:
          sprintf(dset_name, "snapshots/Snap%03d", snap);
          break;
        default:
          mlog_error("Unrecognised input trees identifier (TreesID).");
          break;
      }

      herr_t status = H5LTread_dataset_int(fd, dset_name, snap_counts);
      if (status < 0) {
        continue;
      }

      // loop through non-zero counts
      for (int ii = 0; ii < n_forests; ++ii) {
        int jj = sort_ind[ii];

        if ((snap_counts[jj] > 0) && (final_counts[jj] > 0)) {

          // first appearance!
          int rank = rank_argsort_ind[0];
          rank_counts[rank] += final_counts[jj];
          assigned_ids[rank * max_rank_n_forests + rank_n_assigned[rank]] = forest_ids[jj];
          rank_n_assigned[rank]++;
          if (rank_n_assigned[rank] >= max_rank_n_forests) {
            mlog_error("Forest load-imbalance is above currently allowed threshold.");
            ABORT(EXIT_FAILURE);
          }
          final_counts[jj] = 0;

          rank_max_contemp_halo[rank] += max_contemp_halo[jj];
          rank_max_contemp_fof[rank] += max_contemp_fof[jj];

          // keep the rank_argsort_inds correct
          for (int kk = 1; kk < run_globals.mpi_size; ++kk) {
            if (rank_counts[rank_argsort_ind[kk]] < rank_counts[rank_argsort_ind[kk - 1]]) {
              int tmp = rank_argsort_ind[kk];
              rank_argsort_ind[kk] = rank_argsort_ind[kk - 1];
              rank_argsort_ind[kk - 1] = tmp;
            }
          }
        }
      }
    }

    // Restore previous hdf5 error handler
    H5Eset_auto(estack_id, old_func, old_client_data);

    // Do a quick sanity check
    unsigned long total = 0;
    for (int ii = 0; ii < run_globals.mpi_size; ++ii) {
      total += rank_counts[ii];
    }
    assert(total == true_total);

    free(snap_counts);
    free(rank_argsort_ind);
    free(rank_counts);
    free(sort_ind);
    free(max_contemp_fof);
    free(max_contemp_halo);
    free(final_counts);
    free(forest_ids);

    H5Fclose(fd);
  } // mpi_rank = 0

  // let all ranks know what their forest ID lists are
  MPI_Scatter(rank_n_assigned, 1, MPI_INT, &run_globals.NRequestedForests, 1, MPI_INT, 0, run_globals.mpi_comm);
  if (run_globals.RequestedForestId != NULL)
    run_globals.RequestedForestId =
      realloc(run_globals.RequestedForestId, run_globals.NRequestedForests * sizeof(long));
  else
    run_globals.RequestedForestId = (long*)malloc(sizeof(long) * run_globals.NRequestedForests);

  // NB: displs only valid on rank 0
  int* displs = calloc(run_globals.mpi_size, sizeof(int));
  for (int ii = 0; ii < run_globals.mpi_size; ++ii) {
    displs[ii] = ii * max_rank_n_forests;
  }
  MPI_Scatterv(assigned_ids,
               rank_n_assigned,
               displs,
               MPI_LONG,
               run_globals.RequestedForestId,
               run_globals.NRequestedForests,
               MPI_LONG,
               0,
               run_globals.mpi_comm);
  free(displs);

  // let all ranks know what their max allocation counts are
  MPI_Scatter(rank_max_contemp_halo, 1, MPI_INT, &run_globals.NHalosMax, 1, MPI_INT, 0, run_globals.mpi_comm);
  MPI_Scatter(rank_max_contemp_fof, 1, MPI_INT, &run_globals.NFOFGroupsMax, 1, MPI_INT, 0, run_globals.mpi_comm);

  if (run_globals.mpi_rank == 0) {
    free(rank_max_contemp_fof);
    free(rank_max_contemp_halo);
    free(assigned_ids);
    free(rank_n_assigned);
  }

  // sort the requested forest ids so that they can be bsearch'd later
  qsort(run_globals.RequestedForestId, (size_t)run_globals.NRequestedForests, sizeof(long), compare_longs);

  mlog("...done.", MLOG_MESG | MLOG_TIMERSTOP);
}

trees_info_t read_halos(const int snapshot,
                        halo_t** halos,
                        fof_group_t** fof_groups,
                        int** index_lookup,
                        trees_info_t* snapshot_trees_info)
{
  // if we are doing multiple runs and have already read in this snapshot then we don't need to do read anything else
  if ((run_globals.params.FlagInteractive || run_globals.params.FlagMCMC) &&
      (snapshot_trees_info[snapshot].n_halos != -1)) {
    mlog("Snapshot %d has already been read in... (n_halos = %d)",
         MLOG_MESG,
         snapshot,
         snapshot_trees_info[snapshot].n_halos);
    return snapshot_trees_info[snapshot];
  }

  mlog("Reading snapshot %d (z = %.2f) trees and halos...",
       MLOG_OPEN | MLOG_TIMERSTART,
       snapshot,
       run_globals.ZZ[snapshot]);

  // Read mass ratio modifiers and baryon fraction modifiers if required
  if (run_globals.RequestedMassRatioModifier == 1)
    read_mass_ratio_modifiers(snapshot);

  if (run_globals.RequestedBaryonFracModifier == 1)
    read_baryon_frac_modifiers(snapshot);

  // read in the tree information for this snapshot
  mlog("Reading trees info...", MLOG_MESG | MLOG_TIMERSTART);

  trees_info_t trees_info = { 0 };
  switch (run_globals.params.TreesID) {
    case VELOCIRAPTOR_TREES:
      trees_info = read_trees_info__velociraptor(snapshot);
      break;
    case GBPTREES_TREES:
      trees_info = read_trees_info__gbptrees(snapshot);
      break;
    default:
      mlog_error("Unrecognised input trees identifier (TreesID).");
      break;
  }

  mlog("... done.", MLOG_CONT | MLOG_TIMERSTOP);

  int n_halos = trees_info.n_halos;
  int n_fof_groups = trees_info.n_fof_groups;

  // TODO: Here the original code returned if there were no halos at this
  // snapshot and we were in interactive / MCMC mode...  Not sure why the
  // caveat though...
  if (n_halos < 1) {
    mlog("No halos in this file... skipping...", MLOG_CLOSE | MLOG_TIMERSTOP);
    return trees_info;
  }

  if ((*halos) == NULL) {
    // if required, select forests and calculate the maximum number of halos and
    // fof groups
    if ((run_globals.NRequestedForests > -1) || (run_globals.mpi_size > 1)) {
      if (run_globals.SelectForestsSwitch == true) {
        select_forests();

        // Toggle the SelectForestsSwitch so that we know we don't need to call
        // this again
        run_globals.SelectForestsSwitch = false;
      }
      *index_lookup = malloc(sizeof(int) * run_globals.NHalosMax);
      for (int ii = 0; ii < run_globals.NHalosMax; ii++)
        (*index_lookup)[ii] = -1;
    } else {
      run_globals.NHalosMax = trees_info.n_halos_max;
      run_globals.NFOFGroupsMax = trees_info.n_fof_groups_max;
    }

    mlog("Allocating halo array with %d elements...", MLOG_MESG, run_globals.NHalosMax);
    *halos = malloc(sizeof(halo_t) * run_globals.NHalosMax);
  }

  // Allocate the fof_group array if necessary
  if (*fof_groups == NULL)
    *fof_groups = init_fof_groups();

  // Now actually read in the trees!
  switch (run_globals.params.TreesID) {
    case VELOCIRAPTOR_TREES:
      read_trees__velociraptor(snapshot, *halos, &n_halos, *fof_groups, &n_fof_groups, *index_lookup);
      break;

    case GBPTREES_TREES: {
      int n_halos_kept = 0;
      int n_fof_groups_kept = 0;

      read_trees__gbptrees(snapshot,
                           *halos,
                           n_halos,
                           *fof_groups,
                           run_globals.NRequestedForests,
                           &n_halos_kept,
                           &n_fof_groups_kept,
                           *index_lookup);

      n_halos = n_halos_kept;
      n_fof_groups = n_fof_groups_kept;
    }

    break;

    default:
      mlog_error("Unrecognised input trees identifier (TreesID).");
      break;
  }

  // if subsampling the trees, then update the trees_info to reflect what we now have
  if (run_globals.NRequestedForests > -1) {
    trees_info.n_halos = n_halos;
    trees_info.n_fof_groups = n_fof_groups;
  }

  // if we are doing multiple runs then resize the arrays to save space and store the trees_info
  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC) {
    // Ok - what follows here is hacky as hell.  By calling realloc on these
    // arrays, there is a good chance that the actual array will be moved and
    // there is no way to prevent this.  A side effect will be that all of the
    // pointers that refer to any of these arrays will be broken.  The simplest
    // way to deal with this is to calculate the array offsets which each
    // pointer refers to, and then reset all of the pointers in the realloc'd
    // arrays.
    //
    // In future, this should definitely be changed. All pointers to these
    // array members should be made integer offsets instead.  That will require
    // trawling through the code and making the relevant updates in a number of
    // places (so we'll leave that till later!).

    mlog("Calculating pointer offsets...", MLOG_MESG);

    size_t* halo_FOFGroup_os = malloc(sizeof(size_t) * n_halos);
    size_t* halo_NextHaloInFOFGroup_os = malloc(sizeof(size_t) * n_halos);
    size_t* fof_FirstHalo_os = malloc(sizeof(size_t) * n_fof_groups);

    for (int ii = 0; ii < n_halos; ii++) {
      halo_FOFGroup_os[ii] = (size_t)((*halos)[ii].FOFGroup - (*fof_groups));
      if ((*halos)[ii].NextHaloInFOFGroup != NULL)
        halo_NextHaloInFOFGroup_os[ii] = (size_t)((*halos)[ii].NextHaloInFOFGroup - (*halos));
      else
        halo_NextHaloInFOFGroup_os[ii] = (size_t)-1;
    }
    for (int ii = 0; ii < n_fof_groups; ii++)
      fof_FirstHalo_os[ii] = (size_t)((*fof_groups)[ii].FirstHalo - (*halos));

    mlog("Reallocing halo storage arrays...", MLOG_OPEN);

    *halos = (halo_t*)realloc(*halos, sizeof(halo_t) * n_halos);
    *fof_groups = (fof_group_t*)realloc(*fof_groups, sizeof(fof_group_t) * n_fof_groups);

    // Only realloc `index_lookup` if we are actually using it (`n_proc` > 1
    // or we are subsampling the trees).
    if (*index_lookup)
      *index_lookup = (int*)realloc(*index_lookup, sizeof(int) * n_halos);
    update_pointers_from_offsets(
      n_fof_groups, *fof_groups, fof_FirstHalo_os, n_halos, *halos, halo_FOFGroup_os, halo_NextHaloInFOFGroup_os);

    // save the trees_info for this snapshot as well...
    snapshot_trees_info[snapshot] = trees_info;

    mlog(" ...done (resized to %d halos on rank=0)", MLOG_CLOSE, n_halos);

    free(fof_FirstHalo_os);
    free(halo_NextHaloInFOFGroup_os);
    free(halo_FOFGroup_os);
  }

  MPI_Allreduce(MPI_IN_PLACE, &n_halos, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &n_fof_groups, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
  mlog("Read %d halos in %d fof_groups.", MLOG_MESG, n_halos, n_fof_groups);

  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  return trees_info;
}

void free_halo_storage()
{
  int n_store_snapshots = run_globals.NStoreSnapshots;
  halo_t** snapshot_halo = run_globals.SnapshotHalo;
  fof_group_t** snapshot_fof_group = run_globals.SnapshotFOFGroup;
  int** snapshot_index_lookup = run_globals.SnapshotIndexLookup;
  trees_info_t* snapshot_trees_info = run_globals.SnapshotTreesInfo;

  // Free all of the remaining allocated galaxies, halos and fof groups
  mlog("Freeing FOF groups and halos...", MLOG_MESG);
  for (int ii = 0; ii < n_store_snapshots; ii++) {
    free(snapshot_halo[ii]);
    free(snapshot_fof_group[ii]);
    free(snapshot_index_lookup[ii]);
  }
  free(snapshot_halo);
  free(snapshot_fof_group);
  free(snapshot_index_lookup);
  free(snapshot_trees_info);
}
