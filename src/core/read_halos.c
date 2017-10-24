#include "meraxes.h"
#include <assert.h>
#include <gsl/gsl_sort_int.h>
#include <hdf5_hl.h>

static void inline update_pointers_from_offsets(
    int n_fof_groups_kept,
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

static void reorder_forest_array(int* arr, const size_t* sort_ind, int n_forests, int* temp)
{
    memcpy(temp, arr, sizeof(int) * n_forests);
    for (int ii = 0, jj = n_forests - 1; ii < n_forests; ii++, jj--)
        arr[ii] = temp[sort_ind[jj]];
}

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
            read_halos(i_snap, &((*snapshot_halo)[i_snap]), &((*snapshot_fof_group)[i_snap]), &((*snapshot_index_lookup)[i_snap]), *snapshot_trees_info);
        mlog("...done", MLOG_CLOSE);
    }

    mlog("...done", MLOG_CLOSE);
}

static void select_forests()
{
    // search the input tree files for all unique forest ids, store them, sort
    // them, and then potentially split them amoungst cores
    mlog("Calling select_forests()...", MLOG_MESG);

    // if this is the master rank then read in the forest info
    int *max_contemp_halo = NULL, *max_contemp_fof = NULL, *forest_ids = NULL, *n_halos = NULL;
    int n_forests = 0;
    if (run_globals.mpi_rank == 0) {
        char fname[STRLEN];
        char grp_name[30];

        switch (run_globals.params.TreesID) {
            case VELOCIRAPTOR_TREES:
                sprintf(fname, "%s/trees/meraxes_augmented_stats.h5", run_globals.params.SimulationDir);
                strcpy(grp_name, "forests\0");
                break;
            case GBPTREES_TREES:
                sprintf(fname, "%s/trees/forests_info.hdf5", run_globals.params.SimulationDir);
                strcpy(grp_name, "info\0");
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
        H5LTget_attribute_int(fd, grp_name, "n_forests", &n_forests);

        // allocate the arrays
        max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
        max_contemp_fof = (int*)malloc(sizeof(int) * n_forests);
        forest_ids = (int*)malloc(sizeof(int) * n_forests);
        n_halos = (int*)malloc(sizeof(int) * n_forests);

        // read in the max number of contemporaneous halos and groups and the forest ids
        hid_t grp = H5Gopen2(fd, grp_name, H5P_DEFAULT);
        H5LTread_dataset_int(grp, "max_contemporaneous_halos", max_contemp_halo);
        H5LTread_dataset_int(grp, "max_contemporaneous_fof_groups", max_contemp_fof);
        H5LTread_dataset_int(grp, "forest_ids", forest_ids);
        H5LTread_dataset_int(grp, "n_halos", n_halos);
        H5Gclose(grp);

        // close the file
        H5Fclose(fd);
    }

    // broadcast the forest info
    MPI_Bcast(&n_forests, 1, MPI_INT, 0, run_globals.mpi_comm);
    if (run_globals.mpi_rank > 0) {
        max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
        max_contemp_fof = (int*)malloc(sizeof(int) * n_forests);
        forest_ids = (int*)malloc(sizeof(int) * n_forests);
        n_halos = (int*)malloc(sizeof(int) * n_forests);
    }
    MPI_Bcast(max_contemp_halo, n_forests, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(max_contemp_fof, n_forests, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(forest_ids, n_forests, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(n_halos, n_forests, MPI_INT, 0, run_globals.mpi_comm);

    // if we have read in a list of requested forest IDs then use these to create
    // an array of indices pointing to the elements we want
    int* requested_ind;
    int n_halos_tot = 0;
    if (run_globals.RequestedForestId != NULL) {
        requested_ind = malloc(sizeof(int) * run_globals.NRequestedForests);
        for (int i_forest = 0, i_req = 0; (i_forest < n_forests) && (i_req < run_globals.NRequestedForests); i_forest++)
            if (forest_ids[i_forest] == run_globals.RequestedForestId[i_req]) {
                requested_ind[i_req] = i_forest;
                n_halos_tot += n_halos[i_forest];
                i_req++;
            }
    }
    else {
        // if we haven't asked for any specific forest IDs then just fill the
        // requested ind array sequentially
        requested_ind = (int*)malloc(sizeof(int) * n_forests);
        for (int i_req = 0; i_req < n_forests; i_req++) {
            n_halos_tot += n_halos[i_req];
            requested_ind[i_req] = i_req;
        }
    }

    // sort the forests by the number of halos in each one
    size_t* sort_ind = malloc(sizeof(size_t) * n_forests);
    gsl_sort_int_index(sort_ind, n_halos, 1, (const size_t)n_forests);
    int* temp = malloc(sizeof(int) * n_forests);

    reorder_forest_array(forest_ids, sort_ind, n_forests, temp);
    reorder_forest_array(max_contemp_halo, sort_ind, n_forests, temp);
    reorder_forest_array(max_contemp_fof, sort_ind, n_forests, temp);
    reorder_forest_array(n_halos, sort_ind, n_forests, temp);

    // also rearrange the sampled forest indices if need be
    if (run_globals.RequestedForestId != NULL) {
        int* rank = (int*)malloc(sizeof(int) * n_forests);
        for (int ii = 0; ii < n_forests; ii++)
            rank[sort_ind[ii]] = ii;

        for (int ii = 0; ii < run_globals.NRequestedForests; ii++)
            requested_ind[ii] = n_forests - 1 - rank[requested_ind[ii]];

        free(rank);

        // from this point on we are only concerned with what is on / going to each
        // core, therefore we can reset n_forests appropriately
        n_forests = run_globals.NRequestedForests;
    }

    free(temp);
    free(sort_ind);

    // if we are running the code with more than one core, let's subselect the
    // forests on each one
    if (run_globals.mpi_size > 1) {
        if (n_forests < run_globals.mpi_size) {
            mlog_error(
                "There are fewer processors than there are forests to be processed "
                "(%d).  Try again with fewer cores!",
                n_forests);
            ABORT(EXIT_FAILURE);
        }

        // TODO: Load balancing still seems pretty bad.
        // This should be looked into to see if there are any improvements that can
        // be made.  This would hopefully improve runtimes...

        // malloc the arrays we need for keeping track of the load balancing
        int* rank_first_forest = (int*)malloc(sizeof(int) * run_globals.mpi_size);
        int* rank_last_forest = (int*)malloc(sizeof(int) * run_globals.mpi_size);
        int* rank_n_forests = (int*)malloc(sizeof(int) * run_globals.mpi_size);
        int* rank_n_halos = (int*)malloc(sizeof(int) * run_globals.mpi_size);

        // initialise
        for (int ii = 0; ii < run_globals.mpi_size; ii++) {
            rank_first_forest[ii] = 0;
            rank_last_forest[ii] = 0;
            rank_n_forests[ii] = 0;
            rank_n_halos[ii] = 0;
        }

        // loop through each rank
        int i_forest = 0;
        for (int i_rank = 0, n_halos_used = 0, n_halos_target = 0;
             i_rank < run_globals.mpi_size; i_rank++, i_forest++) {
            // start with the next forest (or first forest if this is the first rank)
            rank_first_forest[i_rank] = i_forest;
            rank_last_forest[i_rank] = i_forest;

            // if we still have forests left to check from the total list of forests
            if (i_forest < n_forests) {
                // add this forest's halos to this rank's total
                rank_n_halos[i_rank] += n_halos[requested_ind[i_forest]];

                // adjust our target to smooth out variations as much as possible
                n_halos_target = (n_halos_tot - n_halos_used) / (run_globals.mpi_size - i_rank);

                // increment the counter of the number of forests on this rank
                rank_n_forests[i_rank]++;

                // keep adding forests until we have reached (or exceeded) the current
                // target
                while ((rank_n_halos[i_rank] < n_halos_target) && (i_forest < (n_forests - 1))) {
                    i_forest++;
                    rank_n_forests[i_rank]++;
                    rank_last_forest[i_rank] = i_forest;
                    rank_n_halos[i_rank] += n_halos[requested_ind[i_forest]];
                }

                // updated the total number of used halos
                n_halos_used += rank_n_halos[i_rank];
            }
        }

        // add any uncounted forests to the last process
        // n.b. i_forest = value from last for loop
        for (i_forest++; i_forest < n_forests; i_forest++) {
            rank_last_forest[run_globals.mpi_size - 1] = i_forest;
            rank_n_halos[run_globals.mpi_size - 1] += n_halos[requested_ind[i_forest]];
            rank_n_forests[run_globals.mpi_size - 1]++;
        }

        assert(rank_n_forests[run_globals.mpi_rank] > 0);

        // create our list of forest_ids for this rank
        run_globals.NRequestedForests = rank_n_forests[run_globals.mpi_rank];
        if (run_globals.RequestedForestId != NULL)
            run_globals.RequestedForestId = realloc(run_globals.RequestedForestId,
                run_globals.NRequestedForests * sizeof(int));
        else
            run_globals.RequestedForestId = (int*)malloc(sizeof(int) * run_globals.NRequestedForests);

        for (int ii = rank_first_forest[run_globals.mpi_rank], jj = 0;
             ii <= rank_last_forest[run_globals.mpi_rank];
             ii++, jj++)
            run_globals.RequestedForestId[jj] = forest_ids[requested_ind[ii]];

        // free the arrays
        free(rank_n_halos);
        free(rank_n_forests);
        free(rank_last_forest);
        free(rank_first_forest);
    }

    assert(run_globals.RequestedForestId != NULL);

    // loop through and tot up the max number of halos and fof_groups we will need
    // to allocate
    int max_halos = 0;
    int max_fof_groups = 0;
    for (int i_forest = 0, i_req = 0;
         (i_forest < n_forests) && (i_req < run_globals.NRequestedForests);
         i_forest++)
        if (forest_ids[requested_ind[i_forest]] == run_globals.RequestedForestId[i_req]) {
            max_halos += max_contemp_halo[requested_ind[i_forest]];
            max_fof_groups += max_contemp_fof[requested_ind[i_forest]];
            i_req++;
        }

    // store the maximum number of halos and fof groups needed at any one snapshot
    run_globals.NHalosMax = max_halos;
    run_globals.NFOFGroupsMax = max_fof_groups;

    // sort the requested forest ids so that they can be bsearch'd later
    qsort(run_globals.RequestedForestId, (size_t)run_globals.NRequestedForests,
        sizeof(int), compare_ints);

    free(requested_ind);
    free(max_contemp_halo);
    free(max_contemp_fof);
    free(forest_ids);
    free(n_halos);
}

trees_info_t read_halos(const int snapshot, halo_t** halos, fof_group_t** fof_groups,
    int** index_lookup, trees_info_t* snapshot_trees_info)
{
    // if we are doing multiple runs and have already read in this snapshot then we don't need to do read anything else
    if ((run_globals.params.FlagInteractive || run_globals.params.FlagMCMC) && (snapshot_trees_info[snapshot].n_halos != -1)) {
        mlog("Snapshot %d has already been read in... (n_halos = %d)", MLOG_MESG, snapshot, snapshot_trees_info[snapshot].n_halos);
        return snapshot_trees_info[snapshot];
    }

    mlog("Reading snapshot %d (z = %.2f) trees and halos...",
            MLOG_OPEN | MLOG_TIMERSTART, snapshot, run_globals.ZZ[snapshot]);

    // Read mass ratio modifiers and baryon fraction modifiers if required
    if (run_globals.RequestedMassRatioModifier == 1)
        read_mass_ratio_modifiers(snapshot);

    if (run_globals.RequestedBaryonFracModifier == 1)
        read_baryon_frac_modifiers(snapshot);

    // read in the tree information for this snapshot
    trees_info_t trees_info;
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
        }
        else {
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

        case GBPTREES_TREES:
            {
                int n_halos_kept = 0;
                int n_fof_groups_kept = 0;

                read_trees__gbptrees(snapshot, *halos, n_halos,
                        *fof_groups, n_fof_groups, run_globals.NRequestedForests,
                        &n_halos_kept, &n_fof_groups_kept, *index_lookup);

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
        update_pointers_from_offsets(n_fof_groups, *fof_groups, fof_FirstHalo_os,
                n_halos, *halos, halo_FOFGroup_os, halo_NextHaloInFOFGroup_os);

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
