#include "meraxes.h"
#include <assert.h>
#include <gsl/gsl_sort_int.h>
#include <hdf5.h>
#include <hdf5_hl.h>

static trees_info_t read_trees_info(const int snapshot)
{
    // TODO: This is wasteful and should probably only ever be done once and stored in run_globals.
    char fname[STRLEN];
    sprintf(fname, "%s/trees/meraxes_augmented_info.h5", run_globals.params.SimulationDir);

    hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) {
        mlog("Failed to open file %s", MLOG_MESG, fname);
        ABORT(EXIT_FAILURE);
    }

    int n_snaps = 0;
    trees_info_t trees_info;
    H5LTget_attribute_int(fd, "/", "n_snaps", &n_snaps);
    H5LTget_attribute_int(fd, "/", "n_halos_max", &(trees_info.n_halos_max));
    H5LTget_attribute_int(fd, "/", "n_fof_groups_max", &(trees_info.n_fof_groups_max));

    int* buffer = malloc(sizeof(int) * n_snaps);
    H5LTread_dataset_int(fd, "n_halos", buffer);
    trees_info.n_halos = buffer[snapshot];
    H5LTread_dataset_int(fd, "n_fof_groups", buffer);
    trees_info.n_fof_groups = buffer[snapshot];
    free(buffer);

    H5Fclose(fd);

    return trees_info;
}

static void reorder_forest_array(int* arr, size_t* sort_ind, int n_forests, int* temp)
{
    memcpy(temp, arr, sizeof(int) * n_forests);
    for (int ii = 0, jj = n_forests - 1; ii < n_forests; ii++, jj--)
        arr[ii] = temp[sort_ind[jj]];
}

static void select_forests()
{
    // search the input tree files for all unique forest ids, store them, sort
    // them, and then potentially split them amoungst cores
    mlog("Calling select_forests()...", MLOG_MESG);

    // are we sampling the forests or just dividing them amoungst cores?
    bool sample_forests = false;
    if (run_globals.RequestedForestId != NULL)
        sample_forests = true;

    // if this is the master rank then read in the forest info
    int *max_contemp_halo, *max_contemp_fof, *forest_id, *n_halos, n_forests;
    if (run_globals.mpi_rank == 0) {
        char fname[STRLEN];
        sprintf(fname, "%s/trees/meraxes_augmented_info.h5", run_globals.params.SimulationDir);

        hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fd < 0) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(EXIT_FAILURE);
        }

        // find out how many forests there are
        H5LTget_attribute_int(fd, "forests", "n_forests", &n_forests);

        // allocate the arrays
        max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
        max_contemp_fof = (int*)malloc(sizeof(int) * n_forests);
        forest_id = (int*)malloc(sizeof(int) * n_forests);
        n_halos = (int*)malloc(sizeof(int) * n_forests);

        // read in the max number of contemporaneous halos and groups and the forest ids
        H5LTread_dataset_int(fd, "forests/max_contemporaneous_halos", max_contemp_halo);
        H5LTread_dataset_int(fd, "forests/max_contemporaneous_fof_groups", max_contemp_fof);
        H5LTread_dataset_int(fd, "forests/forest_id", forest_id);
        H5LTread_dataset_int(fd, "forests/n_halos", n_halos);

        // close the file
        H5Fclose(fd);
    }

    // broadcast the forest info
    MPI_Bcast(&n_forests, 1, MPI_INT, 0, run_globals.mpi_comm);
    if (run_globals.mpi_rank > 0) {
        max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
        max_contemp_fof = (int*)malloc(sizeof(int) * n_forests);
        forest_id = (int*)malloc(sizeof(int) * n_forests);
        n_halos = (int*)malloc(sizeof(int) * n_forests);
    }
    MPI_Bcast(max_contemp_halo, n_forests, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(max_contemp_fof, n_forests, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(forest_id, n_forests, MPI_INT, 0, run_globals.mpi_comm);
    MPI_Bcast(n_halos, n_forests, MPI_INT, 0, run_globals.mpi_comm);

    // if we have read in a list of requested forest IDs then use these to create
    // an array of indices pointing to the elements we want
    int* requested_ind;
    int n_halos_tot = 0;
    if (sample_forests) {
        requested_ind = malloc(sizeof(int) * run_globals.NRequestedForests);
        for (int i_forest = 0, i_req = 0; (i_forest < n_forests) && (i_req < run_globals.NRequestedForests); i_forest++)
            if (forest_id[i_forest] == run_globals.RequestedForestId[i_req]) {
                requested_ind[i_req] = i_forest;
                n_halos_tot += n_halos[i_forest];
                i_req++;
            }
    } else {
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
    gsl_sort_int_index(sort_ind, n_halos, 1, n_forests);
    int* temp = malloc(sizeof(int) * n_forests);

    reorder_forest_array(forest_id, sort_ind, n_forests, temp);
    reorder_forest_array(max_contemp_halo, sort_ind, n_forests, temp);
    reorder_forest_array(max_contemp_fof, sort_ind, n_forests, temp);
    reorder_forest_array(n_halos, sort_ind, n_forests, temp);

    // also rearrange the sampled forest indices if need be
    if (sample_forests) {
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
            run_globals.RequestedForestId[jj] = forest_id[requested_ind[ii]];

        // free the arrays
        free(rank_n_halos);
        free(rank_n_forests);
        free(rank_last_forest);
        free(rank_first_forest);
    }

    // loop through and tot up the max number of halos and fof_groups we will need
    // to allocate
    int max_halos = 0;
    int max_fof_groups = 0;
    for (int i_forest = 0, i_req = 0;
         (i_forest < n_forests) && (i_req < run_globals.NRequestedForests);
         i_forest++)
        if (forest_id[requested_ind[i_forest]] == run_globals.RequestedForestId[i_req]) {
            max_halos += max_contemp_halo[requested_ind[i_forest]];
            max_fof_groups += max_contemp_fof[requested_ind[i_forest]];
            i_req++;
        }

    // store the maximum number of halos and fof groups needed at any one snapshot
    run_globals.NHalosMax = max_halos;
    run_globals.NFOFGroupsMax = max_fof_groups;

    // sort the requested forest ids so that they can be bsearch'd later
    qsort(run_globals.RequestedForestId, run_globals.NRequestedForests,
        sizeof(int), compare_ints);

    free(requested_ind);
    free(max_contemp_halo);
    free(max_contemp_fof);
    free(forest_id);
    free(n_halos);
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

static void read_velociraptor_trees(int snapshot, halo_t* halos, int n_halos, fof_group_t* fof_groups, int n_fof_groups, int* index_lookup)
{
}

trees_info_t read_halos(const int snapshot, halo_t** halos, fof_group_t** fof_groups,
    int** index_lookup)
{
    mlog("Reading snapshot %d (z = %.2f) trees and halos...",
        MLOG_OPEN | MLOG_TIMERSTART, snapshot, run_globals.ZZ[snapshot]);

    // Read mass ratio modifiers and baryon fraction modifiers if required
    if (run_globals.RequestedMassRatioModifier == 1)
        read_mass_ratio_modifiers(snapshot);

    if (run_globals.RequestedBaryonFracModifier == 1)
        read_baryon_frac_modifiers(snapshot);

    // open the file and read in the tree information for this snapshot
    trees_info_t trees_info;
    if (run_globals.mpi_rank == 0)
        trees_info = read_trees_info(snapshot);

    // if necessary, broadcast the snapshot info
    MPI_Bcast(&trees_info, sizeof(trees_info_t), MPI_BYTE, 0, run_globals.mpi_comm);

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
    read_velociraptor_trees(snapshot, *halos, n_halos, *fof_groups, n_fof_groups, *index_lookup);

    return trees_info;
}
