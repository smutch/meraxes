#include "meraxes.h"
#include <gsl/gsl_sort_int.h>
#include <hdf5.h>
#include <hdf5_hl.h>

static trees_info_t read_trees_info(const int snapshot)
{
    char fname[STRLEN];
    sprintf(fname, "%s/trees/meraxes_augmented_info.h5", run_globals.params.SimulationDir);

    hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) {
        mlog("Failed to open file %s", MLOG_MESG, fname);
        ABORT(EXIT_FAILURE);
    }

    char grp_name[STRLEN];
    sprintf(grp_name, "Snap_%03d", snapshot);

    trees_info_t trees_info;
    H5LTget_attribute_int(fd, grp_name, "n_halos", &(trees_info.n_halos));
    H5LTget_attribute_int(fd, grp_name, "n_halos_max", &(trees_info.n_halos_max));
    H5LTget_attribute_int(fd, grp_name, "n_fof_groups", &(trees_info.n_fof_groups));
    H5LTget_attribute_int(fd, grp_name, "n_fof_groups_max", &(trees_info.n_fof_groups_max));

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

    // TODO: Cont here...

    free(requested_ind);
}

trees_info_t read_halos(const int snapshot, halo_t** halos, int** index_lookup)
{

    mlog("Reading snapshot %d (z = %.2f) trees and halos...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, run_globals.ZZ[snapshot]);

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
    // snapshot and we were in interactive / MCMC mode...  Not sure why though...

    if ((*halos) == NULL) {
        // if required, select forests and calculate the maximum number of halos and fof groups
        if ((run_globals.NRequestedForests > -1) || (run_globals.mpi_size > 1)) {
            if (run_globals.SelectForestsSwitch == true) {
                select_forests();

                // Toggle the SelectForestsSwitch so that we know we don't need to call this again
                run_globals.SelectForestsSwitch = false;
            }
            *index_lookup = malloc(sizeof(int) * run_globals.NHalosMax);
        } else {
            run_globals.NHalosMax = trees_info.n_halos_max;
            run_globals.NFOFGroupsMax = trees_info.n_fof_groups_max;
        }

        mlog("Allocating halo array with %d elements...", MLOG_MESG, run_globals.NHalosMax);
        *halos = malloc(sizeof(halo_t) * run_globals.NHalosMax);
    }

    return trees_info;
}
