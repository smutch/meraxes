#include "meraxes.h"
#include <assert.h>
#include <gsl/gsl_sort_int.h>
#include <hdf5_hl.h>
#include <math.h>

trees_info_t read_trees_info__velociraptor(const int snapshot)
{
    trees_info_t trees_info;

    if (run_globals.mpi_rank == 0) {
        // TODO: This is wasteful and should probably only ever be done once and stored in run_globals.
        char fname[STRLEN];
        sprintf(fname, "%s/trees/meraxes_augmented_stats.h5", run_globals.params.SimulationDir);

        hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fd < 0) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(EXIT_FAILURE);
        }

        int n_snaps = 0;
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
    }

    // broadcast the snapshot info
    MPI_Bcast(&trees_info, sizeof(trees_info_t), MPI_BYTE, 0, run_globals.mpi_comm);

    return trees_info;
}

static int id_to_ind(long id)
{
    return (int)((id % (uint64_t)1e12) - 1);
}

static int id_to_snap(long id)
{
    return id / 1e12l;
}

static void inline convert_input_virial_props(double* Mvir, double* Rvir, double* Vvir,
    double* FOFMvirModifier, const int len,
    const int snapshot, const bool fof_flag)
{
    // Update the virial properties for subhalos
    if (*Mvir == -1) {
        assert(len > 0);
        *Mvir = calculate_Mvir(*Mvir, len);
    }
    else {
        if (fof_flag && (run_globals.RequestedMassRatioModifier == 1)) {
            // Modifier the FoF mass and update the virial radius
            assert(FOFMvirModifier != NULL);
            *FOFMvirModifier = interpolate_modifier(run_globals.mass_ratio_modifier, log10(*Mvir / run_globals.params.Hubble_h) + 10.0);
            *Mvir *= *FOFMvirModifier;
        }
    }

    if (*Rvir == -1)
        *Rvir = calculate_Rvir(*Mvir, snapshot);

    if (*Vvir == -1)
        *Vvir = calculate_Vvir(*Mvir, *Rvir);
}

#define READ_TREE_ENTRY_PROP(name, type, h5type)             \
    {                                                        \
        H5LTread_dataset(snap_group, #name, h5type, buffer); \
        for (int ii = 0; ii < n_tree_entries; ii++) {        \
            tree_entries[ii].name = ((type*)buffer)[ii];     \
        }                                                    \
    }

void read_trees__velociraptor(int snapshot, halo_t* halos, int* n_halos, fof_group_t* fof_groups, int* n_fof_groups, int* index_lookup)
{
    // TODO: For the moment, I'll forgo chunking the read.  This will need to
    // be implemented in future though, as we ramp up the size of the
    // simulations...

    //! Tree entry struct
    typedef struct tree_entry_t {
        long ForestID;
        long Head;
        long hostHaloID;
        double Mass_200crit;
        double R_200crit;
        double Vmax;
        double Xc;
        double Yc;
        double Zc;
        double lambda_B;
        unsigned long ID;
        unsigned long npart;
    } tree_entry_t;

    mlog("Reading velociraptor trees for snapshot %d...", MLOG_OPEN, snapshot);

    // analyzer assertions
    assert(run_globals.mpi_rank >= 0);

    tree_entry_t* tree_entries = NULL;
    int n_tree_entries = 0;

    if (run_globals.mpi_rank == 0) {
        char fname[STRLEN];
        sprintf(fname, "%s/trees/VELOCIraptor.tree.t4.unifiedhalotree.withforest.snap.hdf.data",
            run_globals.params.SimulationDir);

        hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fd < 0) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(EXIT_FAILURE);
        }

        char snap_group_name[9];
        sprintf(snap_group_name, "Snap_%03d", snapshot);
        hid_t snap_group = H5Gopen(fd, snap_group_name, H5P_DEFAULT);

        H5LTget_attribute_int(fd, snap_group_name, "NHalos", &n_tree_entries);

        tree_entries = malloc(sizeof(tree_entry_t) * n_tree_entries);

        void* buffer = malloc(n_tree_entries * sizeof(long));

        READ_TREE_ENTRY_PROP(ForestID, long, H5T_NATIVE_LONG);
        READ_TREE_ENTRY_PROP(Head, long, H5T_NATIVE_LONG);
        READ_TREE_ENTRY_PROP(hostHaloID, long, H5T_NATIVE_LONG);
        READ_TREE_ENTRY_PROP(Mass_200crit, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(R_200crit, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(Vmax, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(Xc, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(Yc, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(Zc, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(lambda_B, double, H5T_NATIVE_DOUBLE);
        READ_TREE_ENTRY_PROP(ID, unsigned long, H5T_NATIVE_ULONG);
        READ_TREE_ENTRY_PROP(npart, unsigned long, H5T_NATIVE_ULONG);

        free(buffer);

        // convert units
        for (int ii = 0; ii < n_tree_entries; ii++) {
            tree_entries[ii].Mass_200crit *= 1e-10;
        }

        H5Gclose(snap_group);
        H5Fclose(fd);
    }

    MPI_Bcast(&n_tree_entries, 1, MPI_INT, 0, run_globals.mpi_comm);
    size_t _nbytes = sizeof(tree_entry_t) * n_tree_entries;
    if (run_globals.mpi_rank > 0)
        tree_entries = malloc(_nbytes);
    MPI_Bcast(tree_entries, (int)_nbytes, MPI_BYTE, 0, run_globals.mpi_comm);

    *n_halos = 0;
    *n_fof_groups = 0;
    for (int ii = 0; ii < n_tree_entries; ++ii) {
        bool keep_this_halo = true;

        if ((run_globals.RequestedForestId != NULL)
            && (bsearch(&(tree_entries[ii].ForestID), run_globals.RequestedForestId,
                   (size_t)run_globals.NRequestedForests, sizeof(int), compare_ints))
                == NULL)
            keep_this_halo = false;

        if (keep_this_halo) {
            tree_entry_t tree_entry = tree_entries[ii];
            halo_t* halo = &(halos[*n_halos]);

            halo->ID = tree_entry.ID;
            halo->DescIndex = id_to_ind(tree_entry.Head);
            halo->NextHaloInFOFGroup = NULL;
            halo->Type = tree_entry.hostHaloID == -1 ? 0 : 1;
            halo->SnapOffset = id_to_snap(tree_entry.Head) - snapshot;

            if (index_lookup)
                index_lookup[*n_halos] = ii;

            if (halo->Type == 0) {
                fof_group_t* fof_group = &fof_groups[*n_fof_groups];

                // TODO: What masses and radii should I use for centrals (inclusive vs. exclusive etc.)?
                fof_group->Mvir = tree_entry.Mass_200crit;
                fof_group->Rvir = tree_entry.R_200crit;
                fof_group->Vvir = -1;
                fof_group->FOFMvirModifier = 1.0;

                convert_input_virial_props(&fof_group->Mvir, &fof_group->Rvir, &fof_group->Vvir,
                    &fof_group->FOFMvirModifier, -1, snapshot, true);

                halo->FOFGroup = &(fof_groups[*n_fof_groups]);
                fof_groups[(*n_fof_groups)++].FirstHalo = halo;
            }
            else {
                // We can take advantage of the fact that host halos always seem to appear before their subhalos in the
                // trees to immediately connect FOF group members.
                int host_index = id_to_ind(tree_entry.hostHaloID);

                if (index_lookup)
                    host_index = find_original_index(host_index, index_lookup, *n_halos);

                assert(host_index > -1);
                assert(host_index < *n_halos);

                halo_t* prev_halo = &halos[host_index];
                halo->FOFGroup = prev_halo->FOFGroup;

                while (prev_halo->NextHaloInFOFGroup != NULL)
                    prev_halo = prev_halo->NextHaloInFOFGroup;

                prev_halo->NextHaloInFOFGroup = halo;
            }

            halo->Len = (int)tree_entry.npart;
            halo->Pos[0] = (float)tree_entry.Xc;
            halo->Pos[1] = (float)tree_entry.Yc;
            halo->Pos[2] = (float)tree_entry.Zc;
            halo->Vmax = (float)tree_entry.Vmax;

            // TODO: What masses and radii should I use for satellites (inclusive vs. exclusive etc.)?
            halo->Mvir = -1;
            halo->Rvir = -1;
            halo->Vvir = -1;
            convert_input_virial_props(&halo->Mvir, &halo->Rvir, &halo->Vvir, NULL, halo->Len, snapshot, false);

            // TODO: Ask Pascal for real ang mom vectors
            // N.B. See calculation of spin which has been hacked!
            halo->AngMom[0] = halo->AngMom[1] = 0;
            halo->AngMom[2] = tree_entry.lambda_B;

            halo->Galaxy = NULL;

            (*n_halos)++;
        }
    }

    free(tree_entries);

    mlog("...done", MLOG_CLOSE);
}
