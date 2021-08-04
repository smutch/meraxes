#include <assert.h>
#include <hdf5_hl.h>
#include <math.h>

#include "debug.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "modifiers.h"
#include "read_halos.h"
#include "tree_flags.h"
#include "virial_properties.h"

trees_info_t read_trees_info__velociraptor(const int snapshot)
{
  trees_info_t trees_info;

  if (run_globals.mpi_rank == 0) {
    // TODO: This could maybe only ever be done once and stored in run_globals.
    char fname[STRLEN + 34];
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
  return (int)(((uint64_t)id % (uint64_t)1e12) - 1);
}

static int id_to_snap(long id)
{
  return (int)(id / 1e12l);
}

inline static void convert_input_virial_props(double* Mvir,
                                              double* Rvir,
                                              double* Vvir,
                                              double* FOFMvirModifier,
                                              const int len,
                                              const int snapshot,
                                              const bool fof_flag)
{
  // Update the virial properties for subhalos
  if (*Mvir == -1) {
    assert(len > 0);
    *Mvir = calculate_Mvir(*Mvir, len);
  } else {
    if (fof_flag && (run_globals.RequestedMassRatioModifier == 1)) {
      // Modifier the FoF mass and update the virial radius
      assert(FOFMvirModifier != NULL);
      *FOFMvirModifier =
        interpolate_modifier(run_globals.mass_ratio_modifier, log10(*Mvir / run_globals.params.Hubble_h) + 10.0);
      *Mvir *= *FOFMvirModifier;
    }
  }

  if (*Rvir == -1)
    *Rvir = calculate_Rvir(*Mvir, snapshot);

  if (*Vvir == -1)
    *Vvir = calculate_Vvir(*Mvir, *Rvir);
}

void read_trees__velociraptor(int snapshot,
                              halo_t* halos,
                              int* n_halos,
                              fof_group_t* fof_groups,
                              int* n_fof_groups,
                              int* index_lookup)
{
  // TODO: For the moment, I'll forgo chunking the read.  This will need to
  // be implemented in future though, as we ramp up the size of the trees

  //! Tree entry struct
  typedef struct tree_entry_t
  {
    long ForestID;
    long Head;
    long Tail;
    long hostHaloID;
    double Mass_200crit;
    double Mass_FOF;
    double Mass_tot;
    double R_200crit;
    double Vmax;
    double Xc;
    double Yc;
    double Zc;
    double VXc;
    double VYc;
    double VZc;
    double Lx;
    double Ly;
    double Lz;
    unsigned long ID;
    unsigned long npart;
  } tree_entry_t;

  // simulations...

  mlog("Reading velociraptor trees for snapshot %d...", MLOG_OPEN, snapshot);

  int n_tree_entries = 0;
  hid_t fd = -1;
  hid_t snap_group = -1;
  void* property_buffer = NULL;
  double mass_unit_to_internal = 1.0;
  double scale_factor = -999.;

  *n_halos = 0;
  *n_fof_groups = 0;

  int buffer_size = 10000; // NOTE: Arbitrary. Should be a multiple of the chunk size of arrays in file ideally?
  tree_entry_t* tree_entries = malloc(sizeof(tree_entry_t) * buffer_size);

  if (run_globals.mpi_rank == 0) {
    char fname[STRLEN * 2 + 8];
    sprintf(fname, "%s/trees/%s", run_globals.params.SimulationDir, run_globals.params.CatalogFilePrefix);

    fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(EXIT_FAILURE);
    }

    char snap_group_name[9];
    sprintf(snap_group_name, "Snap_%03d", snapshot);
    snap_group = H5Gopen(fd, snap_group_name, H5P_DEFAULT);

    H5LTget_attribute_int(fd, snap_group_name, "NHalos", &n_tree_entries);

    property_buffer = malloc(buffer_size * sizeof(long));

    // check the units
    H5LTget_attribute_double(fd, "Header/Units", "Mass_unit_to_solarmass", &mass_unit_to_internal);
    mass_unit_to_internal /= 1.0e10;
    H5LTget_attribute_double(fd, snap_group_name, "scalefactor", &scale_factor);
  }

  MPI_Bcast(&n_tree_entries, 1, MPI_INT, 0, run_globals.mpi_comm);

  int n_read = 0;
  int n_to_read = buffer_size > n_tree_entries ? n_tree_entries : buffer_size;
  while (n_read < n_tree_entries) {
    int n_remaining = n_tree_entries - n_read;
    if (n_remaining < n_to_read) {
      n_to_read = n_remaining;
    }

    if (run_globals.mpi_rank == 0) {

      // select a hyperslab in the filespace
      hid_t fspace_id = H5Screate_simple(1, (hsize_t[1]){ n_tree_entries }, NULL);
      H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, (hsize_t[1]){ n_read }, NULL, (hsize_t[1]){ n_to_read }, NULL);
      hid_t memspace_id = H5Screate_simple(1, (hsize_t[1]){ n_to_read }, NULL);

#define READ_TREE_ENTRY_PROP(name, type, h5type)                                                                       \
  {                                                                                                                    \
    hid_t dset_id = H5Dopen(snap_group, #name, H5P_DEFAULT);                                                           \
    herr_t status = H5Dread(dset_id, h5type, memspace_id, fspace_id, H5P_DEFAULT, property_buffer);                    \
    assert(status >= 0);                                                                                               \
    H5Dclose(dset_id);                                                                                                 \
    for (int ii = 0; ii < n_to_read; ii++) {                                                                           \
      tree_entries[ii].name = ((type*)property_buffer)[ii];                                                            \
    }                                                                                                                  \
  }

      // TODO(trees): Read tail.  If head<->tail then first progenitor line, else it's a merger.  We should populate the
      // new halo and then do a standard merger prescription.

      READ_TREE_ENTRY_PROP(ForestID, long, H5T_NATIVE_LONG);
      READ_TREE_ENTRY_PROP(Head, long, H5T_NATIVE_LONG);
      READ_TREE_ENTRY_PROP(Tail, long, H5T_NATIVE_LONG);
      READ_TREE_ENTRY_PROP(hostHaloID, long, H5T_NATIVE_LONG);
      READ_TREE_ENTRY_PROP(Mass_200crit, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Mass_FOF, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Mass_tot, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(R_200crit, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Vmax, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Xc, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Yc, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Zc, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(VXc, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(VYc, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(VZc, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Lx, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Ly, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(Lz, double, H5T_NATIVE_DOUBLE);
      READ_TREE_ENTRY_PROP(ID, unsigned long, H5T_NATIVE_ULONG);
      READ_TREE_ENTRY_PROP(npart, unsigned long, H5T_NATIVE_ULONG);

      H5Sclose(memspace_id);
      H5Sclose(fspace_id);

      double hubble_h = run_globals.params.Hubble_h;
      for (int ii = 0; ii < n_to_read; ii++) {
        tree_entries[ii].Mass_200crit *= hubble_h * mass_unit_to_internal;
        tree_entries[ii].Mass_FOF *= hubble_h * mass_unit_to_internal;
        tree_entries[ii].Mass_tot *= hubble_h * mass_unit_to_internal;
        tree_entries[ii].R_200crit *= hubble_h;
        tree_entries[ii].Xc *= hubble_h / scale_factor;
        tree_entries[ii].Yc *= hubble_h / scale_factor;
        tree_entries[ii].Zc *= hubble_h / scale_factor;
        tree_entries[ii].VXc /= scale_factor;
        tree_entries[ii].VYc /= scale_factor;
        tree_entries[ii].VZc /= scale_factor;
        tree_entries[ii].Lx *= hubble_h * hubble_h * mass_unit_to_internal;
        tree_entries[ii].Ly *= hubble_h * hubble_h * mass_unit_to_internal;
        tree_entries[ii].Lz *= hubble_h * hubble_h * mass_unit_to_internal;

        // TEMPORARY HACK
        double box_size = run_globals.params.BoxSize;
        if (tree_entries[ii].Xc < 0.0)
          tree_entries[ii].Xc = 0.0;
        if (tree_entries[ii].Xc > box_size)
          tree_entries[ii].Xc = box_size;
        if (tree_entries[ii].Yc < 0.0)
          tree_entries[ii].Yc = 0.0;
        if (tree_entries[ii].Yc > box_size)
          tree_entries[ii].Yc = box_size;
        if (tree_entries[ii].Zc < 0.0)
          tree_entries[ii].Zc = 0.0;
        if (tree_entries[ii].Zc > box_size)
          tree_entries[ii].Zc = box_size;

#ifdef DEBUG
        assert((tree_entries[ii].Xc <= box_size) && (tree_entries[ii].Xc >= 0.0));
        assert((tree_entries[ii].Yc <= box_size) && (tree_entries[ii].Yc >= 0.0));
        assert((tree_entries[ii].Zc <= box_size) && (tree_entries[ii].Zc >= 0.0));
#endif
      }
    }

    size_t _nbytes = sizeof(tree_entry_t) * n_to_read;
    MPI_Bcast(tree_entries, (int)_nbytes, MPI_BYTE, 0, run_globals.mpi_comm);

    for (int ii = 0; ii < n_to_read; ++ii) {
      bool keep_this_halo = true;

      if ((run_globals.RequestedForestId != NULL) && (bsearch(&(tree_entries[ii].ForestID),
                                                              run_globals.RequestedForestId,
                                                              (size_t)run_globals.NRequestedForests,
                                                              sizeof(long),
                                                              compare_longs)) == NULL)
        keep_this_halo = false;

      if (keep_this_halo) {
        tree_entry_t tree_entry = tree_entries[ii];
        halo_t* halo = &(halos[*n_halos]);

        halo->ID = tree_entry.ID;
        halo->DescIndex = id_to_ind(tree_entry.Head);

        if (run_globals.params.FlagIgnoreProgIndex)
          halo->ProgIndex = -1;
        else
          halo->ProgIndex = id_to_ind(tree_entry.Tail);

        halo->NextHaloInFOFGroup = NULL;
        halo->Type = tree_entry.hostHaloID == -1 ? 0 : 1;
        halo->SnapOffset = id_to_snap(tree_entry.Head) - snapshot;

        // Any other tree flags need to be set using both the current and
        // progenitor halo information (stored in the galaxy), therefore we
        // need to leave setting those until later...
        if (run_globals.params.FlagIgnoreProgIndex)
          halo->TreeFlags = TREE_CASE_NO_PROGENITORS;
        else
          halo->TreeFlags = (unsigned long)tree_entry.Tail != tree_entry.ID ? 0 : TREE_CASE_NO_PROGENITORS;

        // Here we have a cyclic pointer, indicating that this halo's life ends here
        if ((unsigned long)tree_entry.Head == tree_entry.ID)
          halo->DescIndex = -1;

        if (index_lookup)
          index_lookup[*n_halos] = ii + n_read;

        // TODO: What masses and radii should I use for centrals (inclusive vs. exclusive etc.)?
        if (halo->Type == 0) {
          fof_group_t* fof_group = &fof_groups[*n_fof_groups];

          if (tree_entry.Mass_200crit <= 0) {
            // This "halo" is not above the virial threshold!  Use
            // proxy masses, but flag this fact so we know not to do
            // any or allow any hot halo to exist.
            halo->TreeFlags |= TREE_CASE_BELOW_VIRIAL_THRESHOLD;
            fof_group->Mvir = tree_entry.Mass_FOF;
            fof_group->Rvir = -1;
            // } else if (tree_entry.Mass_200crit < tree_entry.Mass_tot){
            // // The central subhalo has a proxy mass larger than the FOF
            // // group. Entirely possible for non-virialised and relaxed
            // // halos but doesn't really lead to internal consistency.
            // // Let's therefore just set the FOF virial mass to be that
            // // central subhalo proxy mass.
            // fof_group->Mvir = tree_entry.Mass_FOF;
            // fof_group->Rvir = -1;
          } else {
            fof_group->Mvir = tree_entry.Mass_200crit;
            fof_group->Rvir = tree_entry.R_200crit;
          }
          fof_group->Vvir = -1;
          fof_group->FOFMvirModifier = 1.0;

          convert_input_virial_props(
            &fof_group->Mvir, &fof_group->Rvir, &fof_group->Vvir, &fof_group->FOFMvirModifier, -1, snapshot, true);

          halo->FOFGroup = &(fof_groups[*n_fof_groups]);
          fof_groups[(*n_fof_groups)++].FirstHalo = halo;
        } else {
          // We can take advantage of the fact that host halos always
          // seem to appear before their subhalos (checked below) in the
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
        halo->Vel[0] = (float)tree_entry.VXc;
        halo->Vel[1] = (float)tree_entry.VYc;
        halo->Vel[2] = (float)tree_entry.VZc;
        halo->Vmax = (float)tree_entry.Vmax;

        // TODO: What masses and radii should I use for satellites (inclusive vs. exclusive etc.)?
        halo->Mvir = tree_entry.Mass_tot;
        halo->Rvir = -1;
        halo->Vvir = -1;
        convert_input_virial_props(&halo->Mvir, &halo->Rvir, &halo->Vvir, NULL, -1, snapshot, false);

        halo->AngMom[0] = (float)(tree_entry.Lx / tree_entry.Mass_tot);
        halo->AngMom[1] = (float)(tree_entry.Ly / tree_entry.Mass_tot);
        halo->AngMom[2] = (float)(tree_entry.Lz / tree_entry.Mass_tot);

        halo->Galaxy = NULL;

        (*n_halos)++;
      }
    }

    n_read += n_to_read;
  }

  free(tree_entries);

  if (run_globals.mpi_rank == 0) {
    free(property_buffer);
    H5Gclose(snap_group);
    H5Fclose(fd);
  }

  mlog("...done", MLOG_CLOSE);
}
