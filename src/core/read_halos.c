#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>

static trees_info_t read_trees_info(const int snapshot)
{
  char fname[STRLEN];
  sprintf(fname, "%s/trees/VELOCIraptor.tree.t4.unifiedhalotree.withforest.snap.hdf.data", run_globals.params.SimulationDir);

  hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT); 
  if (fd < 0)
  {
    mlog("Failed to open file %s", MLOG_MESG, fname);
    ABORT(EXIT_FAILURE);
  }

  char grp_name[STRLEN];
  sprintf(grp_name, "Snap_%03d", snapshot);

  trees_info_t trees_info;
  H5LTget_attribute_int(fd, grp_name, "NHalos", &(trees_info.n_halos));

  H5Fclose(fd);

  return trees_info;
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
  if (run_globals.mpi_rank == 0)
  {
    char fname[STRLEN];
    sprintf(fname, "%s/trees/VELOCIraptor.tree.t4.unifiedhalotree.withforest.snap.hdf.data", run_globals.params.SimulationDir);

    hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT); 
    if (fd < 0)
    {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(EXIT_FAILURE);
    }

    int *forest_ids = read_forest_ids(fd);

    free(forest_ids);

    H5Fclose(fd);
  }

}


trees_info_t read_halos(const int snapshot, halo_t **halos)
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

  // if necessary, broadcast the tree file info
  MPI_Bcast(&trees_info, sizeof(trees_info_t), MPI_BYTE, 0, run_globals.mpi_comm);

  // TODO: deal with forests and set n_halos for each rank correctly
  if (run_globals.SelectForestsSwitch == true)
  {
    select_forests();
    run_globals.SelectForestsSwitch = false;
  }

  // Allocate the halo array
  if (*halos == NULL)
  {
    mlog("Allocating halo array with %d elements...", MLOG_MESG, trees_info.n_halos);
    *halos = malloc(sizeof(halo_t) * trees_info.n_halos);
  }

}
