#include <assert.h>
#include <complex.h>
#include <dirent.h>
#include <errno.h>
#include <fftw3-mpi.h>
#include <math.h>

#include "meraxes.h"
#include "misc_tools.h"
#include "read_grids.h"

/*
   ==============================================================================
   MAJOR CODE REVISION by Paul Geil (October 2014)
   ==============================================================================

   - A significant numerical accuracy bug was resolved by performing calculations
    using an array of doubles and then casting it as an array of floats (as
    required by 21cmfast)

   - The TIAMAT velocity grids have not been analysed for anomalies.
 */

static inline void read_identifier(FILE* fin, const bool skip_flag, const int n_cell[3])
{
  char identifier[32];

  fread(identifier, sizeof(identifier), 1, fin);

  if (skip_flag) {
    // mlog("Skipping grid: %s...", identifier);
    fseek(fin, sizeof(float) * n_cell[0] * n_cell[1] * n_cell[2], SEEK_CUR);
  } else
    mlog("Reading grid: %s...", MLOG_MESG, identifier);
}

static long read_header(const char* fname,
                        const int snapshot,
                        int n_cell[3],
                        double box_size[3],
                        int* n_grids,
                        int* ma_scheme,
                        const enum grid_prop property)
{
  // Read the header
  if (run_globals.mpi_rank == 0) {
    long start_foffset;

    FILE* fd;
    if ((fd = fopen(fname, "rb")) == NULL) {
      fprintf(stderr, "Failed to open file: %s\n", fname);
      return EXIT_FAILURE;
    }

    fread(n_cell, sizeof(int), 3, fd);
    fread(box_size, sizeof(double), 3, fd);
    fread(n_grids, sizeof(int), 1, fd);
    fread(ma_scheme, sizeof(int), 1, fd);

    mlog("Reading grid for snapshot %d", MLOG_OPEN | MLOG_TIMERSTART, snapshot);
    mlog("n_cell = [%d, %d, %d]", MLOG_MESG, n_cell[0], n_cell[1], n_cell[2]);
    mlog("box_size = [%.2f, %.2f, %.2f] cMpc/h", MLOG_MESG, box_size[0], box_size[1], box_size[2]);
    mlog("ma_scheme = %d", MLOG_MESG, *ma_scheme);

    if (*n_grids != 4) {
      mlog_error("n_grids != 4 as expected...");
      fclose(fd);
      ABORT(EXIT_FAILURE);
    }

    assert((n_cell[0] == n_cell[1]) && (n_cell[1] == n_cell[2]) && "Input grids are not cubic!");

    // Note that we are expecting the first grid to be the density grid
    for (int ii = 0; ii < property; ++ii)
      read_identifier(fd, true, n_cell);

    read_identifier(fd, false, n_cell);

    start_foffset = ftell(fd);
    fclose(fd);

    return start_foffset;
  } else {
    return 0;
  }
}

int read_grid__gbptrees(const enum grid_prop property, const int snapshot, float* slab)
{
  // N.B. We assume in this function that the slab has the fftw3 inplace complex dft padding.

  run_params_t* params = &(run_globals.params);

  // Have we read this slab before?
  if ((params->FlagInteractive || params->FlagMCMC) && !load_cached_slab(slab, snapshot, property))
    return 0;

  if ((property == X_VELOCITY) || (property == Y_VELOCITY) || (property == Z_VELOCITY)) {

    if (run_globals.params.TsVelocityComponent < 1 || run_globals.params.TsVelocityComponent > 3) {
      mlog("Not a valid velocity direction: x=1, y=2, z=3", MLOG_MESG);
      ABORT(EXIT_FAILURE);
    }

    if (run_globals.params.Flag_ConstructLightcone && run_globals.params.TsVelocityComponent != 3) {
      mlog("Light-cone is generated along the z-direction, therefore the velocity component should be in the "
           "z-direction (i.e 3).",
           MLOG_MESG);
      ABORT(EXIT_FAILURE);
    }
  }

  char fname[512];
  int n_cell[3];
  double box_size[3];
  int n_grids;
  int ma_scheme;
  long start_foffset;
  int ReionGridDim = run_globals.params.ReionGridDim;

  // Construct the input filename by first testing to see if there are
  // pre-computed grids of the required resolution.  If not then we will just
  // read the highest res grids available and down sample them.
  char dirname[512];
  sprintf(dirname, "%s/grids/resampled/N%d", params->SimulationDir, run_globals.params.ReionGridDim);
  DIR* dir = opendir(dirname);
  if (dir) {
    closedir(dir);
    sprintf(fname, "%s/snapshot_%03d_dark_grid.dat", dirname, snapshot);
  } else {
    sprintf(fname, "%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, snapshot);
  }

  start_foffset = read_header(fname, snapshot, n_cell, box_size, &n_grids, &ma_scheme, property);

  // share the needed information with all ranks
  MPI_Bcast(n_cell, 3, MPI_INT, 0, run_globals.mpi_comm);
  MPI_Bcast(box_size, 3, MPI_DOUBLE, 0, run_globals.mpi_comm);
  MPI_Bcast(&start_foffset, 1, MPI_LONG, 0, run_globals.mpi_comm);

  // Check if the grid in the file is higher resolution than we require
  double resample_factor = calc_resample_factor(n_cell);

#ifdef DEBUG
  mlog("Resample factor = %.3g", MLOG_MESG, resample_factor);
#endif

  // Malloc the slab
  ptrdiff_t slab_nix = run_globals.reion_grids.slab_nix[run_globals.mpi_rank];
  ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

  ptrdiff_t slab_nix_file, slab_ix_start_file;
  ptrdiff_t slab_n_complex_file = fftwf_mpi_local_size_3d(
    n_cell[0], n_cell[0], n_cell[0] / 2 + 1, run_globals.mpi_comm, &slab_nix_file, &slab_ix_start_file);
  fftwf_complex* slab_file = fftwf_alloc_complex((size_t)slab_n_complex_file);
  ptrdiff_t slab_ni_file = slab_nix_file * n_cell[0] * n_cell[0];

  // Initialise (just in case!)
  for (int ii = 0; ii < slab_n_complex_file; ii++)
    slab_file[ii] = 0 + 0I;
  // N.B. factor of two for fftw padding
  for (int ii = 0; ii < slab_n_complex * 2; ii++)
    slab[ii] = 0.0;

  // Read in the slab for this rank
  MPI_File fin = NULL;
  MPI_Status status;
  MPI_Offset slab_offset = start_foffset / sizeof(float) + (slab_ix_start_file * n_cell[1] * n_cell[2]);

  MPI_File_open(run_globals.mpi_comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
  MPI_File_set_view(fin, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

  ptrdiff_t chunk_size = slab_ni_file;
  int n_reads = 1;
  while (chunk_size > INT_MAX / 16) {
    chunk_size /= 2;
    n_reads *= 2;
  }

  MPI_Offset offset = slab_offset;
  for (int ii = 0; ii < n_reads; ii++, offset += chunk_size) {
    MPI_File_read_at(fin, offset, &(((float*)slab_file)[chunk_size * ii]), (int)chunk_size, MPI_FLOAT, &status);

    int count_check;
    MPI_Get_count(&status, MPI_FLOAT, &count_check);
    if (count_check != chunk_size) {
      mlog_error("Failed to read correct number of elements on rank %d.", run_globals.mpi_rank);
      mlog_error("Expected %d but read %d.", (int)chunk_size, count_check);
      ABORT(EXIT_FAILURE);
    }
  }
  MPI_File_close(&fin);

  // reorder the read slab for inplace fftw padding
  for (int ii = (int)(slab_nix_file - 1); ii >= 0; ii--)
    for (int jj = n_cell[0] - 1; jj >= 0; jj--)
      for (int kk = n_cell[0] - 1; kk >= 0; kk--)
        ((float*)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_PADDED)] =
          ((float*)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_REAL)];

  // smooth the grid and subsample if needed
  smooth_grid(resample_factor, n_cell, slab_file, slab_n_complex_file, slab_ix_start_file, slab_nix_file);
  subsample_grid(resample_factor, n_cell, (int)slab_ix_start_file, (int)slab_nix_file, (float*)slab_file, slab);

  if (property == DENSITY) {
    // N.B. Hubble factor below to account for incorrect units in input DM grids!
    float mean_inv =
      pow(box_size[0], 3) * run_globals.params.Hubble_h / run_globals.params.NPart / run_globals.params.PartMass;

    // At this point grid holds the summed densities in each LR cell
    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    for (int ii = 0; ii < slab_nix; ii++)
      for (int jj = 0; jj < ReionGridDim; jj++)
        for (int kk = 0; kk < ReionGridDim; kk++) {
          float* val = &(slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]);
          // the fmax check here tries to account for negative densities introduced by fftw rounding / aliasing effects
          *val = fmaxf(*val * mean_inv - 1.0, -1.0);
        }
  }

  fftwf_free(slab_file);

  // Do we need to cache this slab?
  if (params->FlagInteractive || params->FlagMCMC) {
    cache_slab(slab, snapshot, property);
  }

  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  return 0;
}
