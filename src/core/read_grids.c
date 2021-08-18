#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>

#include "meraxes.h"
#include "misc_tools.h"
#include "read_grids.h"
#include "reionization.h"

void read_grid(const enum grid_prop property, const int snapshot, float* slab)
{
  // Read in the dark matter density grid
  switch (run_globals.params.TreesID) {
    case VELOCIRAPTOR_TREES:
      read_grid__velociraptor(property, snapshot, slab);
      break;
    case GBPTREES_TREES:
      read_grid__gbptrees(property, snapshot, slab);
      break;
    default:
      mlog_error("Unrecognised input trees identifier (TreesID).");
      break;
  }
}

double calc_resample_factor(int n_cell[3])
{
  int ReionGridDim = run_globals.params.ReionGridDim;

  // Check that n_cell is divisible by ReionGridDim. We need this to ensure
  // that grid points are evenly spaced in the volume after resampling.
  if (n_cell[0] % ReionGridDim != 0) {
    mlog_error(
      "n_cell (%d) is not divisble by ReionGridDim (%d). This is required for downsampling.", n_cell[0], ReionGridDim);
    ABORT(EXIT_FAILURE);
  }

  // Check if the grid in the file is higher resolution than we require
  if (n_cell[0] != ReionGridDim) {
    double resample_factor = (double)ReionGridDim / (double)n_cell[0];
    if (resample_factor > 1.0001) {
      mlog_error("Grid has a resolution less than that required! Aborting!");
      ABORT(EXIT_FAILURE);
    }
    mlog("Using resample factor = %.3f", MLOG_MESG, resample_factor);
    return resample_factor;
  } else
    return 1.0;
}

void smooth_grid(double resample_factor,
                 int n_cell[3],
                 fftwf_complex* slab,
                 ptrdiff_t slab_n_complex,
                 ptrdiff_t slab_ix_start,
                 ptrdiff_t slab_nix)
{
  // TODO: Use wisdom and FFTW_PATIENT planning for this if requested (cf. malloc_reionisation_grids())

  if (resample_factor < 1.0) {
    mlog("Smoothing hi-res grid...", MLOG_OPEN | MLOG_TIMERSTART);
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(
      n_cell[0], n_cell[1], n_cell[2], (float*)slab, slab, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
    // real space to k-space.
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse
    // FFT below
    double total_n_cells = n_cell[0] * n_cell[1] * n_cell[2];
    for (int ii = 0; ii < (int)slab_n_complex; ii++)
      slab[ii] /= total_n_cells;

    filter(slab,
           (int)slab_ix_start,
           (int)slab_nix,
           n_cell[0],
           (float)(run_globals.params.BoxSize / (double)run_globals.params.ReionGridDim / 2.0),
           0); // NOTE: Real space top-hat hard-coded for this

    plan = fftwf_mpi_plan_dft_c2r_3d(
      n_cell[0], n_cell[1], n_cell[2], slab, (float*)slab, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
  }
}

void subsample_grid(double resample_factor, int n_cell[3], int ix_hi_start, int nix_hi, float* slab_file, float* slab)
{
  // we don't need to do anything in this case
  // TODO: This is potentially very wasteful if resample_factor == 1.
  //       Here is would be far better to not allocate slab_file in the first
  //       place in the calling function and then just do nothing here...
  if (resample_factor >= 1.0) {
    memcpy(slab, slab_file, sizeof(fftwf_complex) * run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
    return;
  }

  // Glossary
  // ix_hi: the x-axis index on the high res grid read from the input file
  // ix_lo: the x-axis index on the low res, subsampled, grid which will be used for the reionisation calculation
  // hi_rank: the rank which holds the current ix_hi value
  // lo_rank: the rank which will hold the current ix_lo value
  // slice: a single x-index cut through an array. i.e. one x-value and all of the y and z values

  // gather all of the slab_ix_start_file values
  int* ix_hi_start_allranks = calloc(run_globals.mpi_size, sizeof(int));
  int* nix_hi_allranks = calloc(run_globals.mpi_size, sizeof(int));

  MPI_Allgather(&ix_hi_start, 1, MPI_INT, ix_hi_start_allranks, 1, MPI_INT, run_globals.mpi_comm);
  MPI_Allgather(&nix_hi, 1, MPI_INT, nix_hi_allranks, 1, MPI_INT, run_globals.mpi_comm);

  int ReionGridDim = run_globals.params.ReionGridDim;
  int n_every = n_cell[0] / ReionGridDim;
  int slice_size = ReionGridDim * ReionGridDim;
  float* slice = calloc(slice_size, sizeof(float));
  int nix_lo = (int)run_globals.reion_grids.slab_nix[run_globals.mpi_rank];

  int ix_lo_start = (int)run_globals.reion_grids.slab_ix_start[0];
  int ix_lo_end = ix_lo_start + (int)run_globals.reion_grids.slab_nix[0] - 1;
  ix_hi_start = ix_hi_start_allranks[0];
  int ix_hi_end = ix_hi_start + nix_hi_allranks[0] - 1;

  for (int ix_hi = 0, ix_lo = 0, hi_rank = 0, lo_rank = 0; ix_hi < n_cell[0]; ix_hi += n_every, ++ix_lo) {

    if (ix_lo > ix_lo_end)
      lo_rank++;
    if (ix_hi > ix_hi_end)
      hi_rank++;

    ix_lo_start = (int)run_globals.reion_grids.slab_ix_start[lo_rank];
    ix_lo_end = ix_lo_start + (int)run_globals.reion_grids.slab_nix[lo_rank] - 1;
    ix_hi_start = ix_hi_start_allranks[hi_rank];
    ix_hi_end = ix_hi_start + nix_hi_allranks[hi_rank] - 1;

    int ix_lo_slab = ix_lo - ix_lo_start;
    int ix_hi_slab = ix_hi - ix_hi_start;

    if ((hi_rank == run_globals.mpi_rank) && (lo_rank != run_globals.mpi_rank)) {
      // This rank has the high res x-value but needs to send it somewhere else...

      // pack the slice
      for (int iy_lo = 0; iy_lo < ReionGridDim; iy_lo++) {
        int iy_hi = n_every * iy_lo;
        assert((iy_hi > -1) && (iy_hi < n_cell[0]));

        for (int iz_lo = 0; iz_lo < ReionGridDim; iz_lo++) {
          int iz_hi = n_every * iz_lo;
          assert((iz_hi > -1) && (iz_hi < n_cell[0]));

          slice[grid_index(0, iy_lo, iz_lo, ReionGridDim, INDEX_REAL)] =
            slab_file[grid_index(ix_hi_slab, iy_hi, iz_hi, n_cell[0], INDEX_PADDED)];
        }
      }

      // send it to the lo rank
      int ident = 400000 + 1000 * run_globals.mpi_rank + lo_rank;
      MPI_Send(slice, slice_size, MPI_FLOAT, lo_rank, ident, run_globals.mpi_comm);

    } else if ((hi_rank != run_globals.mpi_rank) && (lo_rank == run_globals.mpi_rank)) {
      // This rank needs the high res x-value but needs to get it from somewhere else...

      // receive the packed slice
      int ident = 400000 + 1000 * hi_rank + run_globals.mpi_rank;
      MPI_Recv(slice, slice_size, MPI_FLOAT, hi_rank, ident, run_globals.mpi_comm, MPI_STATUS_IGNORE);

      assert((ix_lo_slab >= 0) && (ix_lo_slab < nix_lo));

      // store it
      for (int iy_lo = 0; iy_lo < ReionGridDim; iy_lo++) {
        int iy_hi = n_every * iy_lo;
        assert((iy_hi > -1) && (iy_hi < n_cell[0]));

        for (int iz_lo = 0; iz_lo < ReionGridDim; iz_lo++) {
          int iz_hi = n_every * iz_lo;
          assert((iz_hi > -1) && (iz_hi < n_cell[0]));

          slab[grid_index(ix_lo_slab, iy_lo, iz_lo, ReionGridDim, INDEX_PADDED)] =
            slice[grid_index(0, iy_lo, iz_lo, ReionGridDim, INDEX_REAL)];
        }
      }

    } else if ((hi_rank == run_globals.mpi_rank) && (lo_rank == run_globals.mpi_rank)) {
      // This rank needs the high res x-value and already has it...

      assert((ix_lo_slab >= 0) && (ix_lo_slab < nix_lo));

      for (int iy_lo = 0; iy_lo < ReionGridDim; iy_lo++) {
        int iy_hi = n_every * iy_lo;
        assert((iy_hi > -1) && (iy_hi < n_cell[0]));

        for (int iz_lo = 0; iz_lo < ReionGridDim; iz_lo++) {
          int iz_hi = n_every * iz_lo;
          assert((iz_hi > -1) && (iz_hi < n_cell[0]));

          slab[grid_index(ix_lo_slab, iy_lo, iz_lo, ReionGridDim, INDEX_PADDED)] =
            slab_file[grid_index(ix_hi_slab, iy_hi, iz_hi, n_cell[0], INDEX_PADDED)];
        }
      }
    }

    // If we didn't match any of the above cases, we neither have nor want this slice.
  }

  free(slice);
  free(nix_hi_allranks);
  free(ix_hi_start_allranks);
}

int load_cached_slab(float* slab, int snapshot, const enum grid_prop property)
{
  float* cache = NULL;
  switch (property) {
    case DENSITY:
      cache = run_globals.SnapshotDeltax[snapshot];
      break;
    case X_VELOCITY:
    case Y_VELOCITY:
    case Z_VELOCITY:
      cache = run_globals.SnapshotVel[snapshot];
      break;
    default:
      mlog_error("Unrecognised grid property in load_cached_slab!");
      break;
  }

  if (cache != NULL) {
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];
    memcpy(slab, cache, sizeof(float) * slab_n_complex * 2);
    mlog("Loaded deltax slab from cache.", MLOG_MESG);
    return 0;
  } else
    return 1;
}

int cache_slab(float* slab, int snapshot, const enum grid_prop property)
{
  float** cache = NULL;
  switch (property) {
    case DENSITY:
      cache = &run_globals.SnapshotDeltax[snapshot];
      break;
    case X_VELOCITY:
    case Y_VELOCITY:
    case Z_VELOCITY:
      cache = &run_globals.SnapshotVel[snapshot];
      break;
    default:
      mlog_error("Unrecognised grid property in load_cached_slab!");
      break;
  }
  if (cache == NULL) {
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];
    ptrdiff_t mem_size = sizeof(float) * slab_n_complex * 2;

    *cache = fftwf_alloc_real((size_t)mem_size);
    memcpy(*cache, slab, mem_size);
    return 0;
  } else
    return 1;
}

void free_grids_cache()
{
  if (run_globals.params.Flag_PatchyReion) {
    float** snapshot_vel = run_globals.SnapshotVel;
    float** snapshot_deltax = run_globals.SnapshotDeltax;

    if (run_globals.params.FlagInteractive)
      for (int ii = 0; ii < run_globals.NStoreSnapshots; ii++) {
        fftwf_free(snapshot_vel[ii]);
        fftwf_free(snapshot_deltax[ii]);
      }

    free(snapshot_vel);
    free(snapshot_deltax);
  }
}
