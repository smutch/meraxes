#include "meraxes.h"
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

/*
   ==============================================================================
   MAJOR CODE REVISION by Paul Geil (October 2014)
   ==============================================================================

   - A significant numerical accuracy bug was resolved by performing calculations
    using an array of doubles and then casting it as an array of floats (as
    required by 21cmfast)

   - The TIAMAT velocity grids have not been analysed for anomalies.
 */



static inline void read_identifier(FILE *fin, bool skip_flag)
{
  char identifier[32];

  fread(identifier, sizeof(identifier), 1, fin);

  if (skip_flag)
    SID_log_error("Skipping grid: %s...", identifier);
  else
    SID_log("Reading grid: %s...", SID_LOG_COMMENT, identifier);
}


static int load_cached_deltax_slab(float *slab, int snapshot)
{
  if (run_globals.SnapshotDeltax[snapshot] != NULL)
  {
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[SID.My_rank];
    memcpy(slab, run_globals.SnapshotDeltax[snapshot], sizeof(float) * slab_n_complex * 2);
    return 0;
  }
  else
    return 1;
}


static int cache_deltax_slab(float *slab, int snapshot)
{
  if (run_globals.SnapshotDeltax[snapshot] == NULL)
  {
    float   **cache          = &run_globals.SnapshotDeltax[snapshot];
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[SID.My_rank];
    ptrdiff_t mem_size       = sizeof(float) * slab_n_complex * 2;

    *cache = fftwf_alloc_real(mem_size);
    memcpy(*cache, slab, mem_size);
    return 0;
  }
  else
    return 1;
}


int read_dm_grid(
  int    snapshot,
  int    i_grid,
  float *slab)
{
  // N.B. We assume in this function that the slab has the fftw3 inplace complex dft padding.

  run_params_t *params = &(run_globals.params);

  // Have we read this slab before?
  if (params->FlagInteractive && !load_cached_deltax_slab(slab, snapshot))
    return 0;

  char   fname[512];
  int    n_cell[3];
  double box_size[3];
  int    n_grids;
  int    ma_scheme;
  long   start_foffset;
  int    ReionGridDim = run_globals.params.ReionGridDim;

  // Construct the input filename
  sprintf(fname, "%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, snapshot);

  // Read the header
  if (SID.My_rank == 0)
  {
    FILE *fd;
    if((fd = fopen(fname, "rb")) == NULL)
    {
      fprintf(stderr, "Failed to open file: %s\n", fname);
      return EXIT_FAILURE;
    }

    fread(n_cell, sizeof(int), 3, fd);
    fread(box_size, sizeof(double), 3, fd);
    fread(&n_grids, sizeof(int), 1, fd);
    fread(&ma_scheme, sizeof(int), 1, fd);

    SID_log("Reading grid for snapshot %d", SID_LOG_OPEN | SID_LOG_TIMER, snapshot);
    SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, n_cell[0], n_cell[1], n_cell[2]);
    SID_log("box_size = [%.2f, %.2f, %.2f] cMpc/h", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);
    SID_log("ma_scheme = %d", SID_LOG_COMMENT, ma_scheme);

    if (n_grids != 4)
    {
      SID_log_error("n_grids != 4 as expected...");
      fclose(fd);
      ABORT(EXIT_FAILURE);
    }

    assert((n_cell[0] == n_cell[1]) && (n_cell[1] == n_cell[2])
           && "Input grids are not cubic!");

    // Compute the total number of elements in each grid
    int n_total_cell = n_cell[0] * n_cell[1] * n_cell[2];

    // Skip to the grid that we want
    // Note that we are expecting them to be in a particular order here
    for (int ii = 0; ii < i_grid; ii++)
    {
      read_identifier(fd, true);
      fseek(fd, sizeof(float) * n_total_cell, SEEK_CUR);
    }
    read_identifier(fd, false);

    start_foffset = ftell(fd);
    fclose(fd);
  }

  // share the needed information with all ranks
  MPI_Bcast(n_cell, 3, MPI_INT, 0, SID_COMM_WORLD);
  MPI_Bcast(box_size, 3, MPI_DOUBLE, 0, SID_COMM_WORLD);
  MPI_Bcast(&start_foffset, 1, MPI_LONG, 0, SID_COMM_WORLD);

  // Check if the grid in the file is higher resolution than we require
  double resample_factor = 1.;
  if ((n_cell[0] != ReionGridDim) || (n_cell[1] != ReionGridDim) || (n_cell[2] != ReionGridDim))
  {
    resample_factor = (double)ReionGridDim / (double)n_cell[0];
    if (resample_factor > 1.0001)
    {
      SID_log_error("The dark matter density grid in this file has a resolution less than that required! Aborting!");
      ABORT(EXIT_FAILURE);
    }
    SID_log("Using resample factor = %.3f", SID_LOG_COMMENT, resample_factor);
  }
  else
    resample_factor = 1;

  // Malloc the slab
  ptrdiff_t      slab_nix = run_globals.reion_grids.slab_nix[SID.My_rank];
  ptrdiff_t      slab_n_complex = run_globals.reion_grids.slab_n_complex[SID.My_rank];

  ptrdiff_t      slab_nix_file, slab_ix_start_file;
  ptrdiff_t      slab_n_complex_file = fftwf_mpi_local_size_3d(n_cell[0], n_cell[0], n_cell[0] / 2 + 1, SID_COMM_WORLD, &slab_nix_file, &slab_ix_start_file);
  fftwf_complex *slab_file           = fftwf_alloc_complex(slab_n_complex_file);
  ptrdiff_t      slab_ni_file        = slab_nix_file * n_cell[0] * n_cell[0];

  // Initialise (just in case!)
  for(int ii = 0; ii < slab_n_complex_file; ii++)
    slab_file[ii] = 0 + 0 * I;
  // N.B. factor of two for fftw padding
  for(int ii = 0; ii < slab_n_complex * 2; ii++)
    slab[ii] = 0.0;

  // Read in the slab for this rank
  MPI_File   fin         = NULL;
  MPI_Status status;
  MPI_Offset slab_offset = start_foffset / sizeof(float) + (slab_ix_start_file * n_cell[1] * n_cell[2]);

  MPI_File_open(SID_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
  MPI_File_set_view(fin, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

  ptrdiff_t chunk_size = slab_ni_file;
  int       n_reads    = 1;
  while(chunk_size > INT_MAX / 16)
  {
    chunk_size /= 2;
    n_reads    *= 2;
  }

  MPI_Offset offset = slab_offset;
  for(int ii = 0; ii < n_reads; ii++, offset += chunk_size)
  {
    MPI_File_read_at(fin, offset, &(((float *)slab_file)[chunk_size * ii]), chunk_size, MPI_FLOAT, &status);

    int count_check;
    MPI_Get_count(&status, MPI_FLOAT, &count_check);
    if (count_check != chunk_size)
    {
      SID_log_error("Failed to read correct number of elements on rank %d.", SID.My_rank);
      SID_log_error("Expected %d but read %d.", (int)chunk_size, count_check);
      ABORT(EXIT_FAILURE);
    }
  }
  MPI_File_close(&fin);

  // reorder the read slab for inplace fftw padding
  for(int ii = slab_nix_file - 1; ii >= 0; ii--)
    for(int jj = n_cell[0] - 1; jj >= 0; jj--)
      for(int kk = n_cell[0] - 1; kk >= 0; kk--)
        ((float *)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_PADDED)] = ((float *)slab_file)[grid_index(ii, jj, kk, n_cell[0], INDEX_REAL)];

  // smooth the grid if needed
  if (resample_factor < 1.0)
  {
    SID_log("Smoothing hi-res grid...", SID_LOG_OPEN | SID_LOG_TIMER);
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(n_cell[0], n_cell[0], n_cell[0], (float *)slab_file, slab_file, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
    // real space to k-space.
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse
    // FFT below
    long long total_n_cells_file = n_cell[0] * n_cell[0] * n_cell[0];
    for(int ii = 0; ii < slab_n_complex_file; ii++)
      slab_file[ii] /= (double)total_n_cells_file;
    filter(slab_file, slab_ix_start_file, slab_nix_file, n_cell[0], run_globals.params.BoxSize / (double)ReionGridDim / 2.0);

    plan = fftwf_mpi_plan_dft_c2r_3d(n_cell[0], n_cell[0], n_cell[0], slab_file, (float *)slab_file, SID_COMM_WORLD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    SID_log("...done", SID_LOG_CLOSE);
  }

  // Copy the read and smoothed slab into the padded fft slab (already allocated externally)
  int n_every = n_cell[0] / ReionGridDim;
  for (int ii = 0; ii < slab_nix; ii++)
  {
    int i_hr = n_every * ii;
    assert((i_hr > -1) && (i_hr < slab_nix_file));
    for (int jj = 0; jj < ReionGridDim; jj++)
    {
      int j_hr = n_every * jj;
      assert((j_hr > -1) && (j_hr < n_cell[0]));
      for (int kk = 0; kk < ReionGridDim; kk++)
      {
        int k_hr = n_every * kk;
        assert((k_hr > -1) && (k_hr < n_cell[0]));

        slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] = ((float *)slab_file)[grid_index(i_hr, j_hr, k_hr, n_cell[0], INDEX_PADDED)];
      }
    }
  }

  if (i_grid == 0)  // Density grid
  {
    // N.B. Hubble factor below to account for incorrect units in input DM grids!
    double mean = (double)run_globals.params.NPart * run_globals.params.PartMass / pow(box_size[0], 3) / run_globals.params.Hubble_h;

    // At this point grid holds the summed densities in each LR cell
    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    for (int ii = 0; ii < slab_nix; ii++)
      for (int jj = 0; jj < ReionGridDim; jj++)
        for (int kk = 0; kk < ReionGridDim; kk++)
        {
          float *val = &(slab[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]);
          *val = (float)(((double)*val / mean) - 1.);
        }
  }

  fftwf_free(slab_file);

  // Do we need to cache this slab?
  if (params->FlagInteractive)
    cache_deltax_slab(slab, snapshot);

  SID_log("...done", SID_LOG_CLOSE);

  // DEBUG
  // write_single_grid("output/debug.h5", slab, "deltax", true, true);

  return 0;
}


void free_grids_cache()
{
  if (run_globals.params.Flag_PatchyReion)
  {
    float **snapshot_deltax = run_globals.SnapshotDeltax;

    if (run_globals.params.FlagInteractive)
      for(int ii = 0; ii < run_globals.NStoreSnapshots; ii++)
        fftwf_free(snapshot_deltax[ii]);

    SID_free(SID_FARG snapshot_deltax);
  }
}