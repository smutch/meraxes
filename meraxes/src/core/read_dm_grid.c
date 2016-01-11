#ifdef USE_TOCF

#include "meraxes.h"
#include <math.h>


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


int read_dm_grid(
    int            snapshot,
    int            i_grid,
    float         *slab)
{
    // N.B. We assume in this function that the grid has the fftw3 inplace complex dft padding.

    char       fname[512];
    MPI_File   fin = NULL;
    MPI_Status status;
    int        n_cell[3];
    double     box_size[3];
    int        n_grids;
    int        ma_scheme;
    double     resample_factor   = 1.;
    int        HII_dim           = tocf_params.HII_dim;
    long       start_foffset;


    run_params_t *params = &(run_globals.params);

    // Construct the input filename
    sprintf(fname, "%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, snapshot);

    // Read the header
    if (SID.My_rank == 0)
    {
      FILE *fd = fopen(fname, "rb");
      if (!fin)
      {
        SID_log_error("Failed to open file: %s", fname);
        ABORT(EXIT_FAILURE);
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
    MPI_Bcast(n_cell, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(box_size, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&start_foffset, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    // Check if the grid in the file is higher resolution than we require
    if ((n_cell[0] != HII_dim) || (n_cell[1] != HII_dim) || (n_cell[2] != HII_dim))
    {
      resample_factor = (double)HII_dim / (double)n_cell[0];
      if (resample_factor > 1.0001)
      {
        SID_log_error("The dark matter density grid in this file has a resolution less than that required! Aborting!");
        ABORT(EXIT_FAILURE);
      }
      SID_log("Using resample factor = %.3f", SID_LOG_COMMENT, resample_factor);
    }
    else
    {
      resample_factor = 1;
    }

    // Malloc the slab
    ptrdiff_t slab_nix = tocf_params.slab_nix[SID.My_rank];
    
    ptrdiff_t slab_ni_file = slab_nix * n_cell[1] * n_cell[2];
    float *slab_file   = SID_calloc(sizeof(float) * slab_ni_file);

    // Initialise (just in case!)
    for(int ii=0; ii < slab_ni_file; ii++)
      slab_file[ii] = 0.0;
    for(int ii=0; ii < tocf_params.slab_n_complex[SID.My_rank]*2; ii++)
      slab[ii] = 0;


    // Read in the slab for this rank
    long slab_offset = start_foffset + tocf_params.slab_ix_start[SID.My_rank]*sizeof(float);
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
    MPI_File_read_at(fin, (MPI_Offset)slab_offset, slab_file, slab_ni_file, MPI_FLOAT, &status);
    MPI_File_close(&fin);

    // Copy the read slab into the padded fft slab (already allocated externally)
    // Regrid
    for (int i = 0; i < slab_nix; i++)
    {
      int i_lr = (int)(i * resample_factor);
      for (int j = 0; j < n_cell[1]; j++)
      {
        int j_lr = (int)(j * resample_factor);
        for (int k = 0; k < n_cell[2]; k++)
        {
          int k_lr = (int)(k * resample_factor);
          slab[grid_index(i_lr, j_lr, k_lr, HII_dim, INDEX_PADDED)] += slab_file[grid_index(i, j, k, HII_dim, INDEX_REAL)];
        }
      }
    }

    if (i_grid == 0)  // Density grid
    {
      // N.B. this function call destroys the ordering of the input array!
      double mean = accurate_sumf(slab_file, slab_ni_file);

      // Calculate the volume of a single high resolution cell
      double cell_volume = pow(box_size[0] / (double)n_cell[0], 3);

      // Mean density from high res grid
      mean *= cell_volume / pow(box_size[0], 3);

      // At this point grid holds the summed densities in each LR cell
      // Loop through again and calculate the overdensity
      // i.e. (rho - rho_mean)/rho_mean
      double cell_volume_ratio = pow(box_size[0] / (double)HII_dim, 3) / cell_volume;
      for (int i = 0; i < slab_nix; i++)
        for (int j = 0; j < HII_dim; j++)
          for (int k = 0; k < HII_dim; k++)
          {
            float *val = &(slab[grid_index(i, j, k, HII_dim, INDEX_PADDED)]);
            *val = (*val / (cell_volume_ratio * mean)) - 1.;
          }

    }

    SID_free(SID_FARG slab_file);

    SID_log("...done", SID_LOG_CLOSE);

    return 0;
}
#endif
