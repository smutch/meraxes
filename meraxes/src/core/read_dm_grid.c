#ifdef USE_TOCF

#include "meraxes.h"
#include <math.h>

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
  run_globals_t *run_globals,
  int            snapshot,   
  int            i_grid,     
  float         *grid)       
{

  // N.B. We assume in this function that the grid has the fftw3 inplace
  // complex dft padding.

  char fname[512];
  FILE *fin;
  int n_cell[3];
  double box_size[3];
  int n_grids;
  int ma_scheme;
  int n_elem;
  float val;
  float mean = 0.;
  float cell_volume = 0.;
  float resample_factor = 1.;
  run_params_t *params = &(run_globals->params);
  int HII_dim = tocf_params.HII_dim;

  // Construct the input filename 
  sprintf(fname,  "%s/%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, params->SimName, snapshot);
  // ... and open
  fin = fopen(fname, "rb");
  if (!(fin = fopen(fname, "rb")))
  {
    SID_log_error("Failed to open file: %s", fname);
    return(EXIT_FAILURE);
  }

  // Read the header
  fread(n_cell, sizeof(int), 3, fin);
  fread(box_size, sizeof(double), 3, fin);
  fread(&n_grids, sizeof(int), 1, fin);
  fread(&ma_scheme, sizeof(int), 1, fin);

  SID_log("Reading grid for snapshot %d", SID_LOG_OPEN, snapshot);
  SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, n_cell[0], n_cell[1], n_cell[2]);
  SID_log("box_size = [%.2f, %.2f, %.2f]", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);
  SID_log("ma_scheme = %d", SID_LOG_COMMENT, ma_scheme);

  if (n_grids != 4)
  {
    SID_log_error("n_grids != 4 as expected...");
    fclose(fin);
    return -1;
  }

  // Check if the grid in the file is higher resolution than we require
  if ((n_cell[0]!=HII_dim) || (n_cell[1]!=HII_dim) || (n_cell[2]!=HII_dim))
  {
    resample_factor = (float)HII_dim/(float)n_cell[0];
    if(resample_factor > 1.0001)
    {
      SID_log_error("The dark matter density grid in this file has a resolution less than that required! Aborting!");
      fclose(fin);
      ABORT(EXIT_FAILURE);
    }
    SID_log("Using resample factor = %.3f", SID_LOG_COMMENT, resample_factor);
  }
  else
    resample_factor = 1;

  // Compute the total number of elements in each grid
  n_elem = n_cell[0]*n_cell[1]*n_cell[2];

  // Read the grids
  // Note that we are expecting them to be in a particular order here
  for (int ii=0; ii<i_grid; ii++)
  {
    read_identifier(fin, true);
    fseek(fin, sizeof(float)*n_elem, SEEK_CUR);
  }
  read_identifier(fin, false);

  // init the grid
  for (int i=0; i<HII_dim; i++)
    for (int j=0; j<HII_dim; j++)
      for (int k=0; k<HII_dim; k++)
        *(grid + HII_R_FFT_INDEX(i,j,k)) = 0.;

  if(i_grid == 0)  // density grid
  {
    // calculate the volume of a single cell
    cell_volume = powf(box_size[0]/(float)n_cell[0], 3);
  
    // Read in the grid
    for (int i=0; i<n_cell[0]; i++)
      for (int j=0; j<n_cell[1]; j++)
        for (int k=0; k<n_cell[2]; k++)
        {
          fread(&val, sizeof(float), 1, fin);
          val *= cell_volume; // now we have the mass in the cell
          mean += val;
          *(grid + HII_R_FFT_INDEX((int)(i*resample_factor),(int)(j*resample_factor),(int)(k*resample_factor))) += val;
        }

    // In order to calculate the mean density we must actually take the mean of
    // the mass in each cell and then renormalise by the total volume.
    mean /= powf(box_size[0], 3);
    // SID_log("Calculated mean = %.2e :: theory = %.2e", SID_LOG_COMMENT, mean*1.e10*hlittle*hlittle, OMm*RHOcrit);

    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    cell_volume = powf(box_size[0]/(float)HII_dim, 3);
    for (int i=0; i<HII_dim; i++)
      for (int j=0; j<HII_dim; j++)
        for (int k=0; k<HII_dim; k++)
          *(grid + HII_R_FFT_INDEX(i,j,k)) = (*(grid + HII_R_FFT_INDEX(i,j,k))/(cell_volume * mean))-1.;
  } 
  else // velocity component grid
  {
  // Read in the grid
  for (int i=0; i<n_cell[0]; i++)
    for (int j=0; j<n_cell[1]; j++)
      for (int k=0; k<n_cell[2]; k++)
      {
        fread(&val, sizeof(float), 1, fin);
        *(grid + HII_R_FFT_INDEX((int)(i*resample_factor),(int)(j*resample_factor),(int)(k*resample_factor))) += val;
      }
  }
 
  SID_log("...done", SID_LOG_CLOSE);

  // Close the file
  fclose(fin);

  return 0;
}

#endif
