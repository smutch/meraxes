#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "parameter_files/COSMOLOGY.H"
#include "parameter_files/INIT_PARAMS.H"
#include "parameter_files/ANAL_PARAMS.H"
#include "21cmfast.h"

// This file should contain a function that will read a GiggleZ/Tiamat grid and
// complain if the number of cells doesn't equal the resolution set in
// 21cmfast...

static inline void read_identifier(FILE *fin, bool skip_flag)
{
  char identifier[32];
  fread(identifier, sizeof(identifier), 1, fin);
  if (skip_flag)
    fprintf(stderr, "find_HII_bubbles:: Skipping grid: %s...\n", identifier);
  else
    fprintf(stderr, "find_HII_bubbles:: Reading grid: %s...\n", identifier);
}

int read_nbody_grid(
  tocf_params_struct *params,  
  int                 snapshot,
  int                 i_grid,  
  float              *grid)    
{

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

  // Construct the input filename 
  sprintf(fname,  "%s/%s/grids/snapshot_%03d_dark_grid.dat", params->sim_dir, params->sim_name, snapshot);
  // ... and open
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    fprintf(stderr, "find_HII_bubbles:: Failed to open %s\n", fname);
    return(EXIT_FAILURE);
  }

  // Read the header
  fread(n_cell, sizeof(int), 3, fin);
  fread(box_size, sizeof(double), 3, fin);
  fread(&n_grids, sizeof(int), 1, fin);
  fread(&ma_scheme, sizeof(int), 1, fin);

  fprintf(stderr, "find_HII_bubbles:: Reading grid for snapshot %d\n", snapshot);
  fprintf(stderr, "find_HII_bubbles:: n_cell = [%d, %d, %d]\n", n_cell[0], n_cell[1], n_cell[2]);
  fprintf(stderr, "find_HII_bubbles:: box_size = [%.2f, %.2f, %.2f]\n", box_size[0], box_size[1], box_size[2]);
  fprintf(stderr, "find_HII_bubbles:: ma_scheme = %d\n", ma_scheme);

  if (n_grids != 4)
  {
    fprintf(stderr, "find_HII_bubbles:: n_grids != 3 as expected...\n");
    fclose(fin);
    return -1;
  }
  if ((n_cell[0]!=HII_DIM) || (n_cell[1]!=HII_DIM) || (n_cell[2]!=HII_DIM))
  {
    // fprintf(stderr, "find_HII_bubbles:: At least one n_cell axis != INIT_PARAMS.h HII_DIM value...\n");
    // fclose(fin);
    // return -1;
    resample_factor = (float)HII_DIM/(float)n_cell[0];
    fprintf(stderr, "find_HII_bubbles:: Using resample factor = %.3f\n", resample_factor);
  }

  // Compute the total number of elements in each grid
  n_elem = n_cell[0]*n_cell[1]*n_cell[2];

  // and the volume of a single cell
  cell_volume = powf(box_size[0]/(float)n_cell[0], 3);
  
  // Read the grids
  // Note that we are expecting them to be in a particular order here
  for (int ii=0; ii<i_grid; ii++)
  {
    read_identifier(fin, true);
    fseek(fin, sizeof(float)*n_elem, SEEK_CUR);
  }
  read_identifier(fin, false);

  // init the grid
  for (int i=0; i<HII_DIM; i++)
    for (int j=0; j<HII_DIM; j++)
      for (int k=0; k<HII_DIM; k++)
        *(grid + HII_R_FFT_INDEX(i,j,k)) = 0.;

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
  fprintf(stderr, "find_HII_bubbles:: Calculated mean = %.2e :: theory = %.2e\n", mean*1.e10*hlittle*hlittle, OMm*RHOcrit);

  // Loop through again and calculate the overdensity
  // i.e. (rho - rho_mean)/rho_mean
  cell_volume = powf(box_size[0]/(float)HII_DIM, 3);
  for (int i=0; i<HII_DIM; i++)
    for (int j=0; j<HII_DIM; j++)
      for (int k=0; k<HII_DIM; k++)
        *(grid + HII_R_FFT_INDEX(i,j,k)) = (*(grid + HII_R_FFT_INDEX(i,j,k))/(cell_volume * mean))-1.;
 
  fprintf(stderr, "find_HII_bubbles:: ...done\n");

  // Close the file
  fclose(fin);

  return 0;
}
