#include "meraxes.h"
#include <math.h>

void calculate_Mvir_crit(double redshift)
{
  // Calculate the critical Mvir value in each grid cell (ala Sobacchi & Mesinger 2013b)
  float *Mvir_crit = run_globals.reion_grids.Mvir_crit;

  int ReionGridDim = run_globals.params.ReionGridDim;
  float cell_Mvir_crit;
  int local_n_x = (int)(run_globals.reion_grids.slab_nix[SID.My_rank]);
  int local_n_cell = local_n_x * ReionGridDim * ReionGridDim;

  float ReionSMParam_m0 = run_globals.params.physics.ReionSMParam_m0;
  float ReionSMParam_a   = run_globals.params.physics.ReionSMParam_a;
  float ReionSMParam_b   = run_globals.params.physics.ReionSMParam_b;
  float ReionSMParam_c   = run_globals.params.physics.ReionSMParam_c;
  float ReionSMParam_d   = run_globals.params.physics.ReionSMParam_d;

  float *J_21_at_ion = run_globals.reion_grids.J_21_at_ionization;
  float *z_at_ion    = run_globals.reion_grids.z_at_ionization;

  // init
  memset(Mvir_crit, 0, sizeof(float) * local_n_cell);

  // Loop through each cell and calculate the value of Mvir_crit
  for (int ii = 0; ii < local_n_x; ii++)
  {
    for (int jj = 0; jj < ReionGridDim; jj++)
    {
      for (int kk = 0; kk < ReionGridDim; kk++)
      {
        // Initialise critical mass
        cell_Mvir_crit = 0.0;

        // If this cell was ionized in the past then calculate the critical
        // mass using the UVB feedback prescription of Sobacchi & Mesinger
        // 2013b
        if (z_at_ion[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] > redshift)
        {
          cell_Mvir_crit = ReionSMParam_m0 * pow(J_21_at_ion[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)], ReionSMParam_a)
            * pow((1.0 + redshift) / 10.0, ReionSMParam_b)
            * pow((1.0 - pow((1.0 + redshift) / (1.0 + z_at_ion[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)]), ReionSMParam_c)), ReionSMParam_d);
        }

        // Save the critical mass to the grid
        Mvir_crit[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] = cell_Mvir_crit;
      }
    }
  }
}


double tocf_modifier(galaxy_t *gal, double Mvir)
{
  return pow(2.0, -1.0 * gal->MvirCrit / Mvir);
}
