#include "meraxes.h"
#include <math.h>

void calculate_Mvir_crit(double redshift)
{
  // Calculate the critical Mvir value in each grid cell (ala Sobacchi & Mesinger 2013b)
  float *Mvir_crit = run_globals.tocf_grids.Mvir_crit;

  int HII_dim = tocf_params.HII_dim;
  float cell_Mvir_crit;
  float Hubble_h = (float)(run_globals.params.Hubble_h);
  int local_n_cell = (int)(tocf_params.slab_nix[SID.My_rank]) * HII_dim * HII_dim;

  float m_0_sm = tocf_params.m_0_sm;
  float a_sm   = tocf_params.a_sm;
  float b_sm   = tocf_params.b_sm;
  float c_sm   = tocf_params.c_sm;
  float d_sm   = tocf_params.d_sm;

  float *J_21_at_ion = run_globals.tocf_grids.J_21_at_ionization;
  float *z_at_ion    = run_globals.tocf_grids.z_at_ionization;

  // init
  memset(Mvir_crit, 0, sizeof(float) * local_n_cell);

  // Loop through each cell and calculate the value of Mvir_crit
  for (int ii = 0; ii < HII_dim; ii++)
  {
    for (int jj = 0; jj < HII_dim; jj++)
    {
      for (int kk = 0; kk < HII_dim; kk++)
      {
        // Initialise critical mass
        cell_Mvir_crit = 0.0;

        // If this cell was ionized in the past then calculate the critical
        // mass using the UVB feedback prescription of Sobacchi & Mesinger
        // 2013b
        if (z_at_ion[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)] > redshift)
        {
          cell_Mvir_crit = m_0_sm * pow(J_21_at_ion[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)], a_sm)
            * pow((1.0 + redshift) / 10.0, b_sm)
            * pow((1.0 - pow((1.0 + redshift) / (1.0 + z_at_ion[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)]), c_sm)), d_sm);

          // Put the mass back into internal units
          cell_Mvir_crit *= (Hubble_h / 1.0e10);
        }

        // Save the critical mass to the grid
        Mvir_crit[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)] = cell_Mvir_crit;
      }
    }
  }
}


double tocf_modifier(galaxy_t *gal, double Mvir)
{
  return pow(2.0, -1.0 * gal->MvirCrit / Mvir);
}
