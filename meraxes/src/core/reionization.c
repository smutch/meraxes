#include "meraxes.h"
#include <fftw3.h>
#include <math.h>

void malloc_reionization_grids(run_globals_t *run_globals)
{
#ifdef USE_TOCF
  tocf_grids_t *grids = &(run_globals->tocf_grids);
  int n_cell = pow(tocf_params.HII_dim, 3);

  SID_log("Mallocing %.1f GB for required 21cmFAST grids...", SID_LOG_OPEN,
      (float)n_cell * (float)(sizeof(float)) * 6./(1024*1024*1024));

  grids->xH_grid         = (float *)fftwf_malloc(sizeof(float) * (size_t)n_cell);
  grids->stellar_grid    = (float *)fftwf_malloc(sizeof(float) * (size_t)n_cell);
  grids->sfr_grid        = (float *)fftwf_malloc(sizeof(float) * (size_t)n_cell);
  grids->z_at_ionization = (float *)fftwf_malloc(sizeof(float) * (size_t)n_cell);
  grids->J_at_ionization = (float *)fftwf_malloc(sizeof(float) * (size_t)n_cell);
  grids->Mvir_crit       = (float *)fftwf_malloc(sizeof(float) * (size_t)n_cell);

  SID_log(" ...done", SID_LOG_CLOSE);
#endif
}

void free_reionization_grids(run_globals_t *run_globals)
{
#ifdef USE_TOCF
  tocf_grids_t *grids = &(run_globals->tocf_grids);

  fftwf_free(grids->Mvir_crit);
  fftwf_free(grids->J_at_ionization);
  fftwf_free(grids->z_at_ionization);
  fftwf_free(grids->sfr_grid);
  fftwf_free(grids->stellar_grid);
  fftwf_free(grids->xH_grid);
#endif
}

static inline int find_cell(double pos, int xH_dim, double box_size)
{
  return (int)((pos/box_size)*(double)xH_dim);
}

void construct_stellar_grids(run_globals_t *run_globals)
{
#ifdef USE_TOCF
  galaxy_t *gal;
  int i, j, k;
  int xH_dim = tocf_params.HII_dim;
  int n_cell = pow(xH_dim,3);
  double box_size = (double)(run_globals->params.BoxSize);
  double Hubble_h = run_globals->params.Hubble_h;
  double UnitMass_in_g = run_globals->units.UnitMass_in_g;
  double UnitTime_in_s = run_globals->units.UnitTime_in_s;
  float *stellar_grid = run_globals->tocf_grids.stellar_grid;
  float *sfr_grid = run_globals->tocf_grids.sfr_grid;

  // init the grid
  for(int ii=0; ii<n_cell; ii++)
  {
    stellar_grid[ii] = 0.0;
    sfr_grid[ii] = 0.0;
  }

  // Loop through each valid galaxy and add its stellar mass to the appropriate cell
  gal = run_globals->FirstGal;
  while(gal!=NULL)
  {
    if((gal->Type < 3) && (!gal->ghost_flag))
    {
      i = find_cell(gal->Pos[0], xH_dim, box_size);
      j = find_cell(gal->Pos[1], xH_dim, box_size);
      k = find_cell(gal->Pos[2], xH_dim, box_size);
      stellar_grid[HII_R_FFT_INDEX(i,j,k)] += gal->StellarMass;
      sfr_grid[HII_R_FFT_INDEX(i,j,k)] += gal->Sfr;
    }
    gal = gal->Next;
  }

  // Do one final pass to put the grid in the correct units (Msol or Msol/yr)
  for(int ii=0; ii<n_cell; ii++)
  {
    stellar_grid[ii] *= 1.e10/Hubble_h;
    sfr_grid[ii] = sfr_grid[ii] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  }
#endif
}

void assign_ionization_to_halos(run_globals_t *run_globals, halo_t *halo, int n_halos, float *xH_grid, int xH_dim)
{
#ifdef USE_TOCF
  double box_size = (double)(run_globals->params.BoxSize);
  int i, j, k;

  SID_log("Assigning cell ionization values to halos...", SID_LOG_OPEN|SID_LOG_TIMER);
  SID_log("xH_dim = %d", SID_LOG_COMMENT, xH_dim);

  for(int i_halo=0; i_halo<n_halos; i_halo++)
  {
    i = find_cell(halo[i_halo].Pos[0], xH_dim, box_size);
    j = find_cell(halo[i_halo].Pos[1], xH_dim, box_size);
    k = find_cell(halo[i_halo].Pos[2], xH_dim, box_size);
    halo[i_halo].CellIonization = 1.0 - xH_grid[HII_R_INDEX(i,j,k)];
  }

  SID_log("...done", SID_LOG_CLOSE);
#endif
}

