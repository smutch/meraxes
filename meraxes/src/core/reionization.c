#ifdef USE_TOCF

#include "meraxes.h"
#include <fftw3.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>

// NOTE: I am following the indexing conventions off 21cmFAST here.  I have no
// idea what the hell the unsigned long long memory offset stuff is all about
// with the fftwf_complex arrays that have been cast as floats.  Why not still
// just use common indexing in square brackets? 
// i.e. ((float *)deltax)[k+(2*(HII_dim/2+1))*(j+HII_dim*i)]

void call_find_HII_bubbles(run_globals_t *run_globals, int snapshot, int nout_gals)
{
  // Thin wrapper round find_HII_bubbles
  tocf_grids_t *grids = &(run_globals->tocf_grids);

  SID_log("Getting ready to call find_HII_bubbles...", SID_LOG_OPEN);

  // Construct the stellar mass grid
  if(nout_gals>0)
    construct_stellar_grids(run_globals);
  else
  {
    SID_log("No galaxies present - skipping...", SID_LOG_CLOSE);
    return;
  }

  // Read in the dark matter density grid
  read_dm_grid(run_globals, snapshot, 0, (float *)(grids->deltax));

  SID_log("...done", SID_LOG_CLOSE);

  SID_log("Calling find_HII_bubbles...", SID_LOG_OPEN);
  // TODO: Fix if snapshot==0
  find_HII_bubbles(run_globals->ZZ[snapshot], run_globals->ZZ[snapshot-1],
      tocf_params.HII_eff_factor, tocf_params.ion_tvir_min, tocf_params.r_bubble_max, tocf_params.numcores,
      grids->xH,
      grids->stars,
      grids->stars_filtered,
      grids->deltax,
      grids->deltax_filtered,
      grids->z_at_ionization,
      grids->J_21_at_ionization,
      grids->J_21,
      NULL,
      NULL,
      NULL
      );

  SID_log("...done", SID_LOG_CLOSE);
}


void malloc_reionization_grids(run_globals_t *run_globals)
{
  tocf_grids_t *grids = &(run_globals->tocf_grids);

  SID_log("Mallocing %.1f GB for required 21cmFAST grids...", SID_LOG_OPEN,
      ((float)(HII_TOT_NUM_PIXELS * sizeof(float) * 5) +
      (float)(HII_KSPACE_NUM_PIXELS * sizeof(fftwf_complex) * 4))
      /(1024.*1024.*1024.));

  grids->xH                 = (float *)         fftwf_malloc(sizeof(float)        * HII_TOT_NUM_PIXELS);
  grids->stars              = (fftwf_complex *) fftwf_malloc(sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  grids->stars_filtered     = (fftwf_complex *) fftwf_malloc(sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  grids->deltax             = (fftwf_complex *) fftwf_malloc(sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  grids->deltax_filtered    = (fftwf_complex *) fftwf_malloc(sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  grids->z_at_ionization    = (float *)         fftwf_malloc(sizeof(float)        * HII_TOT_NUM_PIXELS);
  grids->J_21_at_ionization = (float *)         fftwf_malloc(sizeof(float)        * HII_TOT_NUM_PIXELS);
  grids->J_21               = (float *)         fftwf_malloc(sizeof(float)        * HII_TOT_NUM_PIXELS);
  grids->Mvir_crit          = (float *)         fftwf_malloc(sizeof(float)        * HII_TOT_NUM_PIXELS);

  SID_log("Initialising grids...", SID_LOG_COMMENT);

  memset(grids->J_21_at_ionization, 0., sizeof(float)*HII_TOT_NUM_PIXELS);
  memset(grids->J_21, 0., sizeof(float)*HII_TOT_NUM_PIXELS);
  memset(grids->Mvir_crit, 0, sizeof(float) * HII_TOT_NUM_PIXELS);
  for(int ii=0; ii<HII_TOT_NUM_PIXELS; ii++)
  {
    grids->z_at_ionization[ii] = -1;
    grids->xH[ii] = 1.0;
  }

  memset(grids->stars, 0, sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  memset(grids->stars_filtered, 0, sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  memset(grids->deltax, 0, sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);
  memset(grids->deltax_filtered, 0, sizeof(fftw_complex) * HII_KSPACE_NUM_PIXELS);

  SID_log(" ...done", SID_LOG_CLOSE);
}


void free_reionization_grids(run_globals_t *run_globals)
{
  tocf_grids_t *grids = &(run_globals->tocf_grids);

  fftwf_free(grids->Mvir_crit);
  fftwf_free(grids->J_21);
  fftwf_free(grids->J_21_at_ionization);
  fftwf_free(grids->z_at_ionization);
  fftwf_free(grids->deltax_filtered);
  fftwf_free(grids->deltax);
  fftwf_free(grids->stars_filtered);
  fftwf_free(grids->stars);
  fftwf_free(grids->xH);
}


int find_cell(double pos, double box_size)
{
  int HII_dim = tocf_params.HII_dim;
  return (int)round(pos/box_size*(double)HII_dim);
}


void construct_stellar_grids(run_globals_t *run_globals)
{
  galaxy_t *gal;
  int i, j, k;
  double box_size = (double)(run_globals->params.BoxSize);
  double Hubble_h = run_globals->params.Hubble_h;
  double UnitMass_in_g = run_globals->units.UnitMass_in_g;
  double UnitTime_in_s = run_globals->units.UnitTime_in_s;
  float *stellar_grid = (float *)(run_globals->tocf_grids.stars);
  int HII_dim = tocf_params.HII_dim;

  SID_log("Constructing stellar mass grid...", SID_LOG_OPEN);

  // init the grid
  for(int ii=0; ii<HII_TOT_FFT_NUM_PIXELS; ii++)
    *(stellar_grid + ii) = 0.0;

  // Loop through each valid galaxy and add its stellar mass to the appropriate cell
  gal = run_globals->FirstGal;
  while(gal!=NULL)
  {
    if(gal->Type < 3)
    {
      i = find_cell(gal->Pos[0], box_size);
      j = find_cell(gal->Pos[1], box_size);
      k = find_cell(gal->Pos[2], box_size);
      *(stellar_grid + HII_R_FFT_INDEX(i,j,k)) += gal->StellarMass;
    }
    gal = gal->Next;
  }

  // Do one final pass to put the grid in the correct units (Msol or Msol/yr)
  for(int i=0; i<HII_dim; i++)
    for(int j=0; j<HII_dim; j++)
      for(int k=0; k<HII_dim; k++)
        *(stellar_grid + HII_R_FFT_INDEX(i,j,k)) *= 1.e10/Hubble_h;

  SID_log("done", SID_LOG_CLOSE);

}


void save_tocf_grids(run_globals_t *run_globals, hid_t group_id)
{
  tocf_grids_t *grids = &(run_globals->tocf_grids);
  hsize_t dims = HII_TOT_NUM_PIXELS;
  int HII_dim = tocf_params.HII_dim;
  float *grid;

  SID_log("Saving tocf grids...", SID_LOG_OPEN);

  // float grids
  H5LTmake_dataset_float(group_id, "xH", 1, &dims, grids->xH);
  H5LTmake_dataset_float(group_id, "J_21", 1, &dims, grids->J_21);
  H5LTmake_dataset_float(group_id, "J_21_at_ionization", 1, &dims, grids->J_21_at_ionization);
  H5LTmake_dataset_float(group_id, "z_at_ionization", 1, &dims, grids->z_at_ionization);
  H5LTmake_dataset_float(group_id, "Mvir_crit", 1, &dims, grids->Mvir_crit);

  // fftw padded grids
  grid = (float *)SID_calloc(HII_TOT_NUM_PIXELS * sizeof(float));

  for(int ii=0; ii<HII_dim; ii++)
    for(int jj=0; jj<HII_dim; jj++)
      for(int kk=0; kk<HII_dim; kk++)
        grid[HII_R_INDEX(ii,jj,kk)] = *((float *)(grids->stars) + HII_R_FFT_INDEX(ii,jj,kk));
  H5LTmake_dataset_float(group_id, "stars", 1, &dims, grid);

  memset((void *)grid, 0, sizeof(float)*HII_TOT_NUM_PIXELS);
  for(int ii=0; ii<HII_dim; ii++)
    for(int jj=0; jj<HII_dim; jj++)
      for(int kk=0; kk<HII_dim; kk++)
        grid[HII_R_INDEX(ii,jj,kk)] = *((float *)(grids->deltax) + HII_R_FFT_INDEX(ii,jj,kk));
  H5LTmake_dataset_float(group_id, "deltax", 1, &dims, grid);

  SID_free(SID_FARG grid);

  SID_log(" done", SID_LOG_CLOSE);
}


// void assign_ionization_to_halos(run_globals_t *run_globals, halo_t *halo, int n_halos)
// {
//   double box_size = (double)(run_globals->params.BoxSize);
//   int i, j, k;

//   SID_log("Assigning cell ionization values to halos...", SID_LOG_OPEN|SID_LOG_TIMER);

//   for(int i_halo=0; i_halo<n_halos; i_halo++)
//   {
//     i = find_cell(halo[i_halo].Pos[0], box_size);
//     j = find_cell(halo[i_halo].Pos[1], box_size);
//     k = find_cell(halo[i_halo].Pos[2], box_size);
//     halo[i_halo].CellIonization = 1.0 - xH_grid[HII_R_INDEX(i,j,k)];
//   }

//   SID_log("...done", SID_LOG_CLOSE);
// }

#endif
