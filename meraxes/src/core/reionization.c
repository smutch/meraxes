#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>

void init_reionization(run_globals_struct *run_globals)
{
#ifdef USE_TOCF
  //! Initialize the 21cmfast parameters structure
  tocf_params_struct *params = &(run_globals->tocf_params);

  // N.B. Currently BOX_LEN and HII_DIM must be set in the INIT_PARAMS.H file
  // of 21cmFAST before it's compilation.  All cosmological parameters must
  // also be set in COSMOLOGY.H.  This should be fixed in future.

  params->meraxes_fname  = run_globals->FNameOut;
  params->logfile_dir    = run_globals->params.TOCF_LogFileDir;
  params->snapshot       = -1;  // Will be updated later
  params->num_th         = run_globals->params.TOCF_NThreads;
  params->ion_eff_factor = -1;  // default
  params->tvir_min       = -1;  // default
  params->mean_free_path = -1;  // default
  params->zlist          = run_globals->ZZ;
  params->sim_dir        = run_globals->params.SimulationDir;
  params->sim_name       = run_globals->params.SimName;
#else
  SID_log_error("TOCF_Flag = 1, but Meraxes has not been compiled with 21cmFAST...", SID_LOG_COMMENT);
  ABORT(EXIT_FAILURE);
#endif
}

int malloc_xH_grid(run_globals_struct *run_globals, int snapshot, float **xH_grid)
{

  herr_t status;
  hid_t file_id;
  hid_t group_id;
  hsize_t dims[1];
  char group_name[128];
  int xH_grid_dim;

  file_id = H5Fopen(run_globals->FNameOut, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(group_name, "Snap%03d", snapshot);
  group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  
  status = H5LTget_dataset_info(group_id, "xH_grid", dims, NULL, NULL);
  if (status<0)
  {
    SID_log("Failed to get xH_grid dimensions.  Aborting...", SID_LOG_COMMENT, snapshot);
    ABORT(EXIT_FAILURE);
  }
  *xH_grid = (float *)SID_malloc(sizeof(float) * (size_t)dims[0]);

  status = H5LTget_attribute_int(group_id, "xH_grid", "HII_dim", &xH_grid_dim); 

  H5Gclose(group_id);
  H5Fclose(file_id);

  return xH_grid_dim;
}

void read_xH_grid(run_globals_struct *run_globals, int snapshot, float *xH_grid)
{
#ifdef USE_TOCF
  hid_t file_id;
  hid_t group_id;
  char group_name[128];

  file_id = H5Fopen(run_globals->FNameOut, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(group_name, "Snap%03d", snapshot);
  group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  
  H5LTread_dataset_float(group_id, "xH_grid", xH_grid);

  H5Gclose(group_id);
  H5Fclose(file_id);

#else
  return -1;
#endif
}

static inline int find_cell(double pos, int xH_dim, double box_size)
{
  return (int)((pos/box_size)*(double)xH_dim);
}

static inline int xH_grid_index(int i, int j, int k, int xH_dim)
{
  return k+xH_dim*(j+xH_dim*i);
}

void assign_ionization_to_galaxies(run_globals_struct *run_globals, float *xH_grid, int xH_dim)
{

  galaxy_struct *gal = run_globals->FirstGal;
  double box_size = (double)(run_globals->params.BoxSize);
  int i, j, k;

  SID_log("Assigning cell ionization values to galaxies...", SID_LOG_OPEN|SID_LOG_TIMER);
  SID_log("xH_dim = %d", SID_LOG_COMMENT, xH_dim);
  
  while(gal != NULL)
  {
    i = find_cell(gal->Pos[0], xH_dim, box_size);
    j = find_cell(gal->Pos[1], xH_dim, box_size);
    k = find_cell(gal->Pos[2], xH_dim, box_size);
    gal->CellIonization = 1.0 - xH_grid[xH_grid_index(i,j,k, xH_dim)];
    gal = gal->Next;
  }

  SID_log("...done", SID_LOG_CLOSE);

}
