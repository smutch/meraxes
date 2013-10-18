#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>

int malloc_xH_grid(run_globals_t *run_globals, int snapshot, float **xH_grid)
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

void read_xH_grid(run_globals_t *run_globals, int snapshot, float *xH_grid)
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

void assign_ionization_to_halos(run_globals_t *run_globals, halo_t *halo, int n_halos, float *xH_grid, int xH_dim)
{

  double box_size = (double)(run_globals->params.BoxSize);
  int i, j, k;

  SID_log("Assigning cell ionization values to halos...", SID_LOG_OPEN|SID_LOG_TIMER);
  SID_log("xH_dim = %d", SID_LOG_COMMENT, xH_dim);

  for(int i_halo=0; i_halo<n_halos; i_halo++)
  {
    i = find_cell(halo[i_halo].Pos[0], xH_dim, box_size);
    j = find_cell(halo[i_halo].Pos[1], xH_dim, box_size);
    k = find_cell(halo[i_halo].Pos[2], xH_dim, box_size);
    halo[i_halo].CellIonization = 1.0 - xH_grid[xH_grid_index(i,j,k, xH_dim)];
  }

  SID_log("...done", SID_LOG_CLOSE);

}
