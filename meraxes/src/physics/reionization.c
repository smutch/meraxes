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

static int read_xH_grid(run_globals_struct *run_globals, int snapshot, float *xH_grid)
{
#ifdef USE_TOCF
  herr_t status;
  hid_t file_id;
  hid_t group_id;
  hsize_t dims[1];
  char group_name[128];
  H5T_class_t class_id;
  size_t type_size;

  file_id = H5Fopen(run_globals->FNameOut, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(group_name, "Snap%03d", snapshot);
  group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  
  status = H5LTget_dataset_info(group_id, "xH_grid", dims, &class_id, &type_size);
  xH_grid = (float *)SID_malloc(sizeof(float) * (size_t)dims[0]);

  status = H5LTread_dataset_float(group_id, "xH_grid", xH_grid);

  H5Gclose(group_id);
  H5Fclose(file_id);

  return (int)(dims[0]);
#endif
}

void do_reionization(run_globals_struct *run_globals, int snapshot)
{
#ifdef USE_TOCF
  tocf_params_struct *params = &(run_globals->tocf_params);
  bool found_flag = false;

  params->snapshot = snapshot;

  find_HII_bubbles(params);
#endif
}

