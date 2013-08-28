#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>

void init_reionization(run_globals_struct *run_globals)
{
  //! Initialize the 21cmfast parameters structure
  tocf_params_struct *params = &(run_globals->tocf_params);

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
}

static int read_xH_grid(run_globals_struct *run_globals, int snapshot, float *xH_grid)
{
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
}

void do_reionization(run_globals_struct *run_globals, int snapshot)
{
  tocf_params_struct *params = &(run_globals->tocf_params);
  bool found_flag = false;
  // float *xH_grid;

  // snapshot -=1;
  params->snapshot = snapshot;

  // for(int ii=0; ii<NOUT; ii++)
  //   if(run_globals->ListOutputSnaps[ii]==(snapshot))
  //   {
  //     found_flag = true;
  //     break;
  //   }

  // if(found_flag == false)
  // {
  //   SID_log("The previous snapshot was not output so we can't do reionization calculation here!", SID_LOG_COMMENT);
  //   ABORT(EXIT_FAILURE);
  // }

  find_HII_bubbles(params);

  // read_xH_grid(run_globals, snapshot, xH_grid);

  // SID_free(SID_FARG xH_grid);
}

