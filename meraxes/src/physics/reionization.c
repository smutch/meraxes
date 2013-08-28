#include "meraxes.h"

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


void do_reionization(run_globals_struct *run_globals, int i_out)
{
  tocf_params_struct *params = &(run_globals->tocf_params);

  params->snapshot = run_globals->ListOutputSnaps[i_out];

  find_HII_bubbles(params);
}

