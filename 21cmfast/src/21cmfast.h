#ifndef _21CMFAST
#define _21CMFAST

/*
 * Structures
 */

typedef struct tocf_params_struct tocf_params_struct;
struct tocf_params_struct{
  char   *meraxes_fname;
  char   *logfile_dir;
  int     snapshot;
  int     num_th;
  float   ion_eff_factor;
  float   tvir_min;
  float   mean_free_path;
  double *zlist;
  char   *sim_dir;
  char   *sim_name;
};

/*
 * Functions
 */

int find_HII_bubbles(tocf_params_struct *params);

#endif
