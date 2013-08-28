#include <stdio.h>
#include <stdlib.h>
#include "parameter_files/INIT_PARAMS.H"
#include "parameter_files/ANAL_PARAMS.H"

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
