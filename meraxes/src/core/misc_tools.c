#include "meraxes.h"
#include <math.h>

void myexit(int signum)
{
  fprintf(stderr, "Task: %d\tnode: %s\tis exiting.\n\n\n", SID.My_rank, SID.My_node);
  cn_quote();
  // if(SID.n_proc > 1)
  // {
  //   fflush(SID.fp_log);
  //   fclose(SID.fp_log);
  // }
  SID_exit(signum);
}


double calc_metallicity(double total_gas, double metals)
{
  double Z;

  if ((total_gas > 0) && (metals > 0))
    Z = metals / total_gas;
  else
    Z = 0.0;

  if (Z < 0)
    Z = 0.0;
  if (Z > 1)
    Z = 1.0;

  return Z;
}


int compare_ints(const void *a, const void *b)
{
  return *((int*)a) - *((int*)b);
}


static float inline apply_pbc(run_globals_t *run_globals, float delta)
{
  float box_size = (float)(run_globals->params.BoxSize);

  if(fabs(delta-box_size)<fabs(delta))
    delta -= box_size;
  if(fabs(delta+box_size)<fabs(delta))
    delta += box_size;

  return delta;
}


float distance(run_globals_t *run_globals, float a[3], float b[3])
{
  float dx = apply_pbc(run_globals, a[0] - b[0]);
  float dy = apply_pbc(run_globals, a[1] - b[1]);
  float dz = apply_pbc(run_globals, a[2] - b[2]);

  return sqrtf(dx*dx + dy*dy + dz*dz);
}

