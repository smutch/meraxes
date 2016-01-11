#include "meraxes.h"
#include <math.h>
#include <assert.h>

void myexit(int signum)
{
  fprintf(stderr, "Task: %d\tnode: %s\tis exiting.\n\n\n", SID.My_rank, SID.My_node);
  cn_quote();
  // if(SID.n_proc > 1)
  // {
  //   fflush(SID.fp_log);
  //   fclose(SID.fp_log);
  // }
  mpi_debug_here();
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


int compare_floats(const void *a, const void *b) {
  float value = *(float *)a - *(float *)b;
  if (value > 0)
    return 1;
  else if (value < 0)
    return -1;
  else
    return 0;
}

static float inline apply_pbc(float delta)
{
  float box_size = (float)(run_globals.params.BoxSize);

  if(fabs(delta-box_size)<fabs(delta))
    delta -= box_size;
  if(fabs(delta+box_size)<fabs(delta))
    delta += box_size;

  return delta;
}


float comoving_distance(float a[3], float b[3])
{
  float dx = apply_pbc(a[0] - b[0]);
  float dy = apply_pbc(a[1] - b[1]);
  float dz = apply_pbc(a[2] - b[2]);

  float dist = sqrtf(dx*dx + dy*dy + dz*dz);
  assert(dist <= (sqrtf(3.0)/2.0 * run_globals.params.BoxSize));

  return dist;
}


double accurate_sumf(float *arr, int n) {
  // inplace reorder and sum
  qsort(arr, n, sizeof(float), compare_floats);

  double total = 0;
  for(int ii=0; ii<n; ii++)
    total += arr[ii];

  return total;
}


int grid_index(int i, int j, int k, int dim, int type)
{
  int ind;
  switch(type)
  {
    case INDEX_PADDED:
      ind  = k + (2*(dim/2 +1)) * (j + dim * i);
      break;
    case INDEX_REAL:
      ind = k + dim * (j + dim *i);
      break;
    case INDEX_COMPLEX_HERM:
      ind = k + (dim/2 + 1) * (j + dim * i);
      break;
  }

  return ind;
}

