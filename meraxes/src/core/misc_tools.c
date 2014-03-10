#include "meraxes.h"

double calc_metallicity(double total_gas, double metals)
{
  double Z;

  if((total_gas > 0) && (metals > 0))
    Z = metals / total_gas;
  else
    Z = 0.0;

  if(Z < 0)
    Z = 0.0;
  if(Z > 1)
    Z= 1.0;

  return Z;
}


int compare_ints(const void *a, const void *b)
{
  return *((int *)a) - *((int *)b);
}

