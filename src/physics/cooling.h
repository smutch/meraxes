#ifndef COOLING_H
#define COOLING_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  double gas_cooling(struct galaxy_t* gal);
  void cool_gas_onto_galaxy(struct galaxy_t* gal, double cooling_mass);

#ifdef __cplusplus
}
#endif

#endif
