#ifndef EVOLVE_H
#define EVOLVE_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

#if USE_MINI_HALOS
  int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof, int* gal_counter_Pop3, int* gal_counter_Pop2, int* gal_counter_enriched);
#else
  int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof);
#endif
  void passively_evolve_ghost(struct galaxy_t* gal, int snapshot);

#ifdef __cplusplus
}
#endif

#endif
