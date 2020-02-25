#ifndef EVOLVE_H
#define EVOLVE_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C" {
#endif

int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof);
void passively_evolve_ghost(struct galaxy_t* gal, int snapshot);

#ifdef __cplusplus
}
#endif

#endif
