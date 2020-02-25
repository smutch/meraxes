#ifndef INFALL_H
#define INFALL_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C" {
#endif

double gas_infall(fof_group_t* FOFgroup, int snapshot);
void add_infall_to_hot(struct galaxy_t* central, double infall_mass);

#ifdef __cplusplus
}
#endif

#endif
