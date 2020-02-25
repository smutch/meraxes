#ifndef PHYS_REIONIZATION_H
#define PHYS_REIONIZATION_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C" {
#endif

void calculate_Mvir_crit(double redshift);
double tocf_modifier(struct galaxy_t* gal, double Mvir);
double reionization_modifier(struct galaxy_t* gal, double Mvir, int snapshot);
double sobacchi2013_modifier(double Mvir, double redshift);
double gnedin2000_modifer(double Mvir, double redshift);

#ifdef __cplusplus
}
#endif

#endif
