#ifndef MAGNITUDES_H
#define MAGNITUDES_H

#include "meraxes.h"

#ifdef CALC_MAGS
#include <sector.h>

#define TOL 1e-30 // Minimum Flux

enum core
{
  MASTER
};

#ifdef __cplusplus
extern "C"
{
#endif

  void init_luminosities(struct galaxy_t* gal);
  void add_luminosities(mag_params_t* miniSpectra,
                        struct galaxy_t* gal,
                        int snapshot,
                        double metals,
                        double sfr,
                        double new_stars);
  void merge_luminosities(struct galaxy_t* target, struct galaxy_t* gal);
  void init_templates_mini(mag_params_t* miniSpectra,
                           char* fName,
                           char* fNameIII,
                           double* LTTime,
                           int* targetSnaps,
                           double* redshifts,
                           double* betaBands,
                           int nBeta,
                           double* restBands,
                           int nRest,
                           double tBC);
  void init_magnitudes(void);
  void cleanup_mags(void);

#ifdef __cplusplus
}
#endif

#endif
#endif
