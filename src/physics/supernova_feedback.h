#ifndef SUPERNOVA_FEEDBACK_H
#define SUPERNOVA_FEEDBACK_H

#include "meraxes.h"

#define EnergySN 1e51 // erg, maybe you can read it from tables
#define N_SN_Pop2 0.01 // number of SN x solar mass, you might add this as an input parameter
#define N_SN_Pop3 0.1 // So far you don't have separation between Pop 2 and Pop 3 but soon...

#ifdef __cplusplus
extern "C"
{
#endif

  void update_reservoirs_from_sn_feedback(struct galaxy_t* gal,
                                          double m_reheat,
                                          double m_eject,
                                          double m_recycled,
                                          double new_metals);
  void delayed_supernova_feedback(struct galaxy_t* gal, int snapshot);
  void calc_metal_bubble(struct galaxy_t* gal, int snapshot);
  void contemporaneous_supernova_feedback(struct galaxy_t* gal,
                                          double* m_stars,
                                          int snapshot,
                                          double* m_reheat,
                                          double* m_eject,
                                          double* m_recycled,
                                          double* new_metals);

#ifdef __cplusplus
}
#endif

#endif
