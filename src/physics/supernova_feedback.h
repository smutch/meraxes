#ifndef SUPERNOVA_FEEDBACK_H
#define SUPERNOVA_FEEDBACK_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C" {
#endif

void update_reservoirs_from_sn_feedback(struct galaxy_t* gal, double m_reheat, double m_eject, double m_recycled, double new_metals);
void delayed_supernova_feedback(struct galaxy_t* gal, int snapshot);
void contemporaneous_supernova_feedback(struct galaxy_t* gal, double* m_stars, int snapshot, double* m_reheat, double* m_eject, double* m_recycled, double* new_metals);

#ifdef __cplusplus
}
#endif

#endif
