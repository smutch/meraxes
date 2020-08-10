#ifndef GALAXIES_H
#define GALAXIES_H

#include "meraxes.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

  struct galaxy_t* new_galaxy(int snapshot, unsigned long halo_ID);
  void copy_halo_props_to_galaxy(struct halo_t* halo, struct galaxy_t* gal);
  void reset_galaxy_properties(struct galaxy_t* gal, int snapshot);
  void connect_galaxy_and_halo(struct galaxy_t* gal, struct halo_t* halo, int* merger_counter);
  void create_new_galaxy(int snapshot, struct halo_t* halo, int* NGal, int* new_gal_counter, int* merger_counter);
  void kill_galaxy(struct galaxy_t* gal, struct galaxy_t* prev_gal, int* NGal, int* kill_counter);

#ifdef __cplusplus
}
#endif

#endif
