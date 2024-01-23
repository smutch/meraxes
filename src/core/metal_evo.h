#ifndef METAL_EVO_H
#define METAL_EVO_H

#if USE_MINI_HALOS
#include <fftw3-mpi.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

  void assign_slabs_metals(void);
  void init_metal_grids(void);
  void malloc_metal_grids(void);
  void free_metal_grids(void);
  int map_galaxies_to_slabs_metals(int ngals);
  void construct_metal_grids(int snapshot, int local_ngals);
  void save_metal_input_grids(int snapshot);
  void gen_metal_grids_fname(const int snapshot, char* name, const bool relative);
  void assign_probability_to_galaxies(int ngals_in_metal_slabs, int snapshot, int flag_property);

#ifdef __cplusplus
}
#endif

#endif
#endif