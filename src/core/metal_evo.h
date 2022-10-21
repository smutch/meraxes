#include <fftw3-mpi.h>
#include <stdbool.h>

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  void assign_slabs_metals(void);
  void init_metal_grids(void);
  void malloc_metal_grids(void);
  void free_metal_grids(void);
  void map_galaxies_to_slab_metals(int ngals);
  void construct_metal_grids(int snapshot, int local_ngals);
  void save_metal_input_grids(int snapshot);
  void save_metal_output_grids(int snapshot);
  void gen_metal_grids_fname(const int snapshot, char* name, const bool relative);
  

#ifdef __cplusplus
}
#endif

#endif