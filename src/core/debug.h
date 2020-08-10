#ifndef DEBUG_H
#define DEBUG_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  void mpi_debug_here(void);
  void check_mhysa_pointer(void);
  void check_counts(fof_group_t* fof_group, int NGal, int NFof);
  void write_single_grid(const char* fname,
                         float* grid,
                         int local_ix_start,
                         int local_nix,
                         int dim,
                         const char* grid_name,
                         bool padded_flag,
                         bool create_file_flag);

#ifdef DEBUG
  void check_pointers(halo_t* halos, fof_group_t* fof_groups, trees_info_t* trees_info);
#endif

#ifdef __cplusplus
}
#endif

#endif
