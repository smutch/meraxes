#ifndef READ_HALOS_H
#define READ_HALOS_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  trees_info_t read_halos(int snapshot,
                          halo_t** halo,
                          fof_group_t** fof_group,
                          int** index_lookup,
                          trees_info_t* snapshot_trees_info);
  void initialize_halo_storage(void);
  void free_halo_storage(void);

  trees_info_t read_trees_info__gbptrees(int snapshot);
  void read_trees__gbptrees(int snapshot,
                            halo_t* halo,
                            int n_halos,
                            fof_group_t* fof_group,
                            int n_fof_groups,
                            int n_requested_forests,
                            int* n_halos_kept,
                            int* n_fof_groups_kept,
                            int* index_lookup);

  void read_trees__velociraptor(int snapshot,
                                halo_t* halos,
                                int* n_halos,
                                fof_group_t* fof_groups,
                                int* n_fof_groups,
                                int* index_lookup);
  trees_info_t read_trees_info__velociraptor(const int snapshot);

#ifdef __cplusplus
}
#endif

#endif
