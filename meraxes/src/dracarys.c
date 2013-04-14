#include "meraxes.h"

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{
  trees_header_struct  trees_header;
  halo_struct         *halo;
  int                  snapshot;
 
  run_params_struct params = run_globals->params;

  // Simple debug test - read in the first snapshot
  snapshot = 0;
  trees_header = read_halos(run_globals, snapshot, &halo);

  int i_halo = 0;
  printf("n_subgroups = %d\n"    , trees_header.n_subgroups);
  printf("halo[%d].id = %d\n"    , i_halo, halo[i_halo].id);
  printf("halo[%d].type = %d\n"  , i_halo, halo[i_halo].type);
  printf("halo[%d].Mvir = %.2e\n", i_halo, halo[i_halo].Mvir);
  printf("halo[%d].Rvir = %.2e\n", i_halo, halo[i_halo].Rvir);

  printf("ListOutputSnaps[0] = %d\n", run_globals->ListOutputSnaps[0]);

  free_halos(&halo);

}
