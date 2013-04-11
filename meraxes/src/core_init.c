#include "meraxes.h"
#include <time.h>

void init_meraxis(run_globals_struct *run_globals)
{
  int i;

  run_globals->random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(run_globals->random_generator, 42);	 // start-up seed 

  // set_units();
  srand((unsigned) time(NULL));

  // read_output_snaps();
  // read_snap_list();

  // for(i = 0; i < (run_globals->params).SnaplistLength; i++)
  // {
  //   (run_globals->ZZ)[i] = 1 / (run_globals->AA)[i] - 1;
  //   (run_globals->Age)[i] = time_to_present((run_globals->ZZ)[i]);
  // }

}
