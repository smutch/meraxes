#include "meraxes.h"

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{
 
  run_params_struct params = run_globals->params;

  // As a simple debug test - just print some of the input parameters
  printf("DEBUG TEST:\n");
  printf("physics function peak = %.2f\n",  params.physics.peak);
  printf("physics function sigma = %.2f\n", params.physics.sigma);
  printf("physics function frac = %.2f\n",  params.physics.stellarfrac);
  printf("BoxSize = %.2f\n", params.BoxSize);
  printf("SimHubble_h = %.2f\n", params.SimHubble_h);
  printf("---------\n");

}
