#include "meraxes.h"

void do_reionization(run_globals_struct *run_globals, int snapshot)
{
#ifdef USE_TOCF
  tocf_params_struct *params = &(run_globals->tocf_params);
  bool found_flag = false;

  params->snapshot = snapshot;

  find_HII_bubbles(params);
#endif
}

