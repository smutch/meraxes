#define _MAIN
#include "meraxes.h"
#include <sys/stat.h>
#include <fenv.h>

int main(int argc, char **argv)
{
  // init SID
  SID_init(&argc, &argv, NULL, NULL);

  struct stat filestatus;

  // deal with any input arguments
  if (argc != 2)
  {
    SID_log("\n  usage: %s <parameterfile>\n\n", SID_LOG_COMMENT, argv[0]);
    ABORT(EXIT_FAILURE);
  }

  // set the rounding mode
  if (!fesetround(1))  // nearest number
  {
    SID_log_error("Failed to set rounding mode!");
    ABORT(EXIT_FAILURE);
  }

#ifdef LOGGER
  dzlog_init("zlog.conf", "default");
  dzlog_notice("Log for rank %d.", SID.My_rank);
#endif

  // read the input parameter file
  read_parameter_file(argv[1], 0);

  // Check to see if the output directory exists and if not, create it
  if (stat(run_globals.params.OutputDir, &filestatus) != 0)
    mkdir(run_globals.params.OutputDir, 02755);

  // initiate meraxes
  init_meraxes();

  // Run the model!
  if ((!run_globals.params.FlagInteractive) & (!run_globals.params.FlagMCMC))
    dracarys();
  else
    while (run_globals.params.FlagInteractive)
    {
      dracarys();
      continue_prompt(argv[1]);
    }

  // cleanup
  cleanup();

#ifdef LOGGER
  zlog_fini();
#endif

  SID_exit(EXIT_SUCCESS);
}