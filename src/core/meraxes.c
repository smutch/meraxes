#define _MAIN
#include "meraxes.h"
#include <sys/stat.h>
#include <fenv.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &run_globals.mpi_comm);
  MPI_Comm_rank(MPI_COMM_WORLD, &run_globals.mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &run_globals.mpi_size);

  // init mlog
  init_mlog(MPI_COMM_WORLD, stdout, stdout, stderr);

  struct stat filestatus;

  // deal with any input arguments
  if (argc != 2)
  {
    mlog("\n  usage: %s <parameterfile>\n\n", MLOG_MESG, argv[0]);
    ABORT(EXIT_FAILURE);
  }

  // set the rounding mode
  if (!fesetround(1))  // nearest number
  {
    mlog_error("Failed to set rounding mode!");
    ABORT(EXIT_FAILURE);
  }

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

  MPI_Finalize();
  return EXIT_SUCCESS;
}
