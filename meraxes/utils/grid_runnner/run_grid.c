#define _MAIN
#include <stdio.h>
#include <meraxes.h>
#include <sys/stat.h>
#include <mpi.h>

typedef struct grid_params_t {
  double SfEfficiency;
  double SnReheatEff;
  double SnEjectionEff;
  double ReincorporationEff;
} grid_params_t;


static void update_params(run_globals_t *run_globals, grid_params_t *grid_params, int i_run)
{
  physics_params_t *params = &(run_globals->params.physics);

  params->SfEfficiency = grid_params[i_run].SfEfficiency;
  params->SnReheatEff = grid_params[i_run].SnReheatEff;
  params->SnEjectionEff = grid_params[i_run].SnEjectionEff;
  params->ReincorporationEff = grid_params[i_run].ReincorporationEff;
}


static int read_grid_params(char *fname, grid_params_t **grid_params)
{

  FILE *fd;
  int n_grid_runs;

  fd = fopen(fname, "r");

  if (fd != NULL)
    fscanf(fd, "%d\n", &n_grid_runs);

  SID_log("Reading %d parameter sets...", SID_LOG_COMMENT, n_grid_runs);

  *grid_params = SID_malloc(sizeof(grid_params_t) * n_grid_runs);

  for (int ii = 0; ii < n_grid_runs; ii++)
    fscanf(fd, "%g %g %g %g\n", (*grid_params)[ii].SfEfficiency,
        (*grid_params)[ii].SnReheatEff, (*grid_params)[ii].SnEjectionEff,
        (*grid_params)[ii].ReincorporationEff);

  fclose(fd);

  return n_grid_runs;

}


int main(int argc, char *argv[])
{

  // init MPI
  MPI_Comm model_comm, world_comm;
  int mpi_world_size;
  MPI_Group model_group;
  int my_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
  MPI_Comm_size(world_comm, &mpi_world_size);
  MPI_Comm_rank(world_comm, &my_rank);
  int analysis_rank = mpi_world_size-1;
  int my_color = (my_rank != analysis_rank) ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(world_comm, my_color, my_rank, &model_comm);

  // init SID
  SID_init(&argc, &argv, NULL, &model_comm);
  if (my_color == 0)
    SID_log_set_fp(stderr);

  struct stat filestatus;
  run_globals_t run_globals;
  grid_params_t *grid_params;
  int n_grid_runs;
  char cmd[STRLEN];
  char file_name_galaxies[STRLEN];

  // deal with any input arguments
  if (argc != 3)
  {
    SID_log("\n  usage: %s <parameterfile> <gridfile>\n\n", SID_LOG_COMMENT, argv[0]);
    ABORT(1);
  }

  // read in the grid parameters
  n_grid_runs = read_grid_params(argv[2], &grid_params);

  if (my_color == 0)
  {
    // read the input parameter file
    read_parameter_file(&run_globals, argv[1], 0);

    // set the interactive flag
    run_globals.params.FlagInteractive = 1;

    // Check to see if the output directory exists and if not, create it
    if (stat(run_globals.params.OutputDir, &filestatus) != 0)
      mkdir(run_globals.params.OutputDir, 02755);

    // initiate meraxes
    init_meraxes(&run_globals);

    // calculate the output hdf5 file properties for later use
    calc_hdf5_props(&run_globals);
  }

  // Run the model!
  for (int ii = 0; ii < n_grid_runs; ii++)
  {
    sprintf(file_name_galaxies, "meraxes_%03d", ii);
    if (my_color == 0)
    {
      strcpy(run_globals.params.FileNameGalaxies, file_name_galaxies);
      update_params(&run_globals, grid_params, ii);
      dracarys(&run_globals);
    }

    // Sync here so that we can ensure we are not going to try and read
    // non-existent files
    MPI_Barrier(world_comm);

    if (my_rank == analysis_rank)
    {
      // N.B. The script called here should delete the output files once it is finished with them
      sprintf(cmd, "/home/smutch/pyenv/bin/python /home/smutch/models/21cm_sam/meraxes/utils/grid_runnner/analyse_run.py %s/%s.hdf5",
          run_globals.params.OutputDir, file_name_galaxies);
      system(cmd);
    }
  }

  // cleanup
  if (my_color == 0)
    cleanup(&run_globals);

  SID_free(SID_FARG grid_params);
  MPI_Comm_free(&world_comm);
  SID_exit(EXIT_SUCCESS);
}
