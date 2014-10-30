#define _MAIN
#include <stdio.h>
#include <string.h>
#include <meraxes.h>
#include <sys/stat.h>
#include <mpi.h>
#include <gbpLib.h>

/*
 * Header info
 */

typedef struct grid_params_t {
  double SfEfficiency;
  double SnReheatEff;
  double SnEjectionEff;
  double ReincorporationEff;
  double MergerTimeFactor;
} grid_params_t;

MPI_Comm world_comm;
int world_rank;


/*
 * Code
 */

static void update_params(run_globals_t *run_globals, grid_params_t *grid_params, int i_run)
{
  physics_params_t *params = &(run_globals->params.physics);

  params->SfEfficiency = grid_params[i_run].SfEfficiency;
  params->SnReheatEff = grid_params[i_run].SnReheatEff;
  params->SnEjectionEff = grid_params[i_run].SnEjectionEff;
  params->ReincorporationEff = grid_params[i_run].ReincorporationEff;
  params->MergerTimeFactor = grid_params[i_run].MergerTimeFactor;
}

static int read_grid_params(char *fname, grid_params_t **grid_params)
{
  int n_grid_runs;

  if (world_rank == 0)
  {
    FILE *fd;
    int ret;

    fd = fopen(fname, "r");

    if (fd != NULL)
      fscanf(fd, " %d ", &n_grid_runs);
    else
    {
      fprintf(stderr, "Failed to open %s\n", fname);
      ABORT(EXIT_FAILURE);
    }

    SID_log("Reading %d parameter sets...", SID_LOG_COMMENT, n_grid_runs);

    *grid_params = malloc(sizeof(grid_params_t) * n_grid_runs);

    for (int ii = 0; ii < n_grid_runs; ii++)
    {
      ret = fscanf(fd, " %lg %lg %lg %lg %lg ", &((*grid_params)[ii].SfEfficiency),
          &((*grid_params)[ii].SnReheatEff), &((*grid_params)[ii].SnEjectionEff),
          &((*grid_params)[ii].ReincorporationEff), &((*grid_params)[ii].MergerTimeFactor));
      if (ret != 5)
      {
        fprintf(stderr, "Failed to read %s correctly\n", fname);
        ABORT(EXIT_FAILURE);
      }
    }

    fclose(fd);
  }

  MPI_Bcast(&n_grid_runs, 1, MPI_INT, 0, world_comm);
  if (world_rank > 0)
    *grid_params = malloc(sizeof(grid_params_t) * n_grid_runs);
  MPI_Bcast(*grid_params, sizeof(grid_params_t)*n_grid_runs, MPI_BYTE, 0, world_comm);

  return n_grid_runs;
}


int main(int argc, char *argv[])
{

  // deal with any input arguments
  if (argc != 3)
  {
    fprintf(stderr, "\n  usage: %s <parameterfile> <gridfile>\n\n", argv[0]);
    ABORT(EXIT_SUCCESS);
  }

  // init MPI
  MPI_Comm model_comm;
  int world_size;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
  MPI_Comm_size(world_comm, &world_size);
  MPI_Comm_rank(world_comm, &world_rank);
  bool analysis_rank = (world_rank == world_size-1);
  int my_color = analysis_rank ? MPI_UNDEFINED : 0;
  MPI_Comm_split(world_comm, my_color, world_rank, &model_comm);

  // init SID
  if (!analysis_rank)
  {
    SID_init(&argc, &argv, NULL, &model_comm);
    SID_log_set_fp(stderr);
  }

  struct stat filestatus;
  run_globals_t run_globals;
  grid_params_t *grid_params;
  int n_grid_runs;
  char cmd[STRLEN];
  char file_name_galaxies[STRLEN];
  char output_dir[STRLEN];

  // read in the grid parameters
  n_grid_runs = read_grid_params(argv[2], &grid_params);

  if (!analysis_rank)
  {

  // debugging
#ifdef DEBUG
    if (SID.My_rank == 0)
      mpi_debug_here();
#endif

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

  // Copy the output dir of the model to the analysis rank
  if (world_rank == 0)
    MPI_Send(run_globals.params.OutputDir, STRLEN, MPI_CHAR, world_size-1, 16, world_comm);
  if (analysis_rank)
    MPI_Recv(output_dir, STRLEN, MPI_CHAR, 0, 16, world_comm, MPI_STATUS_IGNORE);

  // Run the model!
  for (int ii = 0; ii < n_grid_runs; ii++)
  {
    sprintf(file_name_galaxies, "meraxes_%03d", ii);
    if (!analysis_rank)
    {
      strcpy(run_globals.params.FileNameGalaxies, file_name_galaxies);
      SID_log(">>>> New FileNameGalaxies is: %s", SID_LOG_COMMENT, run_globals.params.FileNameGalaxies);
      update_params(&run_globals, grid_params, ii);
      SID_log(">>>> Updated params to: %g %g %g %g %g", SID_LOG_COMMENT,
          grid_params[ii].SfEfficiency, grid_params[ii].SnReheatEff,
          grid_params[ii].SnEjectionEff, grid_params[ii].ReincorporationEff,
          grid_params[ii].MergerTimeFactor);

      dracarys(&run_globals);
    }

    // Sync here so that we can ensure we are not going to try and read
    // non-existent files
    MPI_Barrier(world_comm);

    if (analysis_rank)
    {
      // N.B. The script called here should delete the output files once it is finished with them
      sprintf(cmd, "/home/smutch/pyenv/bin/python analyse_run.py %s/%s.hdf5",
          output_dir, file_name_galaxies);
      printf(" >>>> ANALYSIS : Calling\n\t%s\n", cmd);
      system(cmd);
    }
  }

  // cleanup
  if (!analysis_rank)
    cleanup(&run_globals);

  free(grid_params);
  MPI_Comm_free(&world_comm);

  if (!analysis_rank)
    SID_exit(EXIT_SUCCESS);

  // Only the analysis rank will make it here
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}
