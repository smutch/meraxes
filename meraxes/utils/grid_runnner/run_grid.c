#define _MAIN
#include <meraxes.h>
#include <sys/stat.h>

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
    // init SID
  SID_init(&argc, &argv, NULL);

  struct stat filestatus;
  run_globals_t run_globals;
  grid_params_t *grid_params;
  int n_grid_runs;
  char cmd[STRLEN];

  int real_n_proc = SID.n_proc;
  int real_last_rank = real_n_proc - 1;
  SID.n_proc -= 1;

  // deal with any input arguments
  if (argc != 3)
  {
    SID_log("\n  usage: %s <parameterfile> <gridfile>\n\n", SID_LOG_COMMENT, argv[0]);
    ABORT(1);
  }

  // read in the grid parameters
  n_grid_runs = read_grid_params(argv[2], &grid_params);

  // read the input parameter file
  read_parameter_file(&run_globals, argv[1], 0);

  // set the interactive flag
  run_globals.params.FlagInteractive = 1;

  if (SID.My_rank < real_last_rank)
  {
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
    sprintf(run_globals.params.FileNameGalaxies, "meraxes_%03d", ii);
    if (SID.My_rank < real_last_rank)
    {
      update_params(&run_globals, grid_params, ii);
      dracarys(&run_globals);
    }
    else
    {
      // TODO: Complete this...
      sprintf(cmd, "/home/smutch/pyenv/bin/python /home/smutch/models/21cm_sam/meraxes/output/results/smf_z5.py %s/%s.hdf5",
          run_globals.params.OutputDir, run_globals.params.FileNameGalaxies);
      system(cmd);
    }
  }

  // cleanup
  SID_free(SID_FARG grid_params);
  cleanup(&run_globals);

  SID_exit(EXIT_SUCCESS);

  return 0;
}
