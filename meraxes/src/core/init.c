#include "meraxes.h"
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#ifdef USE_TOCF
#include <21cmfast.h>
#endif

static void read_requested_forest_ids(run_globals_t *run_globals)
{
  if (strlen(run_globals->params.ForestIDFile) == 0)
  {
    run_globals->NRequestedForests = -1;
    run_globals->RequestedForestId = NULL;
    return;
  }

  if (SID.My_rank == 0)
  {
    FILE *fin;
    char *line = NULL;
    size_t len;
    int n_forests = -1;
    int *ids;

    if (!(fin = fopen(run_globals->params.ForestIDFile, "r")))
    {
      SID_log_error("Failed to open file: %s", run_globals->params.ForestIDFile);
      ABORT(EXIT_FAILURE);
    }

    getline(&line, &len, fin);
    n_forests                      = atoi(line);
    run_globals->NRequestedForests = n_forests;

    run_globals->RequestedForestId = SID_malloc(sizeof(int) * n_forests);
    ids                            = run_globals->RequestedForestId;

    for (int ii = 0; ii < n_forests; ii++)
    {
      getline(&line, &len, fin);
      ids[ii] = atoi(line);
    }

    free(line);

    SID_log("Found %d requested halo IDs", SID_LOG_COMMENT, n_forests);

    fclose(fin);
  }


  // broadcast the data to all other ranks
  SID_Bcast(&(run_globals->NRequestedForests), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    run_globals->RequestedForestId = SID_malloc(sizeof(int) * run_globals->NRequestedForests);
  SID_Bcast(run_globals->RequestedForestId, sizeof(int) * run_globals->NRequestedForests, 0, SID.COMM_WORLD);
}

static void read_multiple_runs_params(run_globals_t *run_globals)
{
  int n_params = 6;

  if (strlen(run_globals->params.MultipleRunsParamFile) == 0)
  {
    run_globals->NRuns              = 1;
    run_globals->MultipleRunsParams = NULL;
    return;
  }

  if (SID.My_rank == 0)
  {
    FILE *fin;
    char *line = NULL;
    size_t len;
    int n_runs = -1;
    double **params;

    if (!(fin = fopen(run_globals->params.MultipleRunsParamFile, "r")))
    {
      SID_log_error("Failed to open file: %s", run_globals->params.MultipleRunsParamFile);
      ABORT(EXIT_FAILURE);
    }

    getline(&line, &len, fin);
    n_runs             = atoi(line);
    run_globals->NRuns = n_runs;

    run_globals->MultipleRunsParams = SID_malloc(sizeof(double *) * n_runs);
    params                          = run_globals->MultipleRunsParams;
    for (int ii = 0; ii < n_runs; ii++)
      params[ii] = SID_malloc(sizeof(double) * n_params);

    for (int ii = 0; ii < n_runs; ii++)
    {
      getline(&line, &len, fin);
      sscanf(line, "%le %le %le %le %le %le",
             &(params[ii][0]),
             &(params[ii][1]),
             &(params[ii][2]),
             &(params[ii][3]),
             &(params[ii][4]),
             &(params[ii][5]));
    }

    free(line);

    SID_log("Found %d parameter value sets to run.", SID_LOG_COMMENT, n_runs);

    fclose(fin);
  }


  // broadcast the data to all other ranks
  SID_Bcast(&(run_globals->NRuns), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    run_globals->MultipleRunsParams = SID_malloc(sizeof(double *) * run_globals->NRuns);

  for (int ii = 0; ii < run_globals->NRuns; ii++)
  {
    if (SID.My_rank > 0)
      (run_globals->MultipleRunsParams)[ii] = SID_malloc(sizeof(double) * n_params);

    SID_Bcast((run_globals->MultipleRunsParams)[ii], sizeof(double) * n_params, 0, SID.COMM_WORLD);
  }
}


static void read_snap_list(run_globals_t *run_globals)
{
  if (SID.My_rank == 0)
  {
    FILE *fin;
    int snaplist_len;
    double dummy;
    char fname[STRLEN];
    run_params_t params = run_globals->params;

    sprintf(fname, "%s/a_list.txt", params.SimulationDir);

    if (!(fin = fopen(fname, "r")))
    {
      SID_log_error("failed to read snaplist in file '%s'", fname);
      ABORT(EXIT_FAILURE);
    }

    // Count the number of snapshot list entries
    snaplist_len = 0;
    do
    {
      if (fscanf(fin, " %lg ", &dummy) == 1)
        snaplist_len++;
      else
        break;
    } while (true);

    run_globals->params.SnaplistLength = snaplist_len;
    if (SID.My_rank == 0)
    {
      printf("NOUT = %d\n\n", NOUT);
      printf("found %d defined times in snaplist.\n", snaplist_len);
    }

    // malloc the relevant arrays
    run_globals->AA     = SID_malloc(sizeof(double) * snaplist_len);
    run_globals->ZZ     = SID_malloc(sizeof(double) * snaplist_len);
    run_globals->LTTime = SID_malloc(sizeof(double) * snaplist_len);

    // seek back to the start of the file
    rewind(fin);

    // actually read in the expansion factors
    snaplist_len = 0;
    do
    {
      if (fscanf(fin, " %lg ", &(run_globals->AA[snaplist_len])) == 1)
        snaplist_len++;
      else
        break;
    } while (true);

    // close the file
    fclose(fin);
  }

  // broadcast the read to all other ranks and malloc the necessary arrays
  SID_Bcast(&(run_globals->params.SnaplistLength), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
  {
    run_globals->AA     = SID_malloc(sizeof(double) * run_globals->params.SnaplistLength);
    run_globals->ZZ     = SID_malloc(sizeof(double) * run_globals->params.SnaplistLength);
    run_globals->LTTime = SID_malloc(sizeof(double) * run_globals->params.SnaplistLength);
  }
  SID_Bcast(run_globals->AA, sizeof(double) * run_globals->params.SnaplistLength, 0, SID.COMM_WORLD);
}

double integrand_time_to_present(double a, void *params)
{
  double omega_m      = ((run_params_t*)params)->OmegaM;
  double omega_k      = ((run_params_t*)params)->OmegaK;
  double omega_lambda = ((run_params_t*)params)->OmegaLambda;

  return 1 / sqrt(omega_m / a + omega_k + omega_lambda * a * a);
}

static double time_to_present(run_globals_t *run_globals, double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time;
  double result;
  double abserr;

  workspace  = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;
  F.params   = &(run_globals->params);

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / run_globals->Hubble,
                      1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / run_globals->Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift
  return time;
}

static void set_units(run_globals_t *run_globals)
{
  run_units_t *units = &(run_globals->units);

  units->UnitTime_in_s         = units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s;
  units->UnitTime_in_Megayears = units->UnitTime_in_s / SEC_PER_MEGAYEAR;

  run_globals->G = GRAVITY / pow(units->UnitLength_in_cm, 3) * units->UnitMass_in_g * pow(units->UnitTime_in_s, 2);

  units->UnitDensity_in_cgs     = units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3);
  units->UnitPressure_in_cgs    = units->UnitMass_in_g / units->UnitLength_in_cm / pow(units->UnitTime_in_s, 2);
  units->UnitCoolingRate_in_cgs = units->UnitPressure_in_cgs / units->UnitTime_in_s;

  units->UnitEnergy_in_cgs = units->UnitMass_in_g * pow(units->UnitLength_in_cm, 2) / pow(units->UnitTime_in_s, 2);

  // convert some physical input parameters to internal units
  run_globals->Hubble = HUBBLE * units->UnitTime_in_s;

  // compute a few quantitites
  run_globals->RhoCrit = 3 * run_globals->Hubble * run_globals->Hubble / (8 * M_PI * run_globals->G);

  // debug("UnitTime_in_s = %e\nUnitTime_in_Megayears = %e\nG = %e\nUnitDensity_in_cgs = %e\nUnitPressure_in_cgs = %e\nUnitCoolingRate_in_cgs = %e\nUnitEnergy_in_cgs = %e\n",
  //     units->UnitTime_in_s, units->UnitTime_in_Megayears, units->UnitDensity_in_cgs, units->UnitPressure_in_cgs, units->UnitCoolingRate_in_cgs, units->UnitEnergy_in_cgs);
  // ABORT(EXIT_SUCCESS);
}

static void read_output_snaps(run_globals_t *run_globals)
{
  int *ListOutputSnaps = run_globals->ListOutputSnaps;
  int *LastOutputSnap  = &(run_globals->LastOutputSnap);

  if (SID.My_rank == 0)
  {
    int i;
    char fname[STRLEN];
    FILE *fd;

    strcpy(fname, run_globals->params.FileWithOutputSnaps);

    if (!(fd = fopen(fname, "r")))
    {
      SID_log_error("file `%s' not found.", fname);
      exit(EXIT_FAILURE);
    }

    for (i = 0; i < NOUT; i++)
    {
      if (fscanf(fd, " %d ", &ListOutputSnaps[i]) != 1)
      {
        SID_log_error("I/O error in file '%s'\n", fname);
        exit(EXIT_FAILURE);
      }
    }
    fclose(fd);

    // Loop through the read in snapshot numbers and convert any negative
    // values to positive ones ala python indexing conventions...
    // e.g. -1 -> MAXSNAPS-1 and so on...
    // Also store the last requested output snapnum
    *LastOutputSnap = 0;
    for (i = 0; i < NOUT; i++)
    {
      if (ListOutputSnaps[i] < 0)
        ListOutputSnaps[i] += run_globals->params.SnaplistLength;
      if (ListOutputSnaps[i] > *LastOutputSnap)
        *LastOutputSnap = ListOutputSnaps[i];
    }

    // sort the list from low to high snapnum
    qsort(ListOutputSnaps, NOUT, sizeof(int), compare_ints);
  }

  // broadcast the data to all other ranks
  SID_Bcast(ListOutputSnaps, sizeof(int) * NOUT, 0, SID.COMM_WORLD);
  SID_Bcast(LastOutputSnap, sizeof(int), 0, SID.COMM_WORLD);
}

void init_meraxes(run_globals_t *run_globals)
{
  int i;
  int snaplist_len;

  // initialise the random number generator
  run_globals->random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(run_globals->random_generator, run_globals->params.RandomSeed);

  // set the units
  set_units(run_globals);

  // init the stdlib random number generator (for CN exceptions only)
  srand((unsigned)time(NULL));

  // read the input snaps list
  read_snap_list(run_globals);

  // read the output snap list
  read_output_snaps(run_globals);
  snaplist_len = run_globals->params.SnaplistLength;

  // read in the requested forest IDs (if any)
  read_requested_forest_ids(run_globals);

  // read in the list of parameters for multiple runs (if required)
  read_multiple_runs_params(run_globals);

  // read in the photometric tables if required
  read_photometric_tables(run_globals);

  // read in the cooling functions
  read_cooling_functions(run_globals);

  for (i = 0; i < snaplist_len; i++)
  {
    run_globals->ZZ[i]     = 1 / run_globals->AA[i] - 1;
    run_globals->LTTime[i] = time_to_present(run_globals, run_globals->ZZ[i]);
  }

  // Initialise galaxy pointers
  run_globals->FirstGal = NULL;
  run_globals->LastGal  = NULL;

  // Calculate the sampled LastSnapshotNr value
  run_globals->params.LastSnapshotNr = (int)(run_globals->params.TotalSimSnaps / run_globals->params.NEverySnap);

  // Prep the output file
  sprintf(run_globals->FNameOut, "%s/%s_%d.hdf5", run_globals->params.OutputDir, run_globals->params.FileNameGalaxies, SID.My_rank);
  prep_hdf5_file(run_globals);

  // Set the SelectForestsSwitch
  run_globals->SelectForestsSwitch = true;

#ifdef USE_TOCF
  set_HII_eff_factor(run_globals);
  malloc_reionization_grids(run_globals);

  if (run_globals->params.TOCF_Flag)
  {
    tocf_params.box_len = run_globals->params.BoxSize / run_globals->params.Hubble_h;

    // Set the cosmology parameters of 21cmFAST to match those of Meraxes
    tocf_params.sigma8      = run_globals->params.Sigma8;
    tocf_params.hlittle     = run_globals->params.Hubble_h;
    tocf_params.OMm         = run_globals->params.OmegaM;
    tocf_params.OMl         = run_globals->params.OmegaLambda;
    tocf_params.OMb         = run_globals->params.BaryonFrac * tocf_params.OMm;
    tocf_params.OMn         = 0.0;
    tocf_params.OMk         = run_globals->params.OmegaK;
    tocf_params.OMr         = run_globals->params.OmegaR;
    tocf_params.power_index = run_globals->params.SpectralIndex;
    tocf_params.wl          = run_globals->params.wLambda;
  }
#else
  // Check if we want to use 21cmFAST
  if (run_globals->params.TOCF_Flag)
    SID_log_warning("TOCF_Flag is set but Meraxes has not been compiled against 21cmFAST!", SID_LOG_COMMENT);
#endif
}

