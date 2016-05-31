#include "meraxes.h"
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

static void read_requested_forest_ids()
{
  if (strlen(run_globals.params.ForestIDFile) == 0)
  {
    run_globals.NRequestedForests = -1;
    run_globals.RequestedForestId = NULL;
    return;
  }

  if (SID.My_rank == 0)
  {
    FILE *fin;
    char *line = NULL;
    size_t len;
    int n_forests = -1;
    int *ids;

    if (!(fin = fopen(run_globals.params.ForestIDFile, "r")))
    {
      SID_log_error("Failed to open file: %s", run_globals.params.ForestIDFile);
      ABORT(EXIT_FAILURE);
    }

    getline(&line, &len, fin);
    n_forests                      = atoi(line);
    run_globals.NRequestedForests = n_forests;

    run_globals.RequestedForestId = SID_malloc(sizeof(int) * n_forests);
    ids                            = run_globals.RequestedForestId;

    for (int ii = 0; ii < n_forests; ii++)
    {
      getline(&line, &len, fin);
      ids[ii] = atoi(line);
    }

    free(line);

    SID_log("Found %d requested forest IDs", SID_LOG_COMMENT, n_forests);

    fclose(fin);
  }


  // broadcast the data to all other ranks
  SID_Bcast(&(run_globals.NRequestedForests), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    run_globals.RequestedForestId = SID_malloc(sizeof(int) * run_globals.NRequestedForests);
  SID_Bcast(run_globals.RequestedForestId, sizeof(int) * run_globals.NRequestedForests, 0, SID.COMM_WORLD);
}


static void read_snap_list()
{
  if (SID.My_rank == 0)
  {
    FILE *fin;
    int snaplist_len;
    double dummy;
    char fname[STRLEN];
    run_params_t params = run_globals.params;

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

    run_globals.params.SnaplistLength = snaplist_len;
    if (SID.My_rank == 0)
      printf("found %d defined times in snaplist.\n", snaplist_len);

    // malloc the relevant arrays
    run_globals.AA     = SID_malloc(sizeof(double) * snaplist_len);
    run_globals.ZZ     = SID_malloc(sizeof(double) * snaplist_len);
    run_globals.LTTime = SID_malloc(sizeof(double) * snaplist_len);

    // seek back to the start of the file
    rewind(fin);

    // actually read in the expansion factors
    snaplist_len = 0;
    do
    {
      if (fscanf(fin, " %lg ", &(run_globals.AA[snaplist_len])) == 1)
        snaplist_len++;
      else
        break;
    } while (true);

    // close the file
    fclose(fin);
  }

  // broadcast the read to all other ranks and malloc the necessary arrays
  SID_Bcast(&(run_globals.params.SnaplistLength), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
  {
    run_globals.AA     = SID_malloc(sizeof(double) * run_globals.params.SnaplistLength);
    run_globals.ZZ     = SID_malloc(sizeof(double) * run_globals.params.SnaplistLength);
    run_globals.LTTime = SID_malloc(sizeof(double) * run_globals.params.SnaplistLength);
  }
  SID_Bcast(run_globals.AA, sizeof(double) * run_globals.params.SnaplistLength, 0, SID.COMM_WORLD);
}

double integrand_time_to_present(double a, void *params)
{
  double omega_m      = ((run_params_t*)params)->OmegaM;
  double omega_k      = ((run_params_t*)params)->OmegaK;
  double omega_lambda = ((run_params_t*)params)->OmegaLambda;

  return 1 / sqrt(omega_m / a + omega_k + omega_lambda * a * a);
}

static double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time;
  double result;
  double abserr;

  workspace  = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;
  F.params   = &(run_globals.params);

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / run_globals.Hubble,
                      1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / run_globals.Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift
  return time;
}

void set_units()
{
  run_units_t *units = &(run_globals.units);

  units->UnitTime_in_s         = units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s;
  units->UnitTime_in_Megayears = units->UnitTime_in_s / SEC_PER_MEGAYEAR;

  run_globals.G = GRAVITY / pow(units->UnitLength_in_cm, 3) * units->UnitMass_in_g * pow(units->UnitTime_in_s, 2);

  units->UnitDensity_in_cgs     = units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3);
  units->UnitPressure_in_cgs    = units->UnitMass_in_g / units->UnitLength_in_cm / pow(units->UnitTime_in_s, 2);
  units->UnitCoolingRate_in_cgs = units->UnitPressure_in_cgs / units->UnitTime_in_s;

  units->UnitEnergy_in_cgs = units->UnitMass_in_g * pow(units->UnitLength_in_cm, 2) / pow(units->UnitTime_in_s, 2);

  // convert some physical input parameters to internal units
  run_globals.Hubble = HUBBLE * units->UnitTime_in_s;

  // compute a few quantitites
  run_globals.RhoCrit = 3 * run_globals.Hubble * run_globals.Hubble / (8 * M_PI * run_globals.G);

  // debug("UnitTime_in_s = %e\nUnitTime_in_Megayears = %e\nG = %e\nUnitDensity_in_cgs = %e\nUnitPressure_in_cgs = %e\nUnitCoolingRate_in_cgs = %e\nUnitEnergy_in_cgs = %e\n",
  //     units->UnitTime_in_s, units->UnitTime_in_Megayears, units->UnitDensity_in_cgs, units->UnitPressure_in_cgs, units->UnitCoolingRate_in_cgs, units->UnitEnergy_in_cgs);
  // ABORT(EXIT_SUCCESS);
}

static void read_output_snaps()
{
  int **ListOutputSnaps = &(run_globals.ListOutputSnaps);
  int *LastOutputSnap  = &(run_globals.LastOutputSnap);
  int maxsnaps = run_globals.params.SnaplistLength; 
  int *nout = &(run_globals.NOutputSnaps);

  if (SID.My_rank == 0)
  {
    int i;
    char fname[STRLEN];
    FILE *fd;

    strcpy(fname, run_globals.params.FileWithOutputSnaps);

    if (!(fd = fopen(fname, "r")))
    {
      SID_log_error("file `%s' not found.", fname);
      exit(EXIT_FAILURE);
    }

    // find out how many output snapshots are being requested
    int dummy;
    for (i = 0; i < maxsnaps; i++)
    {
      if (fscanf(fd, " %d ", &dummy) == 1)
        (*nout)++;
      else
        break;
    }
    fseek(fd, 0, SEEK_SET);

    // allocate the ListOutputSnaps array
    *ListOutputSnaps = SID_malloc(sizeof(int) * (*nout));

    for (i = 0; i < (*nout); i++)
    {
      if (fscanf(fd, " %d ", &((*ListOutputSnaps)[i])) != 1)
      {
        SID_log_error("I/O error in file '%s'\n", fname);
        exit(EXIT_FAILURE);
      }
    }
    fclose(fd);

#ifdef CALC_MAGS
    if (*nout != NOUT)
    {
      SID_log_error("Number of entries in output snaplist does not match NOUT!");
      ABORT(EXIT_FAILURE);
    }
#endif

    // Loop through the read in snapshot numbers and convert any negative
    // values to positive ones ala python indexing conventions...
    // e.g. -1 -> MAXSNAPS-1 and so on...
    // Also store the last requested output snapnum
    *LastOutputSnap = 0;
    for (i = 0; i < (*nout); i++)
    {
      if ((*ListOutputSnaps)[i] < 0)
        (*ListOutputSnaps)[i] += run_globals.params.SnaplistLength;
      if ((*ListOutputSnaps)[i] > *LastOutputSnap)
        *LastOutputSnap = (*ListOutputSnaps)[i];
    }

    // DEBUG
    SID_log("nout = %d", SID_LOG_COMMENT, *nout);
    SID_log("LastOutputSnap = %d", SID_LOG_COMMENT, *LastOutputSnap);
    SID_log("ListOutputSnaps = [ ", SID_LOG_CONTINUE);
    for(int ii=0; ii < *nout; ii++)
      SID_log("%d ", SID_LOG_CONTINUE, (*ListOutputSnaps)[ii]);
    SID_log("]", SID_LOG_CONTINUE);

    // sort the list from low to high snapnum
    qsort(*ListOutputSnaps, (*nout), sizeof(int), compare_ints);

  }

  // broadcast the data to all other ranks
  SID_Bcast(nout, sizeof(int), 0, SID.COMM_WORLD);
  
  if(SID.My_rank > 0)
    *ListOutputSnaps = SID_malloc(sizeof(int) * (*nout));

  SID_Bcast(*ListOutputSnaps, sizeof(int) * (*nout), 0, SID.COMM_WORLD);
  SID_Bcast(LastOutputSnap, sizeof(int), 0, SID.COMM_WORLD);

}


static void check_n_history_snaps()
{
  // Check that N_HISTORY_SNAPS is set to a high enough value to allow all
  // SN-II to be tracked across the entire simulation.  This is calculated in
  // an extremely crude fasion!

  double *LTTime = run_globals.LTTime;
  int n_snaps    = run_globals.params.SnaplistLength;
  double min_dt  = LTTime[0] - LTTime[n_snaps-1];
  double m_low;

  for (int ii = 0; ii < n_snaps-N_HISTORY_SNAPS; ii++)
  {
    double diff = LTTime[ii] - LTTime[ii+N_HISTORY_SNAPS];
    if (diff < min_dt)
      min_dt = diff;
  }

  min_dt *= run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h;
  m_low   = sn_m_low(log10(min_dt));

  if (m_low > 8.0)
  {
    SID_log_error("N_HISTORY_SNAPS is likely not set to a high enough value!  Exiting...");
    ABORT(EXIT_FAILURE);
  }
}


void init_meraxes()
{
  int i;
  int snaplist_len;

  // initialise the random number generator
  run_globals.random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(run_globals.random_generator, run_globals.params.RandomSeed);

  // set the units
  set_units();

  // init the stdlib random number generator (for CN exceptions only)
  srand((unsigned)time(NULL));

  // read the input snaps list
  read_snap_list();

  // check to ensure N_HISTORY_SNAPS is set to a high enough value
  check_n_history_snaps();

  // read the output snap list
  read_output_snaps();
  snaplist_len = run_globals.params.SnaplistLength;

  // read in the requested forest IDs (if any)
  read_requested_forest_ids();

  // read in the photometric tables if required
  read_photometric_tables();

  // read in the cooling functions
  read_cooling_functions();

  // read in the mean Mvir_crit table (if needed)
  read_Mcrit_table();

  for (i = 0; i < snaplist_len; i++)
  {
    run_globals.ZZ[i]     = 1 / run_globals.AA[i] - 1;
    run_globals.LTTime[i] = time_to_present(run_globals.ZZ[i]);
  }

  // Initialise galaxy pointers
  run_globals.FirstGal = NULL;
  run_globals.LastGal  = NULL;

  // Initialise some book keeping parameters for the input trees
  run_globals.TreesStep = -1;
  run_globals.TreesScan = -1;

  // Set the SelectForestsSwitch
  run_globals.SelectForestsSwitch = true;

  // Initialize the halo storage arrays
  initialize_halo_storage();

  malloc_reionization_grids();
  set_ReionEfficiency();
}

