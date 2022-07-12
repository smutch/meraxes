#include <assert.h>
#include <gsl/gsl_integration.h>
#include <time.h>

#include "ComputePowerSpectrum.h"
#include "ConstructLightcone.h"
#include "cooling.h"
#include "init.h"
#include "magnitudes.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "parse_paramfile.h"
#include "read_halos.h"
#include "recombinations.h"
#include "reionization.h"
#include "reionization_modifiers.h"
#include "save.h"
#include "stellar_feedback.h"

static void init_gpu()
{
// If we are compiling with CUDA, allocate a structure
//   that will carry information about the device
#ifdef USE_CUDA
  // Alocate the structure that will carry all the information
  //   about the GPU assigned to this thread
  run_globals.gpu = (gpu_info*)malloc(sizeof(gpu_info));

  // This function has all the CUDA device polling calls
  init_CUDA();
#else
  // If we are not compiling with CUDA, set this
  //   pointer to NULL.  This is a good way
  //   to test in the code if a GPU is being used.
  mlog("CPU-only version of Meraxes running.", MLOG_MESG);
  run_globals.gpu = NULL;
#endif
}

static void read_requested_forest_ids()
{
  if (strlen(run_globals.params.ForestIDFile) == 0) {
    run_globals.NRequestedForests = -1;
    run_globals.RequestedForestId = NULL;
    return;
  }

  if (run_globals.mpi_rank == 0) {

    FILE* fin;
    if (!(fin = fopen(run_globals.params.ForestIDFile, "r"))) {
      mlog_error("Failed to open file: %s", run_globals.params.ForestIDFile);
      ABORT(EXIT_FAILURE);
    }

    char* line = NULL;
    size_t len;
    getline(&line, &len, fin);
    int n_forests = atoi(line);
    run_globals.NRequestedForests = n_forests;

    run_globals.RequestedForestId = malloc(sizeof(long) * n_forests);
    long* ids = run_globals.RequestedForestId;

    for (int ii = 0; ii < n_forests; ii++) {
      getline(&line, &len, fin);
      ids[ii] = atol(line);
    }

    free(line);

    mlog("Found %d requested forest IDs", MLOG_MESG, n_forests);

    fclose(fin);
  }

  // broadcast the data to all other ranks
  MPI_Bcast(&(run_globals.NRequestedForests), 1, MPI_INT, 0, run_globals.mpi_comm);
  if (run_globals.mpi_rank > 0)
    run_globals.RequestedForestId = malloc(sizeof(int) * run_globals.NRequestedForests);
  MPI_Bcast(run_globals.RequestedForestId, run_globals.NRequestedForests, MPI_INT, 0, run_globals.mpi_comm);
}

static void read_snap_list()
{
  if (run_globals.mpi_rank == 0) {
    FILE* fin;
    int snaplist_len;
    double dummy;
    char fname[STRLEN + 12];
    run_params_t params = run_globals.params;

    sprintf(fname, "%s/a_list.txt", params.SimulationDir);

    if (!(fin = fopen(fname, "r"))) {
      mlog_error("failed to read snaplist in file '%s'", fname);
      ABORT(EXIT_FAILURE);
    }

    // Count the number of snapshot list entries
    snaplist_len = 0;
    do {
      if (fscanf(fin, " %lg ", &dummy) == 1)
        snaplist_len++;
      else
        break;
    } while (true);

    run_globals.params.SnaplistLength = snaplist_len;
    if (run_globals.mpi_rank == 0)
      mlog("found %d defined times in snaplist.\n", MLOG_MESG, snaplist_len);

    // malloc the relevant arrays
    run_globals.AA = malloc(sizeof(double) * snaplist_len);
    run_globals.ZZ = malloc(sizeof(double) * snaplist_len);
    run_globals.LTTime = malloc(sizeof(double) * snaplist_len);

    // seek back to the start of the file
    rewind(fin);

    // actually read in the expansion factors
    snaplist_len = 0;
    do {
      if (fscanf(fin, " %lg ", &(run_globals.AA[snaplist_len])) == 1)
        snaplist_len++;
      else
        break;
    } while (true);

    // close the file
    fclose(fin);
  }

  // broadcast the read to all other ranks and malloc the necessary arrays
  MPI_Bcast(&(run_globals.params.SnaplistLength), 1, MPI_INT, 0, run_globals.mpi_comm);
  if (run_globals.mpi_rank > 0) {
    run_globals.AA = malloc(sizeof(double) * run_globals.params.SnaplistLength);
    run_globals.ZZ = malloc(sizeof(double) * run_globals.params.SnaplistLength);
    run_globals.LTTime = malloc(sizeof(double) * run_globals.params.SnaplistLength);
  }
  MPI_Bcast(run_globals.AA, run_globals.params.SnaplistLength, MPI_DOUBLE, 0, run_globals.mpi_comm);
}

double integrand_time_to_present(double a, void* params)
{
  double omega_m = ((run_params_t*)params)->OmegaM;
  double omega_k = ((run_params_t*)params)->OmegaK;
  double omega_lambda = ((run_params_t*)params)->OmegaLambda;

  return 1 / sqrt(omega_m / a + omega_k + omega_lambda * a * a);
}

static double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace* workspace;
  double time;
  double result;
  double abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;
  F.params = &(run_globals.params);

  gsl_integration_qag(
    &F, 1.0 / (z + 1), 1.0, 1.0 / run_globals.Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / run_globals.Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift
  return time;
}

void set_units()
{
  run_units_t* units = &(run_globals.units);

  units->UnitTime_in_s = units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s;
  units->UnitTime_in_Megayears = units->UnitTime_in_s / SEC_PER_MEGAYEAR;

  run_globals.G = GRAVITY / pow(units->UnitLength_in_cm, 3) * units->UnitMass_in_g * pow(units->UnitTime_in_s, 2);
  run_globals.Csquare = pow(SPEED_OF_LIGHT / units->UnitVelocity_in_cm_per_s, 2);

  units->UnitDensity_in_cgs = units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3);
  units->UnitPressure_in_cgs = units->UnitMass_in_g / units->UnitLength_in_cm / pow(units->UnitTime_in_s, 2);
  units->UnitCoolingRate_in_cgs = units->UnitPressure_in_cgs / units->UnitTime_in_s;

  units->UnitEnergy_in_cgs = units->UnitMass_in_g * pow(units->UnitLength_in_cm, 2) / pow(units->UnitTime_in_s, 2);

  // convert some physical input parameters to internal units
  run_globals.Hubble = HUBBLE * units->UnitTime_in_s;

  // compute a few quantitites
  run_globals.RhoCrit = 3 * run_globals.Hubble * run_globals.Hubble / (8 * M_PI * run_globals.G);
}

void init_storage()
{
  // Initialize the halo storage arrays
  initialize_halo_storage();

  malloc_reionization_grids();

  // calculate the output hdf5 file properties for later use
  calc_hdf5_props();
}

void init_meraxes()
{
  int i;
  int snaplist_len;

  // initialize GPU
  init_gpu();

  // initialise the random number generator
  run_globals.random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(run_globals.random_generator, (unsigned long)run_globals.params.RandomSeed);

  // set the units
  set_units();

  // init the stdlib random number generator (for CN exceptions only)
  srand((unsigned)time(NULL));

  // read the input snaps list
  read_snap_list();

  // parse the requested output snaps
  parse_output_snaps(run_globals.params.OutputSnapsString);

  snaplist_len = run_globals.params.SnaplistLength;
  for (i = 0; i < snaplist_len; i++) {
    run_globals.ZZ[i] = 1 / run_globals.AA[i] - 1;
    run_globals.LTTime[i] = time_to_present(run_globals.ZZ[i]);
  }

  // validation checks
  if (run_globals.params.Flag_IncludeSpinTemp) {
    if (run_globals.params.physics.ReionMaxHeatingRedshift > run_globals.ZZ[0]) {
      mlog_error("ReionMaxHeatingRedshift > redshift of first snapshot (%.2f) which is not allowed.",
                 run_globals.ZZ[0]);
      ABORT(EXIT_FAILURE);
    }
  }

  // read in the requested forest IDs (if any)
  read_requested_forest_ids();

  // read in the cooling functions
  read_cooling_functions();

  // read in the stellar feedback tables
  read_stellar_feedback_tables();

#ifdef CALC_MAGS
  init_magnitudes();
#endif

  // set RequestedMassRatioModifier and RequestedBaryonFracModifieruto be 1 first
  // it will be set to -1 later if MassRatioModifier or BaryonFracModifier is not specified
  run_globals.RequestedMassRatioModifier = 1;
  run_globals.RequestedBaryonFracModifier = 1;

  // read in the mean Mvir_crit table (if needed 1 for Reio 2 for LW)
  read_Mcrit_table(1);
  read_Mcrit_table(2);

  set_ReionEfficiency();
  set_quasar_fobs();

  // Determine the size of the light-cone for initialising the light-cone grid
  if (run_globals.params.Flag_PatchyReion && run_globals.params.Flag_ConstructLightcone) {
    Initialise_ConstructLightcone();
  }

  if (run_globals.params.Flag_ComputePS) {
    Initialise_PowerSpectrum();
  }

  // Initialise interpolation tables for inhomogeneous recombinations
  // TODO: Should this also depend on Flag_PatchyReion? (apply decision to cleanup.c too)
  if (run_globals.params.Flag_IncludeRecombinations) {
    init_MHR();
  }

  // Initialise galaxy pointers
  run_globals.FirstGal = NULL;
  run_globals.LastGal = NULL;

  // Set the SelectForestsSwitch
  run_globals.SelectForestsSwitch = true;

  // This will be set by Mhysa
  run_globals.mhysa_self = NULL;
}
