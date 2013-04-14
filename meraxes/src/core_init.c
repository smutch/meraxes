#include "meraxes.h"
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

static int read_snap_list(run_globals_struct *run_globals)
{
  FILE *fin;
  int snaplist_len;
  char fname[STRLEN];
  run_params_struct params = run_globals->params;

  sprintf(fname, "%s/%s/trees/%s.a_list",
      params.SimulationDir,
      params.SimName,
      params.SimName);
  
  if(!(fin = fopen(fname, "r")))
  {
    SID_log_error("failed to read snaplist in file '%s'", fname);
    ABORT(EXIT_FAILURE);
  }

  snaplist_len = 0;
  do
  {
    if(fscanf(fin, " %lg ", &(run_globals->AA[snaplist_len])) == 1)
      snaplist_len++;
    else
      break;
  } while(snaplist_len < MAXSNAPS);

  fclose(fin);

  if(SID.My_rank == 0){
    printf("NOUT = %d\n\n", NOUT);
    printf("found %d defined times in snaplist.\n", snaplist_len);
  }

  return snaplist_len;

}

double integrand_time_to_present(double a, void *params)
{
  double omega = ((run_params_struct *)params)->Omega;
  double omega_lambda = ((run_params_struct *)params)->OmegaLambda;
  return 1 / sqrt(omega / a + (1 - omega - omega_lambda) + omega_lambda * a * a);
}

static double time_to_present(run_globals_struct *run_globals, double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;
  F.params = &(run_globals->params);

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / run_globals->Hubble,
    1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / run_globals->Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift 
  return time;
}

static void set_units(run_globals_struct *run_globals)
{
  run_units_struct *units       = &(run_globals->units);

  units->UnitTime_in_s          = units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s;
  units->UnitTime_in_Megayears  = units->UnitTime_in_s / SEC_PER_MEGAYEAR;

  run_globals->G                = GRAVITY / pow(units->UnitLength_in_cm, 3) * units->UnitMass_in_g * pow(units->UnitTime_in_s, 2);

  units->UnitDensity_in_cgs     = units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3);
  units->UnitPressure_in_cgs    = units->UnitMass_in_g / units->UnitLength_in_cm / pow(units->UnitTime_in_s, 2);
  units->UnitCoolingRate_in_cgs = units->UnitPressure_in_cgs / units->UnitTime_in_s;

  units->UnitEnergy_in_cgs      = units->UnitMass_in_g * pow(units->UnitLength_in_cm, 2) / pow(units->UnitTime_in_s, 2);

  // convert some physical input parameters to internal units 
  run_globals->Hubble           = HUBBLE * units->UnitTime_in_s;

  // compute a few quantitites 
  run_globals->RhoCrit          = 3 * run_globals->Hubble * run_globals->Hubble / (8 * M_PI * run_globals->G);
}

static void read_output_snaps(run_globals_struct *run_globals)
{
  int i;

  char fname[STRLEN];
  int  *ListOutputSnaps = run_globals->ListOutputSnaps;
  int  *LastOutputSnap  = &(run_globals->LastOutputSnap);
  FILE *fd;

  strcpy(fname, run_globals->params.FileWithOutputSnaps);

  if(!(fd = fopen(fname, "r")))
  {
    SID_log_error("file `%s' not found.", fname);
    exit(EXIT_FAILURE);
  }

  for(i = 0; i < NOUT; i++)
  {
    if(fscanf(fd, " %d ", &ListOutputSnaps[i]) != 1)
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
  for (i = 0; i < NOUT; i++) {
    if(ListOutputSnaps[i]<0)
      ListOutputSnaps[i] += MAXSNAPS;
    if(ListOutputSnaps[i]>*LastOutputSnap)
      *LastOutputSnap=ListOutputSnaps[i];
  }

}

void init_meraxis(run_globals_struct *run_globals)
{
  int i;
  int snaplist_len;

  run_globals->random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(run_globals->random_generator, 42);	 // start-up seed 

  set_units(run_globals);
  srand((unsigned) time(NULL));

  read_output_snaps(run_globals);
  snaplist_len = read_snap_list(run_globals);

  for(i = 0; i < snaplist_len; i++)
  {
    run_globals->ZZ[i] = 1 / run_globals->AA[i] - 1;
    run_globals->Age[i] = time_to_present(run_globals, run_globals->ZZ[i]);
  }

}

