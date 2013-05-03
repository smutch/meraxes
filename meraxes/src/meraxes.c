#define _MAIN
#include "meraxes.h"
#include <sys/stat.h>
#include <hdf5.h>

static void cleanup(run_globals_struct *run_globals)
{
  H5Tclose(run_globals->hdf5props.array3f_tid);
  SID_free(SID_FARG run_globals->hdf5props.field_types);
  SID_free(SID_FARG run_globals->hdf5props.field_names);
  SID_free(SID_FARG run_globals->hdf5props.dst_field_sizes);
  SID_free(SID_FARG run_globals->hdf5props.dst_offsets);
  gsl_rng_free(run_globals->random_generator);
}

void myexit(int signum)
{
  printf("Task: %d\tnode: %s\tis exiting.\n\n\n", SID.My_rank, SID.My_node);
  SID_exit(signum);
}

static void set_physics_params(
  run_globals_struct *run_globals,
  double             *vals,       
  int                 n_params)   
{

  physics_params_struct *phys_par = &(run_globals->params.physics);

  if ((n_params==3) || (n_params==6)){
    phys_par->peak            = vals[0];
    phys_par->sigma           = vals[1];
    phys_par->stellarfrac     = vals[2];
    if (n_params==6){
      phys_par->peak_evo        = vals[3];
      phys_par->sigma_evo       = vals[4];
      phys_par->stellarfrac_evo = vals[5];
    }
    if (SID.My_rank==0){
      printf("Changed physics_peak to %g\n"        , phys_par->peak);
      printf("Changed physics_sigma to %g\n"       , phys_par->sigma);
      printf("Changed physics_stellarfrac to %g\n" , phys_par->stellarfrac);
      if (n_params==6){
        printf("Changed physics_peak_evo to %g\n"        , phys_par->peak_evo);
        printf("Changed physics_sigma_evo to %g\n"       , phys_par->sigma_evo);
        printf("Changed physics_stellarfrac_evo to %g\n" , phys_par->stellarfrac_evo);
      }
    }
  }

}


int main(int argc, char **argv)
{
  
  struct stat filestatus;

  SID_init(&argc, &argv,NULL);
  
  if(SID.n_proc!=1)
  {
    SID_log_error("Current version of code must be run with ONE CORE (sorry!).");
    ABORT(EXIT_FAILURE);
  }
  
  run_globals_struct run_globals;
  
  if( (argc!=8) && (argc!=4) && (argc!=2) ) 
  {
    if(SID.My_rank==0){
      printf("\n  usage: %s [ -p <paramvals.file> ] <parameterfile> [ <physics.peak> <physics.sigma> <physics.stellarfrac> <physics.peak_evo> <physics.sigma_evo> <physics.stellarfrac_evo> ]\n\n", argv[0]);
      ABORT(1);
    }
  }
  for(int i=1; i<argc-1; i++){
    if (argv[i][0]=='-'){
      switch (argv[i][1]){
        case 'p':
          strcpy(run_globals.params.filename, argv[i+1]);
          break;
        default:
          if(SID.My_rank==0){
            printf("Unrecognised command line argument...\n");
            SID_exit(ERROR_SYNTAX);
          }
          break;
      }
    }
  }

  if(argc==4)
    read_parameter_file(&run_globals, argv[3]);
  else
    read_parameter_file(&run_globals, argv[1]);
  
  // Check to see if the output directory exists and if not, create it
  if (stat(run_globals.params.OutputDir, &filestatus) != 0)
    mkdir(run_globals.params.OutputDir, 02755);
  
  // Deal with any command line parameter values
  if (argc==8){
    double *physics_param_vals = SID_malloc(sizeof(double) * (argc-2));
    for(int i=0; i<argc-2; i++)
      physics_param_vals[i] = atof(argv[i+2]);
    set_physics_params(&run_globals, physics_param_vals, argc-2);
    SID_free(SID_FARG physics_param_vals);
  }

  init_meraxes(&run_globals);
  calc_hdf5_props(&run_globals);

  // Run the model!
  dracarys(&run_globals);

  cleanup(&run_globals);

  SID_exit(EXIT_SUCCESS);
}
