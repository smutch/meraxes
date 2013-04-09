#define _MAIN
#include "meraxes.h"

void myexit(int signum)
{
  printf("Task: %d\tnode: %s\tis exiting.\n\n\n", ThisTask, ThisNode);
  exit(signum);
}

static void set_physics_params(double *vals, int n_params)
{

  if ((n_params==3) || (n_params==6)){
    run_params.physics.peak            = vals[0];
    run_params.physics.sigma           = vals[1];
    run_params.physics.stellarfrac     = vals[2];
    if (n_params==6){
      run_params.physics.peak_evo        = vals[3];
      run_params.physics.sigma_evo       = vals[4];
      run_params.physics.stellarfrac_evo = vals[5];
    }
    if (ThisTask==0){
      printf("Changed physics_peak to %g\n"        , run_params.physics.peak);
      printf("Changed physics_sigma to %g\n"       , run_params.physics.sigma);
      printf("Changed physics_stellarfrac to %g\n" , run_params.physics.stellarfrac);
      if (n_params==6){
        printf("Changed physics_peak_evo to %g\n"        , run_params.physics.peak_evo);
        printf("Changed physics_sigma_evo to %g\n"       , run_params.physics.sigma_evo);
        printf("Changed physics_stellarfrac_evo to %g\n" , run_params.physics.stellarfrac_evo);
      }
    }
  }

}


int main(int argc, char **argv)
{
  
  SID_init(argc, &argv,NULL);
  
  int opt_paramsval_list = 0;
  if( (argc!=8) && (argc!=4) && (argc!=2) ) 
  {
    if(ThisTask==0){
      printf("\n  usage: %s [ -p <paramvals.file> ] <parameterfile> [ <physics.peak> <physics.sigma> <physics.stellarfrac> <physics.peak_evo> <physics.sigma_evo> <physics.stellarfrac_evo> ]\n\n", argv[0]);
      ABORT(1);
    }
  }
  for(int i=1; i<argc-1; i++){
    if (argv[i][0]=='-'){
      switch (argv[i][1]){
        case 'p':
          opt_paramsval_list = 1;
          strcpy(paramvals_file, argv[i+1]);
          break;
        default:
          if(ThisTask==0){
            printf("Unrecognised command line argument...\n");
            SID_exit(ERROR_SYNTAX);
          }
          break;
      }
    }
  }

  if(argc==4)
    read_parameter_file(argv[3]);
  else
    read_parameter_file(argv[1]);
  
  // Check to see if the output directory exists and if not, create it
  if (stat(run_params.output_dir, &filestatus) != 0)
    mkdir(run_params.output_dir, 02755);
  
  // Deal with any command line parameter values
  if (argc==8){
    double *physics_param_vals = malloc(sizeof(double) * (argc-2));
    for(int i=0; i<argc-2; i++)
      physics_param_vals[i] = atof(argv[i+2]);
    set_physics_params(physics_param_vals, argc-2);
    free(physics_param_vals);
  }

  return 0;
}
