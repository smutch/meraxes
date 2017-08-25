#define _MAIN
#include <stdio.h>
#include "meraxes.h"
#include "meraxes_gpu.h"
#include "utils.h"

// The code in this routine is all taken from the meraxes codebase.  It
//    is sufficient to set-up all the input parameters and allocate all the
//    grids we need to run find_HII_bubbles()
static void init_meraxes_globals(int argc,char **argv,char *filename_input_params,FILE *fp_log){
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_dup (MPI_COMM_WORLD,&(run_globals.mpi_comm));
    MPI_Comm_rank(MPI_COMM_WORLD,&(run_globals.mpi_rank));
    MPI_Comm_size(MPI_COMM_WORLD,&(run_globals.mpi_size));

    // Direct log files to /dev/null
    init_mlog(MPI_COMM_WORLD,fp_log,fp_log,fp_log);

    // Read the input parameter file
    read_parameter_file(filename_input_params, 0);
    
    // Initiate a few things in the Meraxes global variables structure
    set_units();

    // read the input snaps list
    read_snap_list();

    // read the output snap list
    read_output_snaps();

    // Set redshifts
    int snaplist_len = run_globals.params.SnaplistLength;
    int i = 0;
    for (; i < snaplist_len; i++)
    {
      run_globals.ZZ[i]     = 1 / run_globals.AA[i] - 1;
      run_globals.LTTime[i] = time_to_present(run_globals.ZZ[i]);
    }

    // Allocate reionization grids
    malloc_reionization_grids();
}

static bool inputs_present(char *reference_directory,int snapshot){
    char  fname[STRLEN];
    FILE *fp_in = NULL;
    bool  rval  = true;
    sprintf(fname,"%s/validation_input-core%03d-z%.2f.h5",reference_directory,run_globals.mpi_rank,run_globals.ZZ[snapshot]);
    if((fp_in = fopen(fname,"r")) != NULL)
        fclose(fp_in);
    else
        rval=false;
    return(rval); 
}

int main(int argc,char *argv[]){

    // Parse the command line arguments
    char reference_directory[256];
    char filename_input_params[256];
    bool flag_write_validation=false;
    if(argc!=2 && argc!=3){
       fprintf(stderr,"syntax: %s reference_dir params_file_in_reference_dir [optional: write validation files? (1=yes)]\n",argv[0]);
       exit(EXIT_FAILURE);
    }
    strcpy(reference_directory,argv[1]);
    if(argc==3)
        flag_write_validation=(atoi(argv[2])==1);
    sprintf(filename_input_params,"%s/input.par",reference_directory);
    fprintf(stdout,"Using parameter file: {%s}\n",filename_input_params);

    // Call several meraxes initialisation routines to set-up
    //    the inputs we need, allocate grid arrays, etc.
    FILE *fp_null=fopen("/dev/null","w");
    init_meraxes_globals(argc,argv,filename_input_params,fp_null);

    // Loop over all snapshots
    float time_total=0;
    int   snapshot  =0;
    for(;snapshot<run_globals.params.SnaplistLength;snapshot++){
        // If the input files are available, call find_HII_bubbles
        if(inputs_present(reference_directory,snapshot)){
            // Call find_HII_bubbles
            timer_info timer;
            find_HII_bubbles_driver(snapshot,reference_directory,flag_write_validation,&timer);
            time_total+=timer_delta(timer);
        }
    }
    fprintf(stdout,"Total time spent in _find_HII_bubbles(): %.2f seconds\n",time_total);

    // Clean-up and exit
    free_reionization_grids();
    fclose(fp_null);
    MPI_Comm_free(&(run_globals.mpi_comm));
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

