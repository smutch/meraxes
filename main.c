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

// Call the appropriate version of find_HII_bubbles(), depending
//    on what compiler flags have been set.
static void call_find_HII_bubbles_gpu(int snapshot,char *reference_directory,timer_info *timer){
#ifdef __NVCC__
#ifndef USE_CUFFT
    fprintf(stdout,"Calling hybrid-GPU/FFTW version of find_HII_bubbles()...");fflush(stdout);
#else
    fprintf(stdout,"Calling pure-GPU version of find_HII_bubbles()...");fflush(stdout);
#endif
    // Run the GPU version of _find_HII_bubbles()
    find_HII_bubbles_driver(snapshot,_find_HII_bubbles_gpu,reference_directory,timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(*timer));fflush(stdout);
#else
    // Run the Meraxes version of _find_HII_bubbles()
    fprintf(stdout,"Calling Meraxes version of find_HII_bubbles()...");fflush(stdout);
    find_HII_bubbles_driver(snapshot,_find_HII_bubbles,reference_directory,timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(*timer));fflush(stdout);
#endif
}

int main(int argc,char *argv[]){

    // Parse the command line arguments
    char reference_directory[256];
    char filename_input_params[256];
    if(argc!=2){
       fprintf(stderr,"You need to pass two arguments (the reference directory and the parameter file residing in it).\n");
       exit(EXIT_FAILURE);
    }
    strcpy(reference_directory,argv[1]);
    sprintf(filename_input_params,"%s/input.par",reference_directory);
    fprintf(stdout,"Using parameter file: {%s}\n",filename_input_params);

    // Call several meraxes initialisation routines to set-up
    //    the inputs we need, allocate grid arrays, etc.
    FILE *fp_null=fopen("/dev/null","w");
    init_meraxes_globals(argc,argv,filename_input_params,fp_null);

    // Create a timer 
    timer_info timer;

    // For the given redshift, which given output are we working with?
    double redshift = 6.99; // xH=0.5
    int    snapshot = 0;
    int    j_out    = 0;
    double max_diff = 1e5;
    for(;j_out<run_globals.params.SnaplistLength;j_out++){
        double diff=fabs(run_globals.ZZ[j_out]-redshift);
        if(diff<max_diff){
            snapshot=j_out;
            max_diff=diff;
        }
    }
    redshift=run_globals.ZZ[snapshot];
    fprintf(stdout,"Processing snapshot %d of %d (z=%lf).\n",snapshot,run_globals.params.SnaplistLength,redshift);fflush(stdout);

    // Run find_HII_bubbles
    call_find_HII_bubbles_gpu(snapshot,reference_directory,&timer);

    // Clean-up and exit
    free_reionization_grids();
    fclose(fp_null);
    MPI_Comm_free(&(run_globals.mpi_comm));
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

