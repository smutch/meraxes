#define _MAIN
#include <stdio.h>
#include "meraxes.h"
#include "meraxes_gpu.h"
#include "utils.h"

int main(int argc,char *argv[]){

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_dup (MPI_COMM_WORLD,&(run_globals.mpi_comm));
    MPI_Comm_rank(MPI_COMM_WORLD,&(run_globals.mpi_rank));
    MPI_Comm_size(MPI_COMM_WORLD,&(run_globals.mpi_size));

    // Direct log files to /dev/null
    FILE *fp_null=fopen("/dev/null","w");
    init_mlog(MPI_COMM_WORLD,fp_null,fp_null,fp_null);

    // Parse the command line arguments
    char reference_directory[256];
    char filename_input_params[256];
    if(argc!=2){
       fprintf(stderr,"You need to pass one argument (the reference directory).\n");
       exit(EXIT_FAILURE);
    }
    strcpy(reference_directory,argv[1]);
    sprintf(filename_input_params,"%s/input.par",reference_directory);
    fprintf(stdout,"Using parameter file: {%s}\n",filename_input_params);

    // Read the input parameter file
    read_parameter_file(filename_input_params, 0);
    
    // Initiate a few things in the Meraxes global variables structure
    set_units();

    // Initialize a timer and the redshift we will process
    timer_info   timer;
    //const double redshift = 6.99f;
    //const double redshift = 11.30f;
    const double redshift =7.145f;

#ifdef __NVCC__
    // Run the GPU version of _find_HII_bubbles()
    fprintf(stdout,"Calling GPU version of find_HII_bubbles()...");fflush(stdout);
    find_HII_bubbles_driver(redshift,_find_HII_bubbles_gpu,reference_directory,&timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(timer));fflush(stdout);
#else
    // Run the Meraxes version of _find_HII_bubbles()
    fprintf(stdout,"Calling Meraxes version of find_HII_bubbles()...");fflush(stdout);
    find_HII_bubbles_driver(redshift,_find_HII_bubbles,reference_directory,&timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(timer));fflush(stdout);
#endif

    // Clean-up and exit
    fclose(fp_null);
    MPI_Comm_free(&(run_globals.mpi_comm));
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

