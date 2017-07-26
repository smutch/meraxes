#define _MAIN
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>
#include "meraxes.h"
#include "meraxes_gpu.h"
#include "utils.h"

int main(int argc,char *argv[]){

    // Initialize MPI and direct log files to 'dev'null
    FILE *fp_null=fopen("/dev/null","w");
    MPI_Init(&argc, &argv);
    init_mlog(MPI_COMM_WORLD,fp_null,fp_null,fp_null);

    // Process z=6.99 snapshot
    timer_info timer;
    double     redshift = 6.99f;
    fprintf(stdout,"Calling Meraxes version of find_HII_bubbles()...");fflush(stdout);
    find_HII_bubbles_driver(redshift,_find_HII_bubbles,&timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(timer));fflush(stdout);
    fprintf(stdout,"Calling GPU version of find_HII_bubbles()...");fflush(stdout);
    find_HII_bubbles_driver(redshift,_find_HII_bubbles_gpu,&timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(timer));fflush(stdout);

    // Clean-up
    MPI_Finalize();
    fclose(fp_null);
    exit(0);
}

