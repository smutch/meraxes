#define _MAIN
#include <stdio.h>
#include "meraxes.h"
#include "meraxes_gpu.h"
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>

int main(int argc,char *argv[]){

    // Initialize MPI and direct log files to 'dev'null
    FILE *fp_null=fopen("/dev/null","w");
    MPI_Init(&argc, &argv);
    init_mlog(MPI_COMM_WORLD,fp_null,fp_null,fp_null);

    // Process z=6.99 snapshot
    double redshift = 6.99f;
    find_HII_bubbles_driver(redshift,_find_HII_bubbles_gpu);

    // Clean-up
    MPI_Finalize();
    fclose(fp_null);
    exit(0);
}

