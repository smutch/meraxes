#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/wait.h>

#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "nbody.h"

int main(int argc, char **argv)
{
  float Z, M, M_MIN, nf;
  char cmnd[1000];
  FILE *LOG;
  time_t start_time, curr_time;
  int status;
  float zlist[MAXSNAPS];

  time(&start_time);

  // make appropriate directories
  system("mkdir ../Log_files");
  system("mkdir ../Boxes");
  system("mkdir ../Output_files");
  system("mkdir ../Output_files/DNDLNM_files");
  system("mkdir ../Output_files/FgtrM_files");
  system("mkdir ../Output_files/Halo_lists");
  system("mkdir ../Output_files/Size_distributions");
  system("mkdir ../Output_files/Deldel_T_power_spec");
  system("mkdir ../Redshift_interpolate_filelists");
  system("mkdir ../Lighttravel_filelists");
  init_ps();

  // open log file
  system("rm ../Log_files/*");
  LOG = fopen("../Log_files/drive_logzscroll_log_file", "w");
  if (!LOG){
    fprintf(stderr, "drive_zscroll_log_file.c: Unable to open log file\n Aborting...\n");
    return -1;
  }

  read_zlist(zlist);

  for (int snapnum=FIRSTSNAP; snapnum<LASTSNAP; snapnum++){

    fprintf(stderr, "*************************************\n");

    //set the minimum source mass
    Z = zlist[snapnum];
    M_MIN = get_M_min_ion(Z);

    if (!USE_HALO_FIELD)
    {
      fprintf("USE_HALO_FIELD not set!\n");
      ABORT(EXIT_FAILURE);
    }

    // if it is the lowest redshift, let's call Ts.c
    if ((USE_TS_IN_21CM && (fabs(Z-ZLOW)/Z < 0.0002) ) || (Z > Z_HEAT_MAX)){ // NEW CONDITIONAL
      //    if (USE_TS_IN_21CM && (fabs(Z-ZLOW)/Z < 0.0002) ){
      sprintf(cmnd, "Ts %d", snapnum);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    } // this will create all of the higher z Ts files in Boxes, provided Ts_verbose is turned on
    // in HEAT_PARAMS.H

    // find bubbles
    sprintf(cmnd, "find_HII_bubbles %d", snapnum);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    status = system(cmnd);
    nf = WEXITSTATUS(status) / 100.0;
    if (nf < 0){
      fprintf(stderr, "find_HII_bubbles exited...\nAborting run...\n");
      fprintf(LOG,  "find_HII_bubbles exited...\nAborting run...\n");
      return -1;
    }

    // do temperature map
    switch(FIND_BUBBLE_ALGORITHM)
    {
      case 2:
        sprintf(cmnd, "delta_T %06.2f ../Boxes/xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc ../Boxes/Ts_z%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_*Pop%i_%i*_%.0fMpc", Z, Z, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN, Z, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, Pop, HII_DIM, BOX_LEN);
        break;
      default:
        sprintf(cmnd, "delta_T %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc ../Boxes/Ts_z%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_*Pop%i_%i*_%.0fMpc", Z, Z, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN, Z, ZETA_X, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, Pop, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    fprintf(stderr, "*************************************\n");
    fflush(NULL);

    }

    fclose(LOG);
    return 0;
  }
