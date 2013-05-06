#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "meraxes.h"

void mpi_debug_here()
{
  int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  printf("Once connected go up stack to 'sleep(5)' and 'set var i=7'\n");
  fflush(stdout);
  while (0 == i)
    sleep(5);
}

void check_counts(run_globals_struct *run_globals, fof_group_struct *fof_group, int NGal, int NFof)
{

  int counter = 0;
  int halo_counter = 0;
  galaxy_struct *gal = NULL;
  halo_struct *halo = NULL;

  SID_log("Running counts check...", SID_LOG_OPEN|SID_LOG_TIMER);

  SID_log("NFof = %d", SID_LOG_COMMENT, NFof);
  SID_log("Ngal = %d", SID_LOG_COMMENT, NGal);

  counter=0;
  gal = run_globals->FirstGal;
  while (gal!=NULL)
  {
    counter++;
    gal = gal->Next;
  }
  SID_log("Counting using gal->Next gives %d gals", SID_LOG_COMMENT, counter);

  counter = 0;
  halo_counter = 0;
  int ii, jj;
  for(int i_fof=0; i_fof<NFof; i_fof++)
  {
    halo = fof_group[i_fof].FirstHalo;
    jj=0;
    while (halo!=NULL) {
      gal = halo->Galaxy;
      ii=0;
      while(gal!=NULL){
        gal = gal->NextGalInHalo;
        counter++;
        ii++;
        if(ii>1000)
          ABORT(EXIT_FAILURE);
      }
      halo = halo->NextHaloInFOFGroup;
      halo_counter++;
      jj++;
      if (jj>1000)
        ABORT(EXIT_FAILURE);
    }
  }
  SID_log("Counting using FOF groups gives %d gals in %d halos", SID_LOG_COMMENT, counter, halo_counter);

  SID_log("...done", SID_LOG_CLOSE);
}
