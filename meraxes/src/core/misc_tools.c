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

static void find_missing_gals(run_globals_struct *run_globals, fof_group_struct *fof_group, int NFof)
{
  galaxy_struct *gal = NULL;
  halo_struct *halo = NULL;
  int counter = 0;
  int missing_counter = 0;
  bool *gal_found;
  galaxy_struct **missing_pointers;

  SID_log("Running `find_missing_galaxies`...", SID_LOG_OPEN);

  // Count the number of galaxies
  gal = run_globals->FirstGal;
  while (gal!=NULL)
  {
    gal->output_index = counter++;
    gal = gal->Next;
  }

  gal_found = SID_calloc(sizeof(bool)*counter);

  // Loop through each FOF halo and mark off each galaxy
  for(int i_fof=0; i_fof<NFof; i_fof++)
  {
    halo = fof_group[i_fof].FirstHalo;
    while (halo!=NULL) {
      gal = halo->Galaxy;
      while(gal!=NULL){
        gal_found[gal->output_index] = true;
        gal = gal->NextGalInHalo;
      }
      halo = halo->NextHaloInFOFGroup;
    }
  }
 
  // Now create an array which holds pointers to the missing galaxies
  for(int ii=0; ii<counter; ii++)
    if(!gal_found[ii])
      missing_counter++;

  missing_pointers = SID_calloc(sizeof(galaxy_struct *)*missing_counter);

  // Loop through the galaxies and store the pointers of the missing ones
  gal = run_globals->FirstGal;
  counter = 0;
  while (gal!=NULL)
  {
    // Note that we only store non-ghost missing pointers here...
    if((!gal_found[gal->output_index]) && (gal->SnapSkipCounter<=0))
      missing_pointers[counter++] = gal;
    gal = gal->Next;
  }

  mpi_debug_here();

  SID_free(SID_FARG missing_pointers);
  SID_free(SID_FARG gal_found);

  SID_log("...done", SID_LOG_CLOSE);
}


void check_counts(run_globals_struct *run_globals, fof_group_struct *fof_group, int NGal, int NFof)
{

  int counter = 0;
  int gal_next_counter = 0;
  int halo_counter = 0;
  int halo_pop_count = 0;
  galaxy_struct *gal = NULL;
  halo_struct *halo = NULL;

  SID_log("Running counts check...", SID_LOG_OPEN|SID_LOG_TIMER);

  SID_log("NFof = %d", SID_LOG_COMMENT, NFof);
  SID_log("NGal = %d", SID_LOG_COMMENT, NGal);
  SID_log("NGhosts = %d", SID_LOG_COMMENT, run_globals->NGhosts);

  counter=0;
  gal = run_globals->FirstGal;
  while (gal!=NULL)
  {
    counter++;
    gal = gal->Next;
  }
  SID_log("Counting using gal->Next gives %d gals (-%d ghosts = %d gals)",
      SID_LOG_COMMENT, counter, run_globals->NGhosts,
      counter-run_globals->NGhosts);
  gal_next_counter = counter;

  halo_pop_count = 0;
  counter = 0;
  halo_counter = 0;
  int ii, jj;
  for(int i_fof=0; i_fof<NFof; i_fof++)
  {
    halo = fof_group[i_fof].FirstHalo;
    jj=0;
    while (halo!=NULL) {
      gal = halo->Galaxy;
      if (gal!=NULL)
        halo_pop_count++;
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
  SID_log("%d halos are populated with at least one galaxy", SID_LOG_COMMENT, halo_pop_count);

  if(gal_next_counter-(run_globals->NGhosts) != counter)
  {
#ifdef DEBUG
    find_missing_gals(run_globals, fof_group, NFof);
#endif
    ABORT(EXIT_FAILURE);
  }

  SID_log("...done", SID_LOG_CLOSE);
}

