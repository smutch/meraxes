#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "meraxes.h"

double calc_metallicity(double total_gas, double metals)
{
  double Z;

  if((total_gas > 0) && (metals > 0))
    Z = metals / total_gas;
  else
    Z = 0.0;

  if(Z < 0)
    Z = 0.0;
  if(Z > 1)
    Z= 1.0;

  return Z;
}


int compare_ints(const void *a, const void *b)
{
  return *((int *)a) - *((int *)b);
}

void mpi_debug_here()
{
#ifdef DEBUG
  int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  printf("Once connected go up stack to 'sleep(5)' and 'set var i=7'\n");
  fflush(stdout);
  while (0 == i)
    sleep(5);
#endif
}

// static void find_missing_gals(run_globals_t *run_globals, fof_group_t *fof_group, int NFof, int flag)
// {
//   galaxy_t *gal = NULL;
//   halo_t *halo = NULL;
//   int counter = 0;
//   int master_counter = 0;
//   int missing_counter = 0;
//   bool *gal_found;
//   galaxy_t **missing_pointers;

//   SID_log("Running `find_missing_galaxies`...", SID_LOG_OPEN);

//   // If flag =0 then we have more galaxies in the global linked list than when
//   // we traverse the FOF groups.  If flag =1 then we have the opposite
//   // situation...

//   if(flag==0)
//   {
//     // Count the number of galaxies
//     gal = run_globals->FirstGal;
//     while (gal!=NULL)
//     {
//       gal->output_index = counter++;
//       gal = gal->Next;
//     }

//     gal_found = SID_calloc(sizeof(bool)*counter);

//     // Loop through each FOF halo and mark off each galaxy
//     for(int i_fof=0; i_fof<NFof; i_fof++)
//     {
//       halo = fof_group[i_fof].FirstHalo;
//       while (halo!=NULL) {
//         gal = halo->Galaxy;
//         while(gal!=NULL){
//           gal_found[gal->output_index] = true;
//           gal = gal->NextGalInHalo;
//         }
//         halo = halo->NextHaloInFOFGroup;
//       }
//     }
//   } else if(flag==1)
//   {
//     // Count the number of galaxies
//     for(int i_fof=0; i_fof<NFof; i_fof++)
//     {
//       halo = fof_group[i_fof].FirstHalo;
//       while (halo!=NULL) {
//         gal = halo->Galaxy;
//         while(gal!=NULL){
//           gal->output_index = counter++;
//           gal = gal->NextGalInHalo;
//         }
//         halo = halo->NextHaloInFOFGroup;
//       }
//     }

//     SID_log("I find counter=%d using FOF traversal...", SID_LOG_COMMENT, counter);

//     gal_found = SID_malloc(sizeof(bool)*counter);
//     for(int ii=0; ii<counter; ii++)
//       gal_found[ii] = false;

//     // Traverse the global linked list and mark off each galaxy
//     gal = run_globals->FirstGal;
//     master_counter = 0;
//     while (gal!=NULL)
//     {
//       if(!gal->ghost_flag)
//         gal_found[gal->output_index] = true;
//       master_counter++;
//       gal = gal->Next;
//     }
//     SID_log("I find %d gals traversing global list...", SID_LOG_COMMENT, master_counter);

//   }

//   // Now create an array which holds pointers to the missing galaxies
//   for(int ii=0; ii<counter; ii++)
//     if(!gal_found[ii])
//     {
//       SID_log("ii = %d", SID_LOG_COMMENT, ii);
//       missing_counter++;
//     }

//   master_counter = counter;

//   // Check the number of gals with ghost_flag=true
//   gal = run_globals->FirstGal;
//   counter = 0;
//   while(gal!=NULL)
//   {
//     if(gal->ghost_flag)
//       counter++;
//   gal = gal->Next;
//   }
//   SID_log("I find %d gals with ghost_flag=true", SID_LOG_COMMENT, counter);
//   counter = 0;
//   for(int i_fof=0; i_fof<NFof; i_fof++)
//   {
//     halo = fof_group[i_fof].FirstHalo;
//     while (halo!=NULL) {
//       gal = halo->Galaxy;
//       while(gal!=NULL){
//         if(gal->ghost_flag)
//           counter++;
//         gal = gal->NextGalInHalo;
//       }
//       halo = halo->NextHaloInFOFGroup;
//     }
//   }
//   SID_log("I find %d gals with ghost_flag=true (FOF traversal)", SID_LOG_COMMENT, counter);

//   missing_pointers = SID_calloc(sizeof(galaxy_t *)*missing_counter);

//   // Loop through the galaxies and store the pointers of the missing ones
//   counter = 0;
//   if(flag==0)
//   {
//     gal = run_globals->FirstGal;
//     while (gal!=NULL)
//     {
//       // Note that we only store non-ghost missing pointers here...
//       if((!gal_found[gal->output_index]) && (gal->SnapSkipCounter<=0))
//         missing_pointers[counter++] = gal;
//       gal = gal->Next;
//     }
//   } else if(flag==1)
//   {
//     for(int ii=0; ii<master_counter; ii++)
//       if(!gal_found[ii])
//       {
//         for(int i_fof=0; i_fof<NFof; i_fof++)
//         {
//           halo = fof_group[i_fof].FirstHalo;
//           while (halo!=NULL) {
//             gal = halo->Galaxy;
//             while(gal!=NULL){
//               if(gal->output_index==ii)
//                 missing_pointers[counter++] = gal;
//               gal = gal->NextGalInHalo;
//             }
//             halo = halo->NextHaloInFOFGroup;
//           }
//         }
//       }
//   }

//   mpi_debug_here();

//   SID_free(SID_FARG missing_pointers);
//   SID_free(SID_FARG gal_found);

//   SID_log("...done", SID_LOG_CLOSE);
// }


void check_counts(run_globals_t *run_globals, fof_group_t *fof_group, int NGal, int NFof)
{

  int counter = 0;
  int gal_next_counter = 0;
  int halo_counter = 0;
  int halo_pop_count = 0;
  int total_NGal = 0;
  int total_NFof = 0;
  int total_NGhosts = 0;
  galaxy_t *gal = NULL;
  halo_t *halo = NULL;

  SID_log("Running counts check...", SID_LOG_OPEN|SID_LOG_TIMER);

  SID_Allreduce(&NFof, &total_NFof, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_Allreduce(&NGal, &total_NGal, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_Allreduce(&(run_globals->NGhosts), &total_NGhosts, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_log("NFof = %d", SID_LOG_COMMENT, total_NFof);
  SID_log("NGal = %d", SID_LOG_COMMENT, total_NGal);
  SID_log("NGhosts = %d", SID_LOG_COMMENT, total_NGhosts);

  counter=0;
  gal = run_globals->FirstGal;
  while (gal!=NULL)
  {
    counter++;
    gal = gal->Next;
  }
  SID_Allreduce(SID_IN_PLACE, &counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_log("Counting using gal->Next gives %d gals (-%d ghosts = %d gals)",
      SID_LOG_COMMENT, counter, total_NGhosts,
      counter-total_NGhosts);
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
  SID_Allreduce(SID_IN_PLACE, &counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE, &halo_counter, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE, &halo_pop_count, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  SID_log("Counting using FOF groups gives %d gals in %d halos", SID_LOG_COMMENT, counter, halo_counter);
  SID_log("%d halos are populated with at least one galaxy", SID_LOG_COMMENT, halo_pop_count);

  // if((gal_next_counter - NGhosts) != counter)
  // {
  //   int flag;
  //   if((gal_next_counter - NGhosts) > counter)
  //     flag = 0;
  //   else
  //     flag = 1;
  //   find_missing_gals(run_globals, fof_group, NFof, flag);
  //   ABORT(EXIT_FAILURE);
  // }

  SID_log("...done", SID_LOG_CLOSE);
}

