#include <math.h>
#include "meraxes.h"
#include "tree_flags.h"

static inline bool check_for_merger(int flags)
{
 if ((flags & TREE_CASE_MERGER)==TREE_CASE_MERGER)
   return true;
 else
   return false;
}

static inline bool check_if_valid_host(int flags)
{
  int invalid_flags = (TREE_CASE_FRAGMENTED_RETURNED
      | TREE_CASE_STRAYED
      | TREE_CASE_SPUTTERED);
  if ((flags & invalid_flags)==0)
    return true;
  else
    return false;
}

//! Actually run the model
void dracarys(run_globals_struct *run_globals)
{

  trees_header_struct   trees_header;
  halo_struct          *halo      = NULL;
  fof_group_struct     *fof_group = NULL;
  galaxy_struct        *gal       = NULL;
  galaxy_struct        *prev_gal  = NULL;
  galaxy_struct        *next_gal  = NULL;
  galaxy_struct        *cur_gal   = NULL;
  int                   i_newhalo;
  int                   NGal         = 0;
  int                   unique_ID    = 0;
  int                   nout_gals;
  int                   last_nout_gals;
  int                   last_snap    = 0;
  double                dt;
 
  // Find what the last requested output snapshot is
  for(int ii=0; ii<NOUT; ii++)
    if (run_globals->ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals->ListOutputSnaps[ii];

  // Loop through each snapshot
  for(int snapshot=0; snapshot<last_snap; snapshot++)
  {

    // Read in the halos for this snapshot
    trees_header = read_halos(run_globals, snapshot, &halo, &fof_group);
    gal      = run_globals->FirstGal;
    prev_gal = NULL;
    dt       = run_globals->LTTime[snapshot-1]-run_globals->LTTime[snapshot];
   
    SID_log("Processing snapshot %d...", SID_LOG_OPEN|SID_LOG_TIMER, snapshot);

    while (gal != NULL) {
      i_newhalo = gal->HaloDescIndex;

      if(i_newhalo>-1)
      {
        if(check_for_merger(gal->TreeFlags))
        {
          // Here we have a new merger...  Mark it and deal with it below.

          // If we have already marked this a type two it has already been
          // processed as a new merger and so we don't need to do this again...
          if (gal->Type != 2) 
            gal->Type = 999;

          gal->Halo = &(halo[i_newhalo]);
          // SID_log("Found a galaxy which now has no halo (merged into halo %d)", SID_LOG_COMMENT, i_newhalo);
        } else if(gal->Type < 2)
        {
          if (halo[i_newhalo].Galaxy == NULL)
            halo[i_newhalo].Galaxy = gal;
          else {
            SID_log("Trying to assign first galaxy to a halo which already has a first galaxy!", SID_LOG_COMMENT);
            ABORT(EXIT_FAILURE);
          }
          
          // Here we have the simplest case where a galaxy continues along in it's halo...
          gal->dM = (halo[i_newhalo]).Mvir - gal->Mvir;
          gal->dMdt = (gal->dM)/dt;

          copy_halo_to_galaxy(run_globals, &(halo[i_newhalo]), gal);

          // SID_log("Assigned existing galaxy to halo %d", SID_LOG_COMMENT, i_newhalo);
        }
      } else
      {
        // This galaxy is done (merged, lost, whatever...) so get rid of it
        if(prev_gal!=NULL)
          prev_gal->Next = gal->Next;
        else
          run_globals->FirstGal = gal->Next;
        cur_gal = gal->FirstGalInHalo;
        while ((cur_gal->NextGalInHalo != gal) && (cur_gal->NextGalInHalo != NULL))
          cur_gal = cur_gal->NextGalInHalo;
        cur_gal->NextGalInHalo = gal->NextGalInHalo;
        SID_free(SID_FARG gal);
        gal = prev_gal;
        NGal--;
        // SID_log("Killed a galaxy and decremented the counter.", SID_LOG_COMMENT);
      }

      // gal may be NULL if we just killed the first galaxy
      if (gal!=NULL)
      {
        prev_gal = gal;
        gal = gal->Next;
      } else
        gal = run_globals->FirstGal;
    }

    // Incase we ended up removing the last galaxy, update the LastGal pointer
    run_globals->LastGal = prev_gal;

    // Find empty (valid) type 0 halos and place new galaxies in them
    for(int i_halo=0; i_halo<trees_header.n_subgroups; i_halo++)
    {
      if((halo[i_halo].Type == 0) && (halo[i_halo].Galaxy == NULL) && check_if_valid_host(halo[i_halo].TreeFlags))
      {
        gal = new_galaxy(&unique_ID);
        copy_halo_to_galaxy(run_globals, &(halo[i_halo]), gal);
        if (run_globals->LastGal != NULL)
          run_globals->LastGal->Next = gal;
        else
          run_globals->FirstGal = gal;
        run_globals->LastGal = gal;
        halo[i_halo].Galaxy = gal;
        gal->FirstGalInHalo = gal;
        // SID_log("Created new galaxy in i_halo=%d", SID_LOG_COMMENT, i_halo);
        NGal++;
      }
    }

    // Loop through each galaxy and deal with mergers now that all other galaxies have been 
    // correctly propogated forwards
    gal = run_globals->FirstGal;
    while (gal != NULL) {
      if(gal->Type == 999)
      {
        if(gal->Halo->Galaxy == NULL)
        {
          // Here we have a halo with a galaxy that has just merged into an
          // empty halo.  From the point of view of the model, this isn't
          // actually a merger and so we need to catch these cases...
          gal->dM = gal->Halo->Mvir - gal->Mvir;
          gal->dMdt = (gal->dM)/dt;
          copy_halo_to_galaxy(run_globals, gal->Halo, gal);
          gal->Halo->Galaxy = gal;
          gal->FirstGalInHalo = gal;
        } else
        {
          // If there is a galaxy in the halo which is being merged into then
          // we actually have a bona fide merger...
          gal->Type = 2;
          cur_gal = gal->Halo->Galaxy;
          while (cur_gal!=NULL) {
            prev_gal = cur_gal;
            cur_gal = cur_gal->NextGalInHalo;
          }
          prev_gal->NextGalInHalo = gal;

          gal->FirstGalInHalo = gal->Halo->Galaxy;
          gal->MergerTarget = gal->FirstGalInHalo;
          // SID_log("Snap %d: Just set merger target of ID %d to ID %d...", SID_LOG_COMMENT, snapshot, gal->ID, gal->MergerTarget->ID);
          gal->MergTime = calculate_merging_time(run_globals, gal, snapshot);
        }
      }
      gal = gal->Next;
    }
    
    // DEBUG
    check_counts(run_globals, fof_group, NGal, trees_header.n_groups);

    // Do the physics
    nout_gals = evolve_galaxies(run_globals, fof_group, snapshot, NGal, trees_header.n_groups);

    // Write the results if this is a requested snapshot
    for(int i_out = 0; i_out < NOUT; i_out++)
      if(snapshot == run_globals->ListOutputSnaps[i_out])
        write_snapshot(run_globals, nout_gals, i_out, &last_nout_gals);
  
    SID_free(SID_FARG halo);
    SID_free(SID_FARG fof_group);

    SID_log("...done", SID_LOG_CLOSE);
  }

  // Free all of the remaining allocated galaxies
  gal = run_globals->FirstGal;
  while (gal != NULL) {
    next_gal = gal->Next;
    SID_free(SID_FARG gal);
    gal = next_gal;
  }

}


