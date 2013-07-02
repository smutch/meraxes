#include <math.h>
#include "meraxes.h"
#include "tree_flags.h"

static inline bool check_for_flag(int flag, int tree_flags)
{
  if ((tree_flags & flag)==flag)
    return true;
  else
    return false;
}

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

  trees_header_struct  trees_header;
  halo_struct         *halo            = NULL;
  fof_group_struct    *fof_group       = NULL;
  galaxy_struct       *gal             = NULL;
  galaxy_struct       *prev_gal        = NULL;
  galaxy_struct       *next_gal        = NULL;
  galaxy_struct       *cur_gal         = NULL;
  int                  i_newhalo;
  int                  NGal            = 0;
  int                  unique_ID       = 0;
  int                  nout_gals;
  int                  last_nout_gals;
  int                  last_snap       = 0;
  double               dt;
  int                  kill_counter    = 0;
  int                  merger_counter  = 0;
  int                  new_gal_counter = 0;
  int                  ghost_counter   = 0;
 
  // Find what the last requested output snapshot is
  for(int ii=0; ii<NOUT; ii++)
    if (run_globals->ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals->ListOutputSnaps[ii];

  // Loop through each snapshot
  for(int snapshot=0; snapshot<=last_snap; snapshot++)
  {

    // Reset book keeping counters
    kill_counter = 0;
    merger_counter = 0;
    new_gal_counter = 0;
    ghost_counter = 0;

    // Read in the halos for this snapshot
    trees_header = read_halos(run_globals, snapshot, &halo, &fof_group);

    // TODO: This should be dependant on the number of snapshots which the galaxy's halo has skipped
    dt       = run_globals->LTTime[snapshot-1]-run_globals->LTTime[snapshot];
   
    SID_log("Processing snapshot %d...", SID_LOG_OPEN|SID_LOG_TIMER, snapshot);

    // Reset the ghost flag on all galaxies
    gal      = run_globals->FirstGal;
    while(gal!=NULL)
    {
      gal->ghost_flag = false;
      gal = gal->Next;
    }

    gal      = run_globals->FirstGal;
    prev_gal = NULL;
    while (gal != NULL) {
      i_newhalo = gal->HaloDescIndex;
      gal->SnapSkipCounter--;  // Decrement the skip counter

      if(gal->SnapSkipCounter<=0)
      {
        if(i_newhalo>-1)
        {
          if(check_for_merger(gal->TreeFlags))
          {
            // Here we have a new merger...  Mark it and deal with it below.

            // If we have already marked this a type two it has already been
            // processed as a new merger and so we don't need to do this again...
            if (gal->Type != 2) 
            {
              gal->Type = 999;
              merger_counter++;
              gal->Halo = &(halo[i_newhalo]);
            }

            // TODO: Would be better just to turn off the gal->TreeFlags MERGER flag at this point...

            // SID_log("Found a galaxy which now has no halo (merged into halo %d)", SID_LOG_COMMENT, i_newhalo);
          } else if(gal->Type < 2)
          {
            // Here we have the simplest case where a galaxy continues along in it's halo...
            gal->dM = (halo[i_newhalo]).Mvir - gal->Mvir;
            gal->dMdt = (gal->dM)/dt;

            copy_halo_to_galaxy(run_globals, &(halo[i_newhalo]), gal);

            // Loop through all of the other galaxies in this halo and set their Halo pointer
            cur_gal = gal->NextGalInHalo;
            while(cur_gal!=NULL)
            {
              cur_gal->Halo = &(halo[i_newhalo]);
              cur_gal = cur_gal->NextGalInHalo;
            }

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

          if (cur_gal != gal)
          {
            while ((cur_gal->NextGalInHalo != gal) && (cur_gal->NextGalInHalo != NULL))
              cur_gal = cur_gal->NextGalInHalo;
            cur_gal->NextGalInHalo = gal->NextGalInHalo;
          } else
          {
            // DEBUG
            SID_log("We have just killed the first galaxy in a halo...(TreeFlags=%d)", SID_LOG_COMMENT, gal->TreeFlags);

            // We have just killed the first galaxy in the halo. If there are any
            // other galaxies left then we must kill them as well...
            // Unfortunately, we don't know if we have alrady processed these
            // remaining galaxies, so we have to just mark them for the moment
            // and then do a check below...
            cur_gal = cur_gal->NextGalInHalo;
            while(cur_gal!=NULL)
            {
              // Note that we mark these galaxies using HaloDescIndex=-999 here.
              // When we check for this below, we want to be able to ignore
              // halos which have no descendant at snapshot = snapshot+1
              cur_gal->HaloDescIndex = -999;
              cur_gal = cur_gal->NextGalInHalo;
            }
          }

          if (prev_gal == NULL)
            run_globals->FirstGal = gal->Next;

          SID_free(SID_FARG gal);
          gal = prev_gal;
          NGal--;
          kill_counter++;
        }
      } else
      {
        // This is a ghost galaxy for this snapshot
        // We need to count all the other galaxies in this halo as ghosts as
        // well since they won't be reachable by traversing the FOF groups
        cur_gal = gal;
        while (cur_gal!=NULL)
        {
          if(cur_gal->HaloDescIndex > -1)
          {
            ghost_counter++;
            cur_gal->ghost_flag = true;
          }
          cur_gal = cur_gal->NextGalInHalo;
        }
      }

      // gal may be NULL if we just killed the first galaxy
      if (gal!=NULL)
      {
        prev_gal = gal;
        gal = gal->Next;
      } else
      {
        //DEBUG
        SID_log("We just killed the first galaxy in the global linked list...", SID_LOG_COMMENT);

        gal = run_globals->FirstGal;
      }
    }

    // Do one more check to make sure that we have killed all galaxies which we
    // should have (i.e. satellites in strayed halos etc.)
    prev_gal = NULL;
    gal = run_globals->FirstGal;
    while(gal!=NULL)
    {
      if(gal->HaloDescIndex == -999)
      {
        // This galaxy is done (merged, lost, whatever...) so get rid of it
        if(prev_gal!=NULL)
          prev_gal->Next = gal->Next;
        else
          run_globals->FirstGal = gal->Next;
      }
      prev_gal = gal;
      gal = gal->Next;
    }

    // Store the number of ghost galaxies present at this snapshot
    run_globals->NGhosts = ghost_counter;

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
        gal->FirstGalInHalo = gal;
        NGal++;
        new_gal_counter++;
      }
    }

    SID_log("Identified %d new merger events.", SID_LOG_COMMENT, merger_counter);
    SID_log("Killed %d galaxies.", SID_LOG_COMMENT, kill_counter);
    SID_log("Created %d new galaxies.", SID_LOG_COMMENT, new_gal_counter);
    SID_log("There are %d galaxies in ghost halos (which are not being evolved).", SID_LOG_COMMENT, ghost_counter);

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
          cur_gal = gal->NextGalInHalo;
          while(cur_gal!=NULL)
          {
            cur_gal->Halo = gal->Halo;
            cur_gal = cur_gal->NextGalInHalo;
          }
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
          
          // Loop through and update the FirstGalInHalo pointers of any other
          // galaxies that are attached to this merged halo
          cur_gal = gal->NextGalInHalo;
          while(cur_gal!=NULL)
          {
            cur_gal->FirstGalInHalo = gal->FirstGalInHalo;
            cur_gal->Halo = gal->Halo;
            cur_gal = cur_gal->NextGalInHalo;
          }

          // DEBUG
          if (gal->FirstGalInHalo == NULL)
            SID_log("Just set gal->FirstGalInHalo = NULL!", SID_LOG_COMMENT);

          gal->MergerTarget = gal->FirstGalInHalo;
          gal->MergTime = calculate_merging_time(run_globals, gal, snapshot);

        }
      }
      gal = gal->Next;
    }
    
#ifdef DEBUG
    check_counts(run_globals, fof_group, NGal, trees_header.n_groups);
#endif

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


