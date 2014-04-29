#include "meraxes.h"
#include <math.h>

double calculate_merging_time(run_globals_t *run_globals, galaxy_t *sat, int snapshot)
{
  galaxy_t *parent;
  galaxy_t *mother;
  galaxy_t *cur_gal;
  double    coulomb;
  double    mergtime;
  double    sat_mass;
  double    sat_rad;

  // Note that we are assuming in this function that the halo properties
  // attached to the galaxies still correspond to the relevant values at the
  // last snapshot the two merging halos were last identified.  i.e. We are
  // assuming that the halo properties of the galaxies have *not* yet been
  // updated to the current snapshot (where the halos have already merged).

  // Find the merger "mother halo".  This is the most massive halo associated
  // with the merger event.  It's possible that there are >2 halos
  // participating in this merger but we want to use the most massive one in
  // the coulomb logarithm.
  cur_gal = sat->FirstGalInHalo;
  mother = cur_gal;
  while(cur_gal!=NULL)
  {
    if((cur_gal->OldType < 2) && (cur_gal->OldType > -1) && (cur_gal->Len > mother->Len))
      mother = cur_gal;
    cur_gal = cur_gal->NextGalInHalo;
  }

  parent = sat->MergerTarget;

  if(parent == sat)
  {
    SID_log_error("Invalid merger...!");
    ABORT(EXIT_FAILURE);
  }

  coulomb = log((double)(mother->Len) / (double)(sat->Len) + 1);

  // Include the baryonic mass in the calculation of the dynamical friction
  // timescale ala Guo+ 2010
  sat_mass = sat->Mvir + sat->ColdGas + sat->StellarMass;

  sat_rad = sqrt(
    pow(parent->Pos[0] - sat->Pos[0], 2.0) +
    pow(parent->Pos[1] - sat->Pos[1], 2.0) +
    pow(parent->Pos[2] - sat->Pos[2], 2.0) );

  // convert to physical length
  // Note that we want to use the redshift corresponding to the previous
  // snapshot (i.e. before the halo merged).  For cases where the halo has
  // skipped snapshots and then next been identified as having merged,
  // `snapshot-1` may not be correct.  However, we don't actually know when
  // during the time the skipped halo is missing from the trees that it last
  // existed unmerged, so `snapshot-1` is as good a time as any to pick.
  sat_rad /= (1 + run_globals->ZZ[snapshot-1]);

  if(sat_rad > mother->Rvir)
    sat_rad = mother->Rvir;

  if(sat_mass > 1e-9)
    mergtime =
    run_globals->params.MergerTimeFactor *
    1.17 * sat_rad * sat_rad * mother->Vvir / (coulomb * run_globals->G * sat_mass);
  else
    // if this is just a gas cloud (i.e. the infalling satellite has no stellar
    // mass), then instantly merge
    mergtime = -9999.9;

  return mergtime;

}


static void merger_driven_starburst(run_globals_t *run_globals, galaxy_t *parent, double merger_ratio, int snapshot)
{

  if(parent->ColdGas > 0)
  {
    // Calculate a merger driven starburst following Guo+ 2010
    double burst_mass;

    burst_mass = 0.56 * pow(merger_ratio, 0.7) * parent->ColdGas;

    // if(burst_mass > parent->ColdGas)
    //   burst_mass = parent->ColdGas;
    if(burst_mass < 0)
      burst_mass = 0.0;

    // apply the supernova feedback scheme and update the baryonic reservoirs
    supernova_feedback(run_globals, parent, burst_mass, merger_ratio, snapshot);
  }

}


void merge_with_target(run_globals_t *run_globals, galaxy_t *gal, int *dead_gals, int snapshot)
{

  galaxy_t *parent         = NULL;
  double    merger_ratio;
  double    parent_baryons;
  double    gal_baryons;

  // Identify the parent galaxy in the merger event.
  // Note that this relies on the merger target coming before this galaxy in
  // the linked list of halo members.  This should be the case but I should
  // confirm that it is always true...
  parent = gal->MergerTarget;
  while (parent->Type==3)
    parent = parent->MergerTarget;
  gal->MergerTarget = parent;

  // use the **baryonic** mass to calculate the merger ratio
  parent_baryons = parent->StellarMass + parent->ColdGas;
  gal_baryons    = gal->StellarMass + gal->ColdGas;
  if(parent_baryons > gal_baryons)
    merger_ratio = gal_baryons/parent_baryons;
  else
    merger_ratio = parent_baryons/gal_baryons;

  // Add galaxies together
  parent->StellarMass       += gal->StellarMass;
  parent->MetalsStellarMass += gal->MetalsStellarMass;
  parent->Sfr               += gal->Sfr;
  parent->HotGas            += gal->HotGas;
  parent->MetalsHotGas      += gal->MetalsHotGas;
  parent->ColdGas           += gal->ColdGas;
  parent->MetalsColdGas     += gal->MetalsColdGas;
  parent->EjectedGas        += gal->EjectedGas;
  parent->MetalsEjectedGas  += gal->MetalsEjectedGas;
  parent->BlackHoleMass     += gal->BlackHoleMass;

  for(int outputbin = 0; outputbin < NOUT; outputbin++)
    sum_luminosities(run_globals, parent, gal, outputbin);

  if(run_globals->params.physics.Flag_BHFeedback)
    merger_driven_BH_growth(run_globals, parent, merger_ratio);

  // merger driven starburst prescription
  merger_driven_starburst(run_globals, parent, merger_ratio, snapshot);

  // Mark the merged galaxy as dead
  gal->Type          = 3;
  gal->HaloDescIndex = -1;
  (*dead_gals)++;

}
