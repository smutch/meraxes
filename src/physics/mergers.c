#include <math.h>

#include "core/misc_tools.h"
#include "blackhole_feedback.h"
#include "mergers.h"
#include "star_formation.h"
#include "supernova_feedback.h"
#include "core/magnitudes.h"
#include "meraxes.h"

double calculate_merging_time(galaxy_t* orphan, int snapshot)
{
    // TODO: What should we do about FOF properties here?  Do we need to use the
    // FOF virial properties for the merger clock calculation if one of the
    // participating galaxies was in a central subhalo?

    galaxy_t* parent = NULL;
    galaxy_t* mother = NULL;
    galaxy_t* sat = NULL;
    galaxy_t* cur_gal = NULL;
    double coulomb;
    double mergtime;
    double sat_mass;
    double sat_rad;
    double min_stellar_mass;

    // Note that we are assuming in this function that the halo properties
    // attached to the galaxies still correspond to the relevant values at the
    // last snapshot the two merging halos were last identified.  i.e. We are
    // assuming that the halo properties of the galaxies have *not* yet been
    // updated to the current snapshot (where the halos have already merged).

    // First check to see if the baryonic mass of the satellite or the parent is
    // zero (this can happen due to reionization).  If true, then just set the
    // merger clock to zero and finish.
    if ((orphan->StellarMass + orphan->ColdGas) < 2e-10)
        return -999;
    if ((orphan->MergerTarget->StellarMass + orphan->MergerTarget->ColdGas) < 2e-10)
        return -999;

    // Next we need to decide, for the purposes of calculating the merger
    // timescale, which galaxy/halo is the parent and which is the satellite.
    // Depending on the construction of the trees, we are not gauranteed that the
    // satellite halo will be more massive thatn the parent halo.  This could
    // introduce crazy merger timescales and so we must explicitly check...
    sat = orphan;
    parent = sat->MergerTarget;
    if (sat->Len > parent->Len) {
        parent = orphan;
        sat = orphan->MergerTarget;
    }

    if (parent == sat) {
        mlog_error("Invalid merger...!");
        ABORT(EXIT_FAILURE);
    }

    min_stellar_mass = (sat->StellarMass <= parent->StellarMass) ? sat->StellarMass : parent->StellarMass;
    if (min_stellar_mass < run_globals.params.physics.MinMergerStellarMass)
        return -999;

    // Find the merger "mother halo".  This is the most massive halo associated
    // with the merger event.  It's possible that there are >2 halos
    // participating in this merger but we want to use the most massive one in
    // the coulomb logarithm.
    cur_gal = sat->FirstGalInHalo;
    mother = cur_gal;
    while (cur_gal != NULL) {
        if ((cur_gal->OldType < 2) && (cur_gal->OldType > -1) && (cur_gal->Len > mother->Len))
            mother = cur_gal;
        cur_gal = cur_gal->NextGalInHalo;
    }

    coulomb = log1p((double)(mother->Len) / (double)(sat->Len));

    // Include the baryonic mass in the calculation of the dynamical friction
    // timescale ala Guo+ 2011.
    // N.B. Hot gas should be zero in the satellite anyway for >99% of cases
    // sat_mass = sat->Mvir + sat->ColdGas + sat->StellarMass + sat->BlackHoleMass + sat->HotGas;
    //
    // OR
    //
    // Don't explicitly include the baryon mass reservoirs for the satellite
    // mass.  The N-body sim used to generate the trees includes the full matter
    // density of the Universe. In other words - the virial mass should already
    // include the baryon mass in it.
    sat_mass = sat->Mvir;

    // TODO: Should this be parent or mother???
    sat_rad = (double)comoving_distance(mother->Pos, sat->Pos);

    // convert to physical length
    // Note that we want to use the redshift corresponding to the previous
    // snapshot (i.e. before the halo merged).  For cases where the halo has
    // skipped snapshots and then next been identified as having merged,
    // `snapshot-1` may not be correct.  However, we don't actually know when
    // during the time the skipped halo is missing from the trees that it last
    // existed unmerged, so `snapshot-1` is as good a time as any to pick.
    sat_rad /= (1 + run_globals.ZZ[snapshot - 1]);

    // TODO: Should this be parent or mother???
    orphan->MergerStartRadius = sat_rad / mother->Rvir;

    if (sat_rad > mother->Rvir)
        sat_rad = mother->Rvir;

    mergtime = run_globals.params.physics.MergerTimeFactor * 1.17 * sat_rad * sat_rad * mother->Vvir / (coulomb * run_globals.G * sat_mass);

    return mergtime;
}

static void merger_driven_starburst(galaxy_t* parent, double merger_ratio, int snapshot)
{
    if ((parent->ColdGas > 0) && (merger_ratio > run_globals.params.physics.MinMergerRatioForBurst)) {
        // Calculate a merger driven starburst following Guo+ 2010
        physics_params_t* params = &(run_globals.params.physics);

        double burst_mass = params->MergerBurstFactor * pow(merger_ratio, params->MergerBurstScaling) * parent->ColdGas;

        if (burst_mass > parent->ColdGas)
            burst_mass = parent->ColdGas;

        if (burst_mass > 0) {
            double m_reheat;
            double m_eject;
            double m_recycled;
            double new_metals;

            contemporaneous_supernova_feedback(parent, &burst_mass, snapshot, &m_reheat, &m_eject, &m_recycled, &new_metals);
            // update the baryonic reservoirs (note that the order we do this in will change the result!)
            update_reservoirs_from_sf(parent, burst_mass, snapshot, MERGER);
            parent->MergerBurstMass += burst_mass;
            update_reservoirs_from_sn_feedback(parent, m_reheat, m_eject, m_recycled, new_metals);
        }
    }
}

void merge_with_target(galaxy_t* gal, int* dead_gals, int snapshot)
{
    galaxy_t* parent = NULL;
    double merger_ratio;
    double parent_baryons;
    double gal_baryons;
    double min_stellar_mass;

    // Identify the parent galaxy in the merger event.
    // Note that this relies on the merger target coming before this galaxy in
    // the linked list of halo members.  This should be the case but I should
    // confirm that it is always true...
    parent = gal->MergerTarget;
    while (parent->Type == 3)
        parent = parent->MergerTarget;
    gal->MergerTarget = parent;

    // use the **baryonic** mass to calculate the merger ratio
    parent_baryons = parent->StellarMass + parent->ColdGas;
    gal_baryons = gal->StellarMass + gal->ColdGas;
    if (parent_baryons > gal_baryons)
        merger_ratio = gal_baryons / parent_baryons;
    else
        merger_ratio = parent_baryons / gal_baryons;

    min_stellar_mass = (gal->StellarMass <= parent->StellarMass) ? gal->StellarMass : parent->StellarMass;

    // Add galaxies together
    parent->StellarMass += gal->StellarMass;
    parent->GrossStellarMass += gal->GrossStellarMass;
    parent->FescWeightedGSM += gal->FescWeightedGSM;
    parent->MetalsStellarMass += gal->MetalsStellarMass;
    parent->Sfr += gal->Sfr;
    parent->HotGas += gal->HotGas;
    parent->MetalsHotGas += gal->MetalsHotGas;
    parent->ColdGas += gal->ColdGas;
    parent->MetalsColdGas += gal->MetalsColdGas;
    parent->EjectedGas += gal->EjectedGas;
    parent->MetalsEjectedGas += gal->MetalsEjectedGas;
    parent->BlackHoleAccretedHotMass += gal->BlackHoleAccretedHotMass;
    parent->BlackHoleAccretedColdMass += gal->BlackHoleAccretedColdMass;
    parent->BlackHoleAccretingColdMass += gal->BlackHoleAccretingColdMass;
    parent->BHemissivity += gal->BHemissivity;
    parent->BlackHoleMass += gal->BlackHoleMass;
    parent->EffectiveBHM += gal->EffectiveBHM;
    parent->mwmsa_num += gal->mwmsa_num;
    parent->mwmsa_denom += gal->mwmsa_denom;
    parent->MergerBurstMass += gal->MergerBurstMass;

    for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
        parent->NewStars[ii] += gal->NewStars[ii];

    for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
        parent->NewMetals[ii] += gal->NewMetals[ii];

#ifdef CALC_MAGS
    merge_luminosities(parent, gal);
#endif

    // merger driven starburst prescription
    if (min_stellar_mass >= run_globals.params.physics.MinMergerStellarMass)
        merger_driven_starburst(parent, merger_ratio, snapshot);

    // TODO: Should this have a stellar mass / baryon limit placed on it?
    if (run_globals.params.physics.Flag_BHFeedback)
        merger_driven_BH_growth(parent, merger_ratio, snapshot);

    // Mark the merged galaxy as dead
    gal->Type = 3;
    gal->HaloDescIndex = -1;
    (*dead_gals)++;
}
