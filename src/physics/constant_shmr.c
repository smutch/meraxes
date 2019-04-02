#include "meraxes.h"
#include <math.h>

//! Evolve existing galaxies forward in time
int constant_shmr(fof_group_t* fof_group, int snapshot, int NGal, int NFof)
{
    galaxy_t* gal = NULL;
    halo_t* halo = NULL;
    int gal_counter = 0;
    int dead_gals = 0;
    double infalling_gas = 0;
    double cooling_mass = 0;
    int NSteps = run_globals.params.NSteps;
    bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

    mlog("Using constant SHMR...", MLOG_OPEN | MLOG_TIMERSTART);

    for (int i_fof = 0; i_fof < NFof; i_fof++) {
        // First check to see if this FOF group is empty.  If it is then skip it.
        if (fof_group[i_fof].FirstOccupiedHalo == NULL)
            continue;
    
        for (int i_step = 0; i_step < NSteps; i_step++) {
            halo = fof_group[i_fof].FirstHalo;
            while (halo != NULL) {
                gal = halo->Galaxy;

                while (gal != NULL) {
                    if (gal->Type == 0) {
                        double old_stellar_mass = gal->StellarMass;
                        gal->StellarMass = pow(fof_group[i_fof].Mvir * run_globals.units.UnitMass_in_g / SOLAR_MASS, 1.4) * pow(10., -6.3) / run_globals.units.UnitMass_in_g * SOLAR_MASS;
                        
                        double new_stars = gal->StellarMass - old_stellar_mass;
                        if (new_stars < 0.0)
                            new_stars = 0.0;
                        gal->Sfr = new_stars / gal->dt;
                        gal->NewStars[0] = new_stars;
                    } else {
                        gal->StellarMass = 0.0;
                        gal->Sfr = 0.0;
                    }

                    // If this is a type 2 then decrement the merger clock
                    if (gal->Type == 2)
                        gal->MergTime -= gal->dt;

                    if (i_step == NSteps - 1)
                        gal_counter++;

                    gal = gal->NextGalInHalo;
                }

                halo = halo->NextHaloInFOFGroup;
            }

            // Check for mergers
            halo = fof_group[i_fof].FirstHalo;
            while (halo != NULL) {
                gal = halo->Galaxy;
                while (gal != NULL) {
                    if (gal->Type == 2)
                        // If the merger clock has run out or our target halo has already
                        // merged then process a merger event.
                        if ((gal->MergTime < 0) || (gal->MergerTarget->Type == 3)) {
                            galaxy_t* parent = gal->MergerTarget;
                            while (parent->Type == 3)
                                parent = parent->MergerTarget;
                            gal->MergerTarget = parent;

                            // Mark the merged galaxy as dead
                            gal->Type = 3;
                            gal->HaloDescIndex = -1;
                            dead_gals++;
                        }

                    // At this point, we've done all of the evolution for this galaxy,
                    // so we can update the escape fractions for stars and
                    // black holes based on the current properties and newly
                    // formed stars.
                    // calculate_galaxy_fesc_vals(gal, gal->NewStars[0], snapshot);

                    gal = gal->NextGalInHalo;
                }
                halo = halo->NextHaloInFOFGroup;
            }
        }
    }

    if (gal_counter + (run_globals.NGhosts) != NGal) {
        mlog_error("We have not processed the expected number of galaxies...");
        mlog("gal_counter = %d but NGal = %d", MLOG_MESG, gal_counter, NGal);
        ABORT(EXIT_FAILURE);
    }

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

    return gal_counter - dead_gals;
}
