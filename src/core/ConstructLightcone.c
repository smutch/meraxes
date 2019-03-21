#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

// DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>

/*
 * This code is a re-write of the light-cone (cuboid) construction from 21cmFAST.
 * Modified for usage within Meraxes by Bradley Greig.
 */

void Initialise_ConstructLightcone()
{
    int i;
    int closest_snapshot;

    double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc
    int ReionGridDim = run_globals.params.ReionGridDim;

    double dR = (box_size / (double)ReionGridDim) * MPC; // cell size in comoving cm

    double z1_LC, z2_LC, z_LC;
    double start_redshift = run_globals.ZZ[0];

    long long total_slice_i = 0;
    i = 0;

    // Find the first snapshot beyond the user provided beginning of the light-cone (defined from the lowest redshift)
    while(run_globals.ZZ[i] > run_globals.params.End_Lightcone) {
        i++;
    }
    closest_snapshot = i;

    // Determine if the beginning of the light-cone (lower redshift) is lower than the final snapshot provided by the snapshot list
    // If it is, set the lowest redshift to be the final snapshot provided by the snapshot list
    if(closest_snapshot > run_globals.LastOutputSnap) {
        closest_snapshot = run_globals.LastOutputSnap;
    }

    mlog("final snapshot for light-cone = %d",MLOG_MESG,closest_snapshot);

    // End of the light-cone, i.e. low redshift part
    z1_LC = z_LC = run_globals.ZZ[closest_snapshot];

    i = 0;
    // Below incrementes backwards to match 21cmFAST
    while(z1_LC < start_redshift) {

        z2_LC = run_globals.ZZ[closest_snapshot-1-i];

        while( z_LC < z2_LC ) {

            total_slice_i++;
            z_LC -= dR / drdz(z_LC);
        }
        z1_LC = z2_LC;
        i += 1;
    }

    // Store the length of the light-cone to be able to allocate the array to hold the light-cone
    run_globals.params.LightconeLength = total_slice_i;
    run_globals.params.End_Lightcone_snapshot = closest_snapshot;

}

void ConstructLightcone(int snapshot)
{
    int iz;

    float* LightconeBox = run_globals.reion_grids.LightconeBox;

    float* delta_T = run_globals.reion_grids.delta_T;
    float* delta_T_prev = run_globals.reion_grids.delta_T_prev;

    float* Lightcone_redshifts = run_globals.reion_grids.Lightcone_redshifts;

    int ReionGridDim = run_globals.params.ReionGridDim;
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
    int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
    int slab_n_real_LC = local_nix * ReionGridDim * run_globals.params.LightconeLength;

    double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc
    double dR = (box_size / (double)ReionGridDim) * MPC; // cell size in comoving cm

    double z_LC, z1_LC, z2_LC, t_z1_LC, t_z2_LC, z_slice, t_z_slice, fz1, fz2, fz, stored_z_LC;

    long long slice_ct = 0;
    long long slice_ct_snapshot = 0;

    int i_real, i_real_LC, i, closest_snapshot;

    if(snapshot>0) {

        // Set the redshift (and time) for the upper redshift end of the light-cone
        z2_LC = run_globals.ZZ[snapshot-1];
        t_z2_LC = gettime(z2_LC);

        // Set the redshift (and time) for the lower redshift end of the light-cone
        z1_LC = z_LC  = run_globals.ZZ[snapshot];
        t_z1_LC = gettime(z1_LC);

        // Repeat some of the light-cone initialisation just to avoid passing things around

        // Set the lowest redshift snapshot
        closest_snapshot = run_globals.params.End_Lightcone_snapshot;

        z_LC = run_globals.ZZ[closest_snapshot];

        // Loop through to find the size of lightcone added at this time step (solely used for indexing the light-cone properly across snapshots)
        while (z_LC < z2_LC) {

            if(z_LC >= run_globals.ZZ[snapshot] && z_LC < run_globals.ZZ[snapshot-1]) {
                slice_ct_snapshot++;
            }
            z_LC -= dR / drdz(z_LC);
        }



        z_LC = run_globals.ZZ[closest_snapshot];

        // Determine the starting indexes for this co-eval snapshot of the light-cone (incrementing from lower to upper redshifts)
        run_globals.params.CurrentLCPos = run_globals.params.CurrentLCPos - slice_ct_snapshot;
        slice_ct = run_globals.params.CurrentLCPos;
        iz = slice_ct % ReionGridDim;

        // Now do the interpolation of the light-cone
        while (z_LC < z2_LC) {

            if(z_LC >= run_globals.ZZ[snapshot] && z_LC < run_globals.ZZ[snapshot-1]) {
                z_slice = z_LC;
                t_z_slice = gettime(z_slice);

                Lightcone_redshifts[slice_ct] = z_slice;

                // Ensure we don't overstep the gridsize of the co-eval box
                if(iz >= ReionGridDim) {
                    iz = 0;
                }

                for (int ii = 0; ii < local_nix; ii++) {
                    for (int jj=0;jj<ReionGridDim; jj++){
                        i_real_LC =  grid_index_LC(ii, jj, slice_ct, ReionGridDim, run_globals.params.LightconeLength);
                        i_real = grid_index(ii, jj, iz, ReionGridDim, INDEX_REAL);

                        fz1 = delta_T[i_real];
                        fz2 = delta_T_prev[i_real];
                        run_globals.reion_grids.LightconeBox[i_real_LC] = (fz2 - fz1) / (t_z2_LC - t_z1_LC) * (t_z_slice - t_z1_LC) + fz1; // linearly interpolate in z (time actually)
                    }
                }

                iz++;
                slice_ct++;
            }
            z_LC -= dR / drdz(z_LC);
        }
    }
    else {
        // Used to correctly index the starting point of the co-eval boxes for the light-cone
        run_globals.params.CurrentLCPos = run_globals.params.LightconeLength;
    }

    // Update the previous delta_T box with the one we just finished using
    memcpy(delta_T_prev, delta_T, sizeof(float) * slab_n_real);

}

