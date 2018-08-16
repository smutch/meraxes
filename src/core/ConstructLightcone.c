#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

// DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>

//#include "XRayHeatingFunctions.c"

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
    
}

void ConstructLightcone(int snapshot)
{

    

}
