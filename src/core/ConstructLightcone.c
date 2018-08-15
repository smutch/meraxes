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

    double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc
    int ReionGridDim = run_globals.params.ReionGridDim;

    double dR = (box_size / (double)ReionGridDim) * MPC; // cell size in comoving cm

    double start_redshift = run_globals.ZZ[0];
    double z1_LC, z2_LC, z_LC;

    z1_LC = z_LC = run_globals.params.End_Lightcone;

    long long total_slice_i = 0;
    i = 0;

    while(z1_LC < start_redshift) {
  
        z2_LC = run_globals.ZZ[i];      
        
        while( z_LC < z2_LC ) {

            total_slice_i++;
            z_LC -= dR / drdz(z_LC);
        }
        z1_LC = z2_LC;
        i += 1; 
    }
    
    run_globals.params.LightconeLength = total_slice_i;
    
}

void ConstructLightcone(int snapshot)
{

    

}
