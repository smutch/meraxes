#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

/*
 * This code computes the 21cm brightness temperature, which is somewhat
 * of a re-write of the first part of delta_T.c from 21cmFAST. Also includes
 * a prescription for line-of-sight redshift-space distortions, taken from
 * 21CMMC (Greig & Mesinger 2018). 
 * Written by Bradley Greig.
 */

void ComputeBrightnessTemperatureBox(int snapshot) {


    float* xH = run_globals.reion_grids.xH;
    float* deltax = run_globals.reion_grids.deltax;

    float* delta_T = run_globals.reion_grids.delta_T;



}
