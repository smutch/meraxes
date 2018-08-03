#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

/*
 * A generic function to compute the 21cm PS of any field. Algorithm taken from delta_T.c
 * from 21cmFAST, but generalised to take any input cubic field.
 * Written by Bradley Greig
 * 
 */

void Compute_PS(int snapshot)
{

    // A blank function for now.
    // Will need to be written in a way to be able to determine how to compute the power spectrum given the field
    // as some fields are already fluctuations, whereas others will need to be converted to a fluctuations quantity.
    // Fields of interest:
    // 1. Density field (already a fluctuation)
    // 2. Ionisation field (not a fluctuating quantity)
    // 3. 21cm brightness temperature (not a fluctuating quantity)

}
