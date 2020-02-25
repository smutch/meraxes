#include <hdf5_hl.h>
#include <math.h>

#include "cooling.h"
#include "meraxes.h"

// *** This code is taken from the model of Croton+ 2006 with minimal modifications. ***
//
// Future todo list:
//  - Read in the below hard coded values from the HDF5 file itself and malloc appropriately
//  - Set up more flexible interpolation (e.g. without assumption of fixed temp step)
//  - Use GSL routines for interpolation

void read_cooling_functions()
{
    if (run_globals.mpi_rank == 0) {
        hid_t fd;
        char dset_name[30];
        char fname[STRLEN+11];

        sprintf(fname, "%s/SD93.hdf5", run_globals.params.CoolingFuncsDir);
        fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

        for (int i_m = 0; i_m < N_METALLICITIES; i_m++) {
            sprintf(dset_name, "%s/log(lambda_norm)", group_name[i_m]);
            H5LTread_dataset_double(fd, dset_name, cooling_rate[i_m]);
        }

        H5Fclose(fd);
    }

    // broadcast the values to all cores
    MPI_Bcast(&cooling_rate, sizeof(cooling_rate), MPI_BYTE, 0, run_globals.mpi_comm);

    // add solar metallicity to all metallicity values
    for (int i_m = 0; i_m < N_METALLICITIES; i_m++)
        metallicities[i_m] += log10(0.02);
}

static double interpolate_temp_dependant_cooling_rate(int i_m, double logTemp)
{
    int i_t;
    double temp_step;
    double rate, rate_below, rate_above;
    double logT_below;

    // First deal with boundary conditions
    if (logTemp < MIN_TEMP)
        return -27.0;

    // Now find the index of the tabulated temp immediately below our value
    temp_step = (MAX_TEMP - MIN_TEMP) / (double)(N_TEMPS - 1);
    i_t = (int)((logTemp - MIN_TEMP) / temp_step);
    if (i_t >= (N_TEMPS - 1))
        i_t = N_TEMPS - 2;

    // Now grab the cooling rates for the temp immediately above and below
    rate_below = cooling_rate[i_m][i_t];
    rate_above = cooling_rate[i_m][i_t + 1];

    // Calculate the tabulated temperature immediately below our value
    logT_below = MIN_TEMP + temp_step * i_t;

    // Now linearly interpolate the cooling rate
    rate = rate_below + (rate_above - rate_below) / temp_step * (logTemp - logT_below);

    return rate;
}

double interpolate_cooling_rate(double logTemp, double logZ)
{
    int i_m;
    double rate_below, rate_above, rate;

    // First deal with boundary conditions
    if (logZ < metallicities[0])
        logZ = metallicities[0];
    if (logZ > metallicities[N_METALLICITIES - 1])
        logZ = metallicities[N_METALLICITIES - 1];

    // Now find the indices of the metallicity values which bound our input
    i_m = 0;
    while (logZ > metallicities[i_m + 1])
        i_m++;

    // Get the cooling rates for this temperature value
    rate_below = interpolate_temp_dependant_cooling_rate(i_m, logTemp);
    rate_above = interpolate_temp_dependant_cooling_rate(i_m + 1, logTemp);

    // Finally, linearly interpolate the cooling rates
    rate = rate_below + (rate_above - rate_below) / (metallicities[i_m + 1] - metallicities[i_m]) * (logZ - metallicities[i_m]);

    return pow(10, rate);
}
