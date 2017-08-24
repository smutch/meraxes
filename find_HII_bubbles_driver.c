#include "meraxes.h"
#include "meraxes_gpu.h"
#include "utils.h"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>

#include <hdf5.h>
#include <hdf5_hl.h>

void find_HII_bubbles_driver(
    int        snapshot,
    void       (*_find_HII_bubbles_passed)(double redshift),
    const char *reference_directory,
    timer_info *timer)
{
    // Get the snapshot & redshift for this output
    double redshift = run_globals.ZZ[snapshot];

    // Read input datasets
    char fname[STRLEN];
    sprintf(fname, "%s/validation_input-core%03d-z%.2f.h5",reference_directory,run_globals.mpi_rank,redshift);
    hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    float *deltax             = run_globals.reion_grids.deltax;
    float *stars              = run_globals.reion_grids.stars;
    float *sfr                = run_globals.reion_grids.sfr;
    float *z_at_ionization    = run_globals.reion_grids.z_at_ionization;
    float *J_21_at_ionization = run_globals.reion_grids.J_21_at_ionization;
    H5LTread_dataset_float(file_id, "deltax", deltax);
    H5LTread_dataset_float(file_id, "sfr",  sfr);
    H5LTread_dataset_float(file_id, "stars", stars);
    H5LTread_dataset_float(file_id, "z_at_ionization",    z_at_ionization);
    H5LTread_dataset_float(file_id, "J_21_at_ionization", J_21_at_ionization);
    H5Fclose(file_id);

    // Set a couple reionization parameters
    set_fesc(snapshot);
    set_ReionEfficiency();

    // Call the version of find_HII_bubbles we've been passed (and time it)
    timer_start(timer);
    _find_HII_bubbles_passed(redshift);
    timer_stop(timer);
    
}

