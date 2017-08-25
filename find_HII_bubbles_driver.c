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
    const char *reference_directory,
    const bool  flag_write_validation_data,
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
    #ifdef __NVCC__
        #ifndef USE_CUFFT
            fprintf(stdout,"Calling hybrid-GPU/FFTW version of find_HII_bubbles() for snap=%d/z=%.2lf...",snapshot,redshift);fflush(stdout);
        #else
            fprintf(stdout,"Calling pure-GPU version of find_HII_bubbles() for snap=%d/z=%.2lf...",snapshot,redshift);fflush(stdout);
        #endif
        // Run the GPU version of _find_HII_bubbles()
        timer_start(timer);
        _find_HII_bubbles_gpu(redshift,flag_write_validation_data);
    #else
        // Run the Meraxes version of _find_HII_bubbles()
        fprintf(stdout,"Calling Meraxes version of find_HII_bubbles() for snap=%d/z=%.2lf...",snapshot,redshift);fflush(stdout);
        timer_start(timer);
        _find_HII_bubbles(redshift,flag_write_validation_data);
    #endif
    timer_stop(timer);
    fprintf(stdout,"Done. (%ld seconds)\n",timer_delta(*timer));fflush(stdout);

    // Write final output
    char fname_out[STRLEN];
    sprintf(fname_out,"validation_output-core%03d-z%.2f.h5",run_globals.mpi_rank, redshift);
    hid_t file_id_out = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const int slab_n_real = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]) * run_globals.params.ReionGridDim * run_globals.params.ReionGridDim;
    H5LTmake_dataset_float(file_id_out, "xH",                 1, (hsize_t []){slab_n_real}, run_globals.reion_grids.xH);
    H5LTmake_dataset_float(file_id_out, "z_at_ionization",    1, (hsize_t []){slab_n_real}, run_globals.reion_grids.z_at_ionization);
    H5LTmake_dataset_float(file_id_out, "J_21_at_ionization", 1, (hsize_t []){slab_n_real}, run_globals.reion_grids.J_21_at_ionization);
    H5LTset_attribute_double(file_id_out, "/", "volume_weighted_global_xH", &(run_globals.reion_grids.volume_weighted_global_xH), 1);
    H5LTset_attribute_double(file_id_out, "/", "mass_weighted_global_xH",   &(run_globals.reion_grids.mass_weighted_global_xH),   1);
    H5Fclose(file_id_out);
}

