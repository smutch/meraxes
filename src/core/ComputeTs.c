#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

// DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>

#include "XRayHeatingFunctions.c"

/*
 * This code is a re-write of the spin temperature calculation (Ts.c) within 21cmFAST.
 * Modified for usage within Meraxes by Bradley Greig.
 */

void _ComputeTs(int snapshot)
{

    double box_size = run_globals.params.BoxSize; // Mpc/h
    int ReionGridDim = run_globals.params.ReionGridDim;
    double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
    double total_n_cells = pow((double)ReionGridDim, 3);
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
    int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
    int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
    double ReionEfficiency = run_globals.params.physics.ReionEfficiency;
    run_units_t* units = &(run_globals.units);

    double redshift = run_globals.ZZ[snapshot];
    double prev_redshift;
    if(snapshot==0) {
        prev_redshift = run_globals.ZZ[snapshot];
    }
    else {
        prev_redshift = run_globals.ZZ[snapshot-1];
    }

    int i_real, i_padded, R_ct, x_e_ct, n_ct, m_xHII_low, m_xHII_high, NO_LIGHT;

    double prev_zpp, prev_R, zpp, zp, lower_int_limit, filling_factor_of_HI_zp, R_factor, R, nuprime, dzp, Luminosity_converstion_factor;
    double collapse_fraction, total_SFR, density_over_mean;

    float curr_xalpha;

    double freq_int_heat[NUM_FILTER_STEPS_FOR_Ts], freq_int_ion[NUM_FILTER_STEPS_FOR_Ts], freq_int_lya[NUM_FILTER_STEPS_FOR_Ts];
    double freq_int_heat_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];
    double R_values[NUM_FILTER_STEPS_FOR_Ts];

    double *evolve_ans, ans[2], dansdz[5], xHII_call;
    double curr_delNL0[NUM_FILTER_STEPS_FOR_Ts];

    float* deltax = run_globals.reion_grids.deltax;

    float* x_e_box = run_globals.reion_grids.x_e_box;
    float* x_e_box_prev = run_globals.reion_grids.x_e_box_prev;
    float* Tk_box = run_globals.reion_grids.Tk_box;
    float* Tk_box_prev = run_globals.reion_grids.Tk_box_prev;
    float* TS_box = run_globals.reion_grids.TS_box;

    fftwf_complex* deltax_unfiltered = (fftwf_complex*)deltax; // WATCH OUT!
    fftwf_complex* deltax_filtered = run_globals.reion_grids.deltax_filtered;
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax, deltax_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float* stars = run_globals.reion_grids.stars;
    fftwf_complex* stars_unfiltered = (fftwf_complex*)stars; // WATCH OUT!
    fftwf_complex* stars_filtered = run_globals.reion_grids.stars_filtered;
    plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars, stars_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float* sfr = run_globals.reion_grids.sfr;
    fftwf_complex* sfr_unfiltered = (fftwf_complex*)sfr; // WATCH OUT!
    fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
    plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr, sfr_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    // TODO: Double check that looping over correct number of elements here
    for (int ii = 0; ii < slab_n_complex; ii++) {
        deltax_unfiltered[ii] /= total_n_cells;
        stars_unfiltered[ii] /= total_n_cells;
        sfr_unfiltered[ii] /= total_n_cells;
    }

    int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

    double** delNL0 = (double**)calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double*));
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        delNL0[R_ct] = (double*)calloc(total_n_cells, sizeof(double));
    }    


    // Initialise the RECFAST, electron rate tables
    init_heat();

    x_e_ave = 0.0;

    

    // Place current redshift in 21cmFAST nomenclature (zp), delta zp (dzp) and delta z in seconds (dt_dzp)
    zp = redshift;
    dzp = zp - prev_redshift;
    dt_dzp = dtdz(zp);


    // Check redshift against Z_HEAT_MAX. If zp > Z_HEAT_MAX assume the x_e (electron fraction) and gas temperatures are homogenous 
    // Equivalent to the default setup of 21cmFAST. 
    if(zp >= run_globals.params.physics.Z_HEAT_MAX) {

        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                    x_e_box_prev[i_real] = xion_RECFAST(zp,0);
                    Tk_box_prev[i_real] = T_RECFAST(zp,0);

                    TS_box[i_real] = get_Ts(zp, deltax[i_real], Tk_box_prev[i_real], x_e_box_prev[i_real],1, &curr_xalpha);

                }
    }



    double Ave_Ts = 0.0;
    double Ave_x_e = 0.0;
    double Ave_Tk = 0.0;

    for (int ix = 0; ix < local_nix; ix++)
        for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
                i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                Ave_Ts += (double)TS_box[i_real];
                Ave_Tk += (double)Tk_box_prev[i_real];
                Ave_x_e += (double)x_e_box_prev[i_real];
            }

    MPI_Allreduce(MPI_IN_PLACE, &Ave_Ts, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Ave_Tk, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Ave_x_e, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

    Ave_Ts /= total_n_cells;
    Ave_Tk /= total_n_cells;
    Ave_x_e /= total_n_cells;

    mlog("zp = %e Ts_ave = %e Tk_ave = %e x_e_ave = %e", MLOG_MESG, zp, Ave_Ts, Ave_Tk, Ave_x_e);

}

// This function makes sure that the right version of ComputeTs() gets called.
// Note: Only the CPU version works for now
void ComputeTs(int snapshot, timer_info* timer_total)
{
    // Call the version of ComputeTs we've been passed (and time it)
    int flag_write_validation_data = false;

    timer_info timer;
#ifdef USE_CUDA
#ifdef USE_CUFFT
    mlog("Calling pure-GPU version of ComputeTs() for snap=%d/z=%.2lf...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, redshift);
#else
    mlog("Calling hybrid-GPU/FFTW version of ComputeTs() for snap=%d/z=%.2lf...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, redshift);
#endif
    // Run the GPU version of _ComputeTs()
    timer_start(&timer);
//    _ComputeTs_gpu(redshift, flag_write_validation_data);
#else
    // Run the Meraxes version of _ComputeTs()
    mlog("Calling pure-CPU version of ComputeTs() for snap=%d/z=%.2lf...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, run_globals.ZZ[snapshot]);
    timer_start(&timer);
    _ComputeTs(snapshot);
#endif
    timer_stop(&timer);
    timer_stop(timer_total);
    timer_gpu += timer_delta(timer);
    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
    mlog("Total time spent in ComputeTs vs. total run time (snapshot %d ): %.2f of %.2f s", MLOG_MESG, snapshot, timer_gpu, timer_delta(*timer_total));
}
