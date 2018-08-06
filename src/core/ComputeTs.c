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

    double prev_zpp, prev_R, zpp, zp, lower_int_limit_GAL, lower_int_limit_QSO, filling_factor_of_HI_zp, R_factor, R, nuprime, dzp, Luminosity_converstion_factor;
    double collapse_fraction, total_SFR, density_over_mean;

    float curr_xalpha;

    double freq_int_heat_GAL[NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_GAL[NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_GAL[NUM_FILTER_STEPS_FOR_Ts];
    double freq_int_heat_QSO[NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_QSO[NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_QSO[NUM_FILTER_STEPS_FOR_Ts];

    double freq_int_heat_tbl_GAL[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl_GAL[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl_GAL[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];
    double freq_int_heat_tbl_QSO[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl_QSO[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl_QSO[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];

    double R_values[NUM_FILTER_STEPS_FOR_Ts];

    double *evolve_ans, ans[2], dansdz[5], xHII_call;
    double curr_delNL0[NUM_FILTER_STEPS_FOR_Ts];

    float* deltax = run_globals.reion_grids.deltax;
    float* stars = run_globals.reion_grids.stars;

    float* x_e_box = run_globals.reion_grids.x_e_box;
    float* x_e_box_prev = run_globals.reion_grids.x_e_box_prev;
    float* Tk_box = run_globals.reion_grids.Tk_box;
    float* Tk_box_prev = run_globals.reion_grids.Tk_box_prev;
    float* TS_box = run_globals.reion_grids.TS_box;

    float* sfr = run_globals.reion_grids.sfr;
    fftwf_complex* sfr_unfiltered = (fftwf_complex*)sfr; // WATCH OUT!
    fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr, sfr_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    // TODO: Double check that looping over correct number of elements here
    for (int ii = 0; ii < slab_n_complex; ii++) {
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

	
	// Below I calculate the collapse fraction for all sources.
        // This should be zero (especially for the default high redshift Z_HEAT_MAX = 35). However, I compute it anyway in case it is non-zero.
        // In principle I think this should probably be used instead of Z_HEAT_MAX to switch between homogeneous/inhomogeneous.
        // However, I do not think it'll matter too much. Will look into this later.         
        
        collapse_fraction = 0.;

        R = ( L_FACTOR*box_size/(float)ReionGridDim ) / run_globals.params.Hubble_h;

        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                    density_over_mean = 1.0 + deltax[i_real];

                    collapse_fraction += stars[i_real] / (RtoM(R) * density_over_mean)
                                * (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;
                }

        MPI_Allreduce(MPI_IN_PLACE, &collapse_fraction, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        collapse_fraction = collapse_fraction/total_n_cells;

//        mlog("zp = %e collapse_fraction = %e", MLOG_MESG, collapse_fraction);
    }
    else {

        collapse_fraction = 0.;

        // Setup starting radius (minimum) and scaling to obtaining the maximum filtering radius for the X-ray background
        R = ( L_FACTOR*box_size/(float)ReionGridDim ) / run_globals.params.Hubble_h;
        R_factor = pow(R_XLy_MAX/R, 1/(float)NUM_FILTER_STEPS_FOR_Ts);

        // Smooth the density, stars and SFR fields over increasingly larger filtering radii (for evaluating the heating/ionisation integrals)
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){

            R_values[R_ct] = R;

            memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex);

            if(R_ct > 0) {
                int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

                filter(sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.HeatingFilterType);
            }
 
            // inverse fourier transform back to real space
            plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr_filtered, (float*)sfr_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);

            // Compute and store the collapse fraction and average electron fraction. Necessary for evaluating the integrals back along the light-cone.
            // Need the non-smoothed version, hence this is only done for R_ct == 0.
            if(R_ct == 0) {

                for (int ix = 0; ix < local_nix; ix++)
                    for (int iy = 0; iy < ReionGridDim; iy++)
                        for (int iz = 0; iz < ReionGridDim; iz++) {
                            i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                            i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                            ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

                            delNL0[R_ct][i_real] = ( ((float*)sfr_filtered)[i_padded] / pixel_volume )
                                     * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm, -3. ) * pow( run_globals.params.Hubble_h , -3. )/ SOLAR_MASS;

                            density_over_mean = 1.0 + deltax[i_real];

                            collapse_fraction += stars[i_real] / (RtoM(R) * density_over_mean)
                                * (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;

                            x_e_ave += x_e_box_prev[i_real];

                    }

                    MPI_Allreduce(MPI_IN_PLACE, &collapse_fraction, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
                    MPI_Allreduce(MPI_IN_PLACE, &x_e_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

                    collapse_fraction = collapse_fraction/total_n_cells;
                    x_e_ave = x_e_ave/total_n_cells;

                    stored_fcoll[snapshot] = collapse_fraction;

            }
            else {

                // Perform sanity checks to account for aliasing effects
                for (int ix = 0; ix < local_nix; ix++)
                    for (int iy = 0; iy < ReionGridDim; iy++)
                        for (int iz = 0; iz < ReionGridDim; iz++) {
                            i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

                            ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

                            delNL0[R_ct][i_real] = (((float*)sfr_filtered)[i_padded] / pixel_volume )
                                * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm, -3. ) * pow( run_globals.params.Hubble_h , -3. )/ SOLAR_MASS;

                        }

            }

            R *= R_factor;

        }
 
        // A condition (defined by whether or not there are stars) for evaluating the heating/ionisation integrals
        if(collapse_fraction > 0.0) {
            NO_LIGHT = 0;
        }
        else {
            NO_LIGHT = 1;
        }

        // Populate the initial ionisation/heating tables
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){

            if (R_ct==0){
                prev_zpp = zp;
                prev_R = 0;
            }
            else{
                prev_zpp = zpp_edge[R_ct-1];
                prev_R = R_values[R_ct-1];
            }

            zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*MPC / drdz(prev_zpp); // cell size
            zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''

            filling_factor_of_HI_zp = 1. - ReionEfficiency * collapse_fraction / (1.0 - x_e_ave);

            lower_int_limit_GAL = fmax(nu_tau_one(zp, zpp, x_e_ave, collapse_fraction, filling_factor_of_HI_zp, snapshot), run_globals.params.physics.NU_X_GAL_THRESH);

            if(run_globals.params.SEP_QSO_XRAY) {
                lower_int_limit_QSO = fmax(nu_tau_one(zp, zpp, x_e_ave, collapse_fraction, filling_factor_of_HI_zp, snapshot), run_globals.params.physics.NU_X_QSO_THRESH);
            }

            if (filling_factor_of_HI_zp < 0) filling_factor_of_HI_zp = 0; // for global evol; nu_tau_one above treats negative (post_reionization) inferred filling factors properly

            for (x_e_ct = 0; x_e_ct < x_int_NXHII; x_e_ct++){
                freq_int_heat_tbl_GAL[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit_GAL, run_globals.params.physics.NU_X_GAL_THRESH, run_globals.params.physics.X_RAY_SPEC_INDEX_GAL, 0);
                freq_int_ion_tbl_GAL[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit_GAL, run_globals.params.physics.NU_X_GAL_THRESH, run_globals.params.physics.X_RAY_SPEC_INDEX_GAL, 1);
                freq_int_lya_tbl_GAL[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit_GAL, run_globals.params.physics.NU_X_GAL_THRESH, run_globals.params.physics.X_RAY_SPEC_INDEX_GAL, 2);
            
                if(run_globals.params.SEP_QSO_XRAY)	{
                    freq_int_heat_tbl_QSO[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit_QSO, run_globals.params.physics.NU_X_QSO_THRESH, run_globals.params.physics.X_RAY_SPEC_INDEX_QSO, 0);
                    freq_int_ion_tbl_QSO[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit_QSO, run_globals.params.physics.NU_X_QSO_THRESH, run_globals.params.physics.X_RAY_SPEC_INDEX_QSO, 1);
                    freq_int_lya_tbl_QSO[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit_QSO, run_globals.params.physics.NU_X_QSO_THRESH, run_globals.params.physics.X_RAY_SPEC_INDEX_QSO, 2);
                }
            }            

            // and create the sum over Lya transitions from direct Lyn flux
            sum_lyn[R_ct] = 0;
            for (n_ct=NSPEC_MAX; n_ct>=2; n_ct--){
            if (zpp > zmax(zp, n_ct))
                continue;

                nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
                sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0);
            }
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
