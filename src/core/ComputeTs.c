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
 * 
 * Note: 21cmFAST has little_h included, therefore below I explicitly convert (using little_h) where appropriate. Be careful with units!
 */

void _ComputeTs(int snapshot)
{

    double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h ; // Mpc
    int ReionGridDim = run_globals.params.ReionGridDim;
    double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc)^3
    double total_n_cells = pow((double)ReionGridDim, 3);
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
    int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
    int slab_n_real_LC = local_nix * run_globals.params.NUM_FILTER_STEPS_FOR_Ts * ReionGridDim * ReionGridDim;
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

    int i_real, i_padded, i_smoothedSFR, R_ct, x_e_ct, n_ct, m_xHII_low, m_xHII_high, NO_LIGHT;

    double prev_zpp, prev_R, zpp, zp, lower_int_limit_GAL, lower_int_limit_QSO, filling_factor_of_HI_zp, R_factor, R, nuprime, dzp, Luminosity_converstion_factor_GAL, Luminosity_converstion_factor_QSO;
    double collapse_fraction, total_SFR, density_over_mean;

    float curr_xalpha;

    double freq_int_heat_GAL[run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_GAL[run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_GAL[run_globals.params.NUM_FILTER_STEPS_FOR_Ts];
    double freq_int_heat_QSO[run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_QSO[run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_QSO[run_globals.params.NUM_FILTER_STEPS_FOR_Ts];

    double freq_int_heat_tbl_GAL[x_int_NXHII][run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl_GAL[x_int_NXHII][run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl_GAL[x_int_NXHII][run_globals.params.NUM_FILTER_STEPS_FOR_Ts];
    double freq_int_heat_tbl_QSO[x_int_NXHII][run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl_QSO[x_int_NXHII][run_globals.params.NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl_QSO[x_int_NXHII][run_globals.params.NUM_FILTER_STEPS_FOR_Ts];

    double R_values[run_globals.params.NUM_FILTER_STEPS_FOR_Ts];

    double dt_dzpp_list[run_globals.params.NUM_FILTER_STEPS_FOR_Ts];

    double *evolve_ans, ans[2], dansdz[18], xHII_call;
    double SFR_GAL[run_globals.params.NUM_FILTER_STEPS_FOR_Ts], SFR_QSO[run_globals.params.NUM_FILTER_STEPS_FOR_Ts];

    float* deltax = run_globals.reion_grids.deltax;
    float* stars = run_globals.reion_grids.stars;

    float* x_e_box = run_globals.reion_grids.x_e_box;
    float* x_e_box_prev = run_globals.reion_grids.x_e_box_prev;
    float* Tk_box = run_globals.reion_grids.Tk_box;
    float* Tk_box_prev = run_globals.reion_grids.Tk_box_prev;
    float* TS_box = run_globals.reion_grids.TS_box;

    float* sfr = run_globals.reion_grids.sfr;
    float* sfr_temp = run_globals.reion_grids.sfr_temp;
    
    // Make a copy of the box for FFT'ing
    memcpy(sfr_temp, sfr, sizeof(fftwf_complex) * slab_n_complex);    

    fftwf_complex* sfr_unfiltered = (fftwf_complex*)sfr_temp; // WATCH OUT!
    fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr_temp, sfr_unfiltered, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    // TODO: Double check that looping over correct number of elements here
    for (int ii = 0; ii < slab_n_complex; ii++) {
        sfr_unfiltered[ii] /= total_n_cells;
    }

    int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

    double* SMOOTHED_SFR_GAL = run_globals.reion_grids.SMOOTHED_SFR_GAL;
    double* SMOOTHED_SFR_QSO; 
    if(run_globals.params.SEP_QSO_XRAY) {
        SMOOTHED_SFR_QSO = run_globals.reion_grids.SMOOTHED_SFR_QSO;
    }

    // Initialise the RECFAST, electron rate tables
    init_heat();

    x_e_ave = 0.0;

    
    double J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave;
    J_alpha_ave = xalpha_ave = Xheat_ave = Xion_ave = 0.0;

    double quantity1, quantity2, quantity3, quantity4, quantity5, quantity6, quantity7, quantity8, quantity9;
    double quantity10, quantity11, quantity12, quantity13, quantity14, quantity15, quantity16;
    quantity1 = quantity2 = quantity3 = quantity4 = quantity5 = quantity6 = quantity7 = quantity8 = quantity9 = 0.0;
    quantity10 = quantity11 = quantity12 = quantity13 = quantity14 = quantity15 = quantity16 = 0.0;

    // Place current redshift in 21cmFAST nomenclature (zp), delta zp (dzp) and delta z in seconds (dt_dzp)
    zp = redshift;
    dzp = zp - prev_redshift;
    dt_dzp = dtdz(zp);


    // Check redshift against Z_HEAT_MAX. If zp > Z_HEAT_MAX assume the x_e (electron fraction) and gas temperatures are homogenous 
    // Equivalent to the default setup of 21cmFAST. 
//    if( zp >= run_globals.params.physics.Z_HEAT_MAX && ( fabs(zp - run_globals.params.physics.Z_HEAT_MAX) >= FRACT_FLOAT_ERR) ) {
//    if( zp >= 34.0 ) {
    if( (zp - run_globals.params.physics.Z_HEAT_MAX) >= -0.0001) {
        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

                    x_e_box_prev[i_padded] = xion_RECFAST(zp,0);
                    Tk_box_prev[i_real] = T_RECFAST(zp,0);

                    TS_box[i_real] = get_Ts(zp, run_globals.reion_grids.deltax[i_padded], Tk_box_prev[i_real], x_e_box_prev[i_padded],0, &curr_xalpha);
          
                }

	
	// Below I calculate the collapse fraction for all sources.
        // This should be zero (especially for the default high redshift Z_HEAT_MAX = 35). However, I compute it anyway in case it is non-zero.
        // In principle I think this should probably be used instead of Z_HEAT_MAX to switch between homogeneous/inhomogeneous.
        // However, I do not think it'll matter too much. Will look into this later.         
        
        collapse_fraction = 0.;

        R = L_FACTOR*box_size/(float)ReionGridDim;  // Mpc

        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

                    density_over_mean = 1.0 + run_globals.reion_grids.deltax[i_padded];

                    // Multiplied by h^2 as RtoM(R) uses RhoCrit which doesn't include h factors
                    collapse_fraction += run_globals.reion_grids.stars[i_padded] / ( RtoM(R) * run_globals.params.Hubble_h * run_globals.params.Hubble_h * density_over_mean)
                                * (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;
                }

        MPI_Allreduce(MPI_IN_PLACE, &collapse_fraction, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        mlog("zp = %e collapse_fraction = %e",MLOG_MESG,zp,collapse_fraction);

        collapse_fraction = collapse_fraction/total_n_cells;

    }
    else {

        collapse_fraction = 0.;

        // Setup starting radius (minimum) and scaling to obtaining the maximum filtering radius for the X-ray background
        R = L_FACTOR*box_size/(float)ReionGridDim;
        R_factor = pow(R_XLy_MAX/R, 1/(float)run_globals.params.NUM_FILTER_STEPS_FOR_Ts);

        // Smooth the density, stars and SFR fields over increasingly larger filtering radii (for evaluating the heating/ionisation integrals)
        for (R_ct=0; R_ct<run_globals.params.NUM_FILTER_STEPS_FOR_Ts; R_ct++){

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

                quantity1 = 0.;

                for (int ix = 0; ix < local_nix; ix++)
                    for (int iy = 0; iy < ReionGridDim; iy++)
                        for (int iz = 0; iz < ReionGridDim; iz++) {
                            i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                            i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                            i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, run_globals.params.NUM_FILTER_STEPS_FOR_Ts, ReionGridDim);

                            ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

                            SMOOTHED_SFR_GAL[i_smoothedSFR] = ( ((float*)sfr_filtered)[i_padded] / pixel_volume )
                                     * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm, -3. )/ SOLAR_MASS;

                            quantity1 += SMOOTHED_SFR_GAL[i_smoothedSFR];

                            if(run_globals.params.SEP_QSO_XRAY) {
                                SMOOTHED_SFR_QSO[i_smoothedSFR] = ( ((float*)sfr_filtered)[i_padded] / pixel_volume )
                                     * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm, -3. )/ SOLAR_MASS;
                            }

                            density_over_mean = 1.0 + run_globals.reion_grids.deltax[i_padded];

                            collapse_fraction += run_globals.reion_grids.stars[i_padded] / (RtoM(R) * run_globals.params.Hubble_h * run_globals.params.Hubble_h * density_over_mean)
                                * (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;

                            x_e_ave += x_e_box_prev[i_padded];

                    }

                MPI_Allreduce(MPI_IN_PLACE, &collapse_fraction, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
                MPI_Allreduce(MPI_IN_PLACE, &x_e_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
                MPI_Allreduce(MPI_IN_PLACE, &quantity1, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

                collapse_fraction = collapse_fraction/total_n_cells;
                x_e_ave = x_e_ave/total_n_cells;
                quantity1 = quantity1/total_n_cells;

                stored_fcoll[snapshot] = collapse_fraction;

                mlog("zp = %e collapse_fraction = %e", MLOG_MESG, zp,collapse_fraction);
                mlog("zp = %e SFR_density = %e", MLOG_MESG, zp, quantity1);
 
            }
            else {

                // Perform sanity checks to account for aliasing effects
                for (int ix = 0; ix < local_nix; ix++)
                    for (int iy = 0; iy < ReionGridDim; iy++)
                        for (int iz = 0; iz < ReionGridDim; iz++) {
                            i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                            i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                            i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, run_globals.params.NUM_FILTER_STEPS_FOR_Ts, ReionGridDim);

                            ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

                            SMOOTHED_SFR_GAL[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume )
                                * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm, -3. )/ SOLAR_MASS;

                            if(run_globals.params.SEP_QSO_XRAY) {
                                SMOOTHED_SFR_QSO[i_smoothedSFR] = ( ((float*)sfr_filtered)[i_padded] / pixel_volume )
       	       	       	       	     * (units->UnitMass_in_g / units->UnitTime_in_s) * pow( units->UnitLength_in_cm, -3. )/ SOLAR_MASS;
       	       	       	    }

                        }

            }

            R *= R_factor;

        }

        mlog("zp = %e collapse_fraction = %e", MLOG_MESG, zp,collapse_fraction);
 
        // A condition (defined by whether or not there are stars) for evaluating the heating/ionisation integrals
        if(collapse_fraction > 0.0) {
            NO_LIGHT = 0;
        }
        else {
            NO_LIGHT = 1;
        }

        // Populate the initial ionisation/heating tables
        for (R_ct=0; R_ct<run_globals.params.NUM_FILTER_STEPS_FOR_Ts; R_ct++){

            if (R_ct==0){
                prev_zpp = zp;
                prev_R = 0;
            }
            else{
                prev_zpp = zpp_edge[R_ct-1];
                prev_R = R_values[R_ct-1];
            }

            zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*MPC / ( drdz(prev_zpp) ); // cell size
            zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''

            dt_dzpp_list[R_ct] = dtdz(zpp);

            filling_factor_of_HI_zp = 1. - ReionEfficiency * collapse_fraction / (1.0 - x_e_ave);

            lower_int_limit_GAL = fmax(nu_tau_one(zp, zpp, x_e_ave, collapse_fraction, filling_factor_of_HI_zp, snapshot), run_globals.params.physics.NU_X_GAL_THRESH*NU_over_EV);

            if(run_globals.params.SEP_QSO_XRAY) {
                lower_int_limit_QSO = fmax(nu_tau_one(zp, zpp, x_e_ave, collapse_fraction, filling_factor_of_HI_zp, snapshot), run_globals.params.physics.NU_X_QSO_THRESH*NU_over_EV);
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

        mlog("zp = %e collapse_fraction = %e", MLOG_MESG, zp,collapse_fraction);

        growth_factor_zp = dicke(zp);
        dgrowth_factor_dzp = ddicke_dz(zp);
        dt_dzp = dtdz(zp);
        
        // Below is the converstion of the soft-band X_ray luminosity into number of X-ray photons produced. This is the code taken from 21CMMC, which somewhat
        // uses the 21cmFAST nomenclature (to ease flipping between old/new parameterisation), so isn't necessarily the most intuitive way to express this.

        // Conversion of the input bolometric luminosity (new) to a ZETA_X (old) to be consistent with Ts.c from 21cmFAST
        // Conversion here means the code otherwise remains the same as the original Ts.c
        if(fabs(run_globals.params.physics.X_RAY_SPEC_INDEX_GAL - 1.0) < 0.000001) {
            Luminosity_converstion_factor_GAL = (run_globals.params.physics.NU_X_GAL_THRESH*NU_over_EV ) * log( run_globals.params.physics.NU_X_BAND_MAX/run_globals.params.physics.NU_X_GAL_THRESH );
            Luminosity_converstion_factor_GAL = 1./Luminosity_converstion_factor_GAL;
        }
        else {
            Luminosity_converstion_factor_GAL = pow( run_globals.params.physics.NU_X_BAND_MAX*NU_over_EV , 1. - run_globals.params.physics.X_RAY_SPEC_INDEX_GAL ) - pow( run_globals.params.physics.NU_X_GAL_THRESH*NU_over_EV , 1. - run_globals.params.physics.X_RAY_SPEC_INDEX_GAL ) ;
            Luminosity_converstion_factor_GAL = 1./Luminosity_converstion_factor_GAL;
            Luminosity_converstion_factor_GAL *= pow( run_globals.params.physics.NU_X_GAL_THRESH*NU_over_EV, - run_globals.params.physics.X_RAY_SPEC_INDEX_GAL )*(1 - run_globals.params.physics.X_RAY_SPEC_INDEX_GAL);
        }
	// Finally, convert to the correct units. NU_over_EV*hplank as only want to divide by eV -> erg (owing to the definition of Luminosity)
        Luminosity_converstion_factor_GAL *= (SEC_PER_YEAR)/(PLANCK);

        // Leave the original 21cmFAST code for reference. Refer to Greig & Mesinger (2017) for the new parameterisation.
//        const_zp_prefactor_GAL = (1.0/0.59)*( run_globals.params.physics.L_X_GAL * Luminosity_converstion_factor_GAL ) / (run_globals.params.physics.NU_X_GAL_THRESH*NU_over_EV) * C * pow(1+zp, run_globals.params.physics.X_RAY_SPEC_INDEX_GAL+3);
        const_zp_prefactor_GAL = ( run_globals.params.physics.L_X_GAL * Luminosity_converstion_factor_GAL ) / (run_globals.params.physics.NU_X_GAL_THRESH*NU_over_EV) * C * pow(1+zp, run_globals.params.physics.X_RAY_SPEC_INDEX_GAL+3);
        mlog("luminosity pre-factor = %e", MLOG_MESG, const_zp_prefactor_GAL);       
        // Note the factor of 0.59 appears to be required to match 21cmFAST

        // I believe it arises from differing definitions of a stellar baryon mass
       	// 21cmFAST appears to define a stellar baryon as 0.59*m_p, whereas Meraxes defines it as m_p
       	// Had issues with normalisation factors comparing the codes, adding 0.59 here rectified the normalisation somewhat
        // The Lya background appeared to be a factor of ~ 2 higher than 21cmFAST for the same luminosity. Spent days searching for the difference, this seems to explain it.
        // Can either boost the luminosity conversion, or lower the lya background by the same factor in XRayHeatingFunctions.c (evolveInt).
        // Note: When comparing to 21cmFAST it is important that this factor is included!
        // Will mean the normalisation within Meraxes is (1/0.59) higher than 21cmFAST, which can be trivially compensated for by reducing L_X.
        // Ultimately the backgrounds in Meraxes will be this same factor higher than 21cmFAST, but at least it is understood why and trivially accounted for.

        if(run_globals.params.SEP_QSO_XRAY) {
            
            if(fabs(run_globals.params.physics.X_RAY_SPEC_INDEX_QSO - 1.0) < 0.000001) {
                Luminosity_converstion_factor_QSO = (run_globals.params.physics.NU_X_QSO_THRESH*NU_over_EV) * log( run_globals.params.physics.NU_X_BAND_MAX/run_globals.params.physics.NU_X_QSO_THRESH );
                Luminosity_converstion_factor_QSO = 1./Luminosity_converstion_factor_QSO;
            }
            else {
                Luminosity_converstion_factor_QSO = pow( run_globals.params.physics.NU_X_BAND_MAX*NU_over_EV , 1. - run_globals.params.physics.X_RAY_SPEC_INDEX_QSO ) - pow( run_globals.params.physics.NU_X_QSO_THRESH*NU_over_EV , 1. - run_globals.params.physics.X_RAY_SPEC_INDEX_QSO ) ;
                Luminosity_converstion_factor_QSO = 1./Luminosity_converstion_factor_QSO;
                Luminosity_converstion_factor_QSO *= pow( run_globals.params.physics.NU_X_QSO_THRESH*NU_over_EV, - run_globals.params.physics.X_RAY_SPEC_INDEX_QSO )*(1 - run_globals.params.physics.X_RAY_SPEC_INDEX_QSO);
            }
            Luminosity_converstion_factor_QSO *= (SEC_PER_YEAR)/(PLANCK);

            // Leave the original 21cmFAST code for reference. Refer to Greig & Mesinger (2017) for the new parameterisation.
//            const_zp_prefactor_QSO = (1.0/0.59)*( run_globals.params.physics.L_X_QSO * Luminosity_converstion_factor_QSO ) / (run_globals.params.physics.NU_X_QSO_THRESH*NU_over_EV) * C * pow(1+zp, run_globals.params.physics.X_RAY_SPEC_INDEX_QSO+3);
            const_zp_prefactor_QSO = ( run_globals.params.physics.L_X_QSO * Luminosity_converstion_factor_QSO ) / (run_globals.params.physics.NU_X_QSO_THRESH*NU_over_EV) * C * pow(1+zp, run_globals.params.physics.X_RAY_SPEC_INDEX_QSO+3);
        }

        mlog("zp = %e collapse_fraction = %e", MLOG_MESG, zp,collapse_fraction);

        //interpolate to correct nu integral value based on the cell's ionization state
        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

                    ans[0] = x_e_box_prev[i_padded];
                    ans[1] = Tk_box_prev[i_real];

                    for (R_ct=0; R_ct<run_globals.params.NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                        i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, run_globals.params.NUM_FILTER_STEPS_FOR_Ts, ReionGridDim);

                        SFR_GAL[R_ct] = SMOOTHED_SFR_GAL[i_smoothedSFR];

                        if(run_globals.params.SEP_QSO_XRAY) {
                            SFR_QSO[R_ct] = SMOOTHED_SFR_QSO[i_smoothedSFR];
                        }

                        xHII_call = x_e_box_prev[i_padded];

                        dt_dzpp = dt_dzpp_list[R_ct];

                        // Check if ionized fraction is within boundaries; if not, adjust to be within
                        if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
                            xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
                        } else if (xHII_call < x_int_XHII[0]) {
                            xHII_call = 1.001*x_int_XHII[0];
                        }

                        m_xHII_low = locate_xHII_index(xHII_call);
                        m_xHII_high = m_xHII_low + 1;

                        // heat
                        freq_int_heat_GAL[R_ct] = (freq_int_heat_tbl_GAL[m_xHII_high][R_ct] - freq_int_heat_tbl_GAL[m_xHII_low][R_ct]) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                        freq_int_heat_GAL[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
                        freq_int_heat_GAL[R_ct] += freq_int_heat_tbl_GAL[m_xHII_low][R_ct];

                        // ionization
                        freq_int_ion_GAL[R_ct] = (freq_int_ion_tbl_GAL[m_xHII_high][R_ct] - freq_int_ion_tbl_GAL[m_xHII_low][R_ct]) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                        freq_int_ion_GAL[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
                        freq_int_ion_GAL[R_ct] += freq_int_ion_tbl_GAL[m_xHII_low][R_ct];

                        // lya
                        freq_int_lya_GAL[R_ct] = (freq_int_lya_tbl_GAL[m_xHII_high][R_ct] - freq_int_lya_tbl_GAL[m_xHII_low][R_ct]) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                        freq_int_lya_GAL[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
                        freq_int_lya_GAL[R_ct] += freq_int_lya_tbl_GAL[m_xHII_low][R_ct];

                        if(run_globals.params.SEP_QSO_XRAY) {

                            // heat
                            freq_int_heat_QSO[R_ct] = (freq_int_heat_tbl_QSO[m_xHII_high][R_ct] - freq_int_heat_tbl_QSO[m_xHII_low][R_ct]) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                            freq_int_heat_QSO[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
                            freq_int_heat_QSO[R_ct] += freq_int_heat_tbl_QSO[m_xHII_low][R_ct];

                            // ionization
                            freq_int_ion_QSO[R_ct] = (freq_int_ion_tbl_QSO[m_xHII_high][R_ct] - freq_int_ion_tbl_QSO[m_xHII_low][R_ct]) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                            freq_int_ion_QSO[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
                            freq_int_ion_QSO[R_ct] += freq_int_ion_tbl_QSO[m_xHII_low][R_ct];

                            // lya
                            freq_int_lya_QSO[R_ct] = (freq_int_lya_tbl_QSO[m_xHII_high][R_ct] - freq_int_lya_tbl_QSO[m_xHII_low][R_ct]) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                            freq_int_lya_QSO[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
                            freq_int_lya_QSO[R_ct] += freq_int_lya_tbl_QSO[m_xHII_low][R_ct];
                        }
                    }
 
                    // Perform the calculation of the heating/ionisation integrals, updating relevant quantities etc.
                    evolveInt(zp, run_globals.reion_grids.deltax[i_padded], SFR_GAL, SFR_QSO, freq_int_heat_GAL, freq_int_ion_GAL, freq_int_lya_GAL, freq_int_heat_QSO, freq_int_ion_QSO, freq_int_lya_QSO, NO_LIGHT, ans, dansdz);

                    quantity1 += dansdz[5];
       	       	    quantity2 += dansdz[6];
       	       	    quantity3 += dansdz[7];
       	       	    quantity4 += dansdz[8];

       	       	    quantity5 += dansdz[9];
       	       	    quantity6 += dansdz[10];
       	       	    quantity7 += dansdz[11];
       	       	    quantity8 += dansdz[12];
                    
                    quantity9 += dansdz[13];

                    quantity10 += dansdz[0];
                    quantity11 += dansdz[1];

                    quantity12 += dansdz[14];
                    quantity13 += dansdz[15];
                    quantity14 += dansdz[16];                    
                    quantity15 += dansdz[17];
                    quantity16 += dansdz[18];

                    x_e_box_prev[i_padded] += dansdz[0] * dzp; // remember dzp is negative
                    if (x_e_box_prev[i_padded] > 1) // can do this late in evolution if dzp is too large
                        x_e_box_prev[i_padded] = 1 - FRACT_FLOAT_ERR;
                    else if (x_e_box_prev[i_padded] < 0)
                        x_e_box_prev[i_padded] = 0;
                    if (Tk_box_prev[i_real] < MAX_TK)
                        Tk_box_prev[i_real] += dansdz[1] * dzp;

                    if (Tk_box_prev[i_real]<0){ // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
                        Tk_box_prev[i_real] = TCMB*(1+zp);
                    }

                    TS_box[i_real] = get_Ts(zp, run_globals.reion_grids.deltax[i_padded], Tk_box_prev[i_real], x_e_box_prev[i_padded], dansdz[2], &curr_xalpha);
                
                    J_alpha_ave += dansdz[2];
                    xalpha_ave += curr_xalpha;
                    Xheat_ave += dansdz[3];
                    Xion_ave += dansdz[4];

                } 

        MPI_Allreduce(MPI_IN_PLACE, &J_alpha_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &xalpha_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &Xheat_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &Xion_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        MPI_Allreduce(MPI_IN_PLACE, &quantity1, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &quantity2, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &quantity3, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &quantity4, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

       	MPI_Allreduce(MPI_IN_PLACE, &quantity5, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity6, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity7, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity8, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        MPI_Allreduce(MPI_IN_PLACE, &quantity9, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        MPI_Allreduce(MPI_IN_PLACE, &quantity10, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity11, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity12, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity13, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity14, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
       	MPI_Allreduce(MPI_IN_PLACE, &quantity15, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &quantity16, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        J_alpha_ave /= total_n_cells;
        xalpha_ave /= total_n_cells;
        Xheat_ave /= total_n_cells;
        Xion_ave /= total_n_cells;

        quantity1 /= total_n_cells;
       	quantity2 /= total_n_cells;
       	quantity3 /= total_n_cells;
       	quantity4 /= total_n_cells;

       	quantity5 /= total_n_cells;
        quantity6 /= total_n_cells;
        quantity7 /= total_n_cells;
        quantity8 /= total_n_cells;

        quantity9 /= total_n_cells;

        quantity10 /= total_n_cells;
        quantity11 /= total_n_cells;
        quantity12 /= total_n_cells;
        quantity13 /= total_n_cells;
        quantity14 /= total_n_cells;
        quantity15 /= total_n_cells;
        quantity16 /= total_n_cells;

        mlog("zp = %e collapse_fraction = %e", MLOG_MESG, zp,collapse_fraction);

        mlog("dxion_source_dt = %e dxheat_dt = %e dxlya_dt = %e dstarlya_dt = %e", MLOG_MESG, quantity1, quantity2, quantity3, quantity4);
        mlog("tables: heat = %e ion = %e = lya = %e star = %e", MLOG_MESG, quantity5, quantity6, quantity7, quantity8);
        mlog("SFR dz = %e",MLOG_MESG,quantity9);

        mlog("dxe_dzp = %e dxheat_dzp + dcomp_dzp + dspec_dzp + dadia_dzp = %e",MLOG_MESG,quantity10, quantity11);
        mlog("dxheat_dzp = %e dcomp_dzp = %e dspec_dzp = %e dadia_dzp = %e dxion_sink_dt = %e",MLOG_MESG,quantity16,quantity12, quantity13, quantity14, quantity15);

    }

    memcpy(x_e_box, x_e_box_prev, sizeof(fftwf_complex) * slab_n_complex);

    double Ave_Ts = 0.0;
    double Ave_x_e = 0.0;
    double Ave_Tk = 0.0;

    for (int ix = 0; ix < local_nix; ix++)
        for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
                i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                Ave_Ts += (double)TS_box[i_real];
                Ave_Tk += (double)Tk_box_prev[i_real];
                Ave_x_e += (double)x_e_box_prev[i_padded];
            }

    MPI_Allreduce(MPI_IN_PLACE, &Ave_Ts, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Ave_Tk, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Ave_x_e, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

    Ave_Ts /= total_n_cells;
    Ave_Tk /= total_n_cells;
    Ave_x_e /= total_n_cells;

    mlog("zp = %e collapse_fraction = %e ionising efficiency = %e product = %e", MLOG_MESG, zp, collapse_fraction, run_globals.params.physics.ReionEfficiency, run_globals.params.physics.ReionEfficiency * collapse_fraction);

    mlog("zp = %e Ts_ave = %e Tk_ave = %e x_e_ave = %e", MLOG_MESG, zp, Ave_Ts, Ave_Tk, Ave_x_e);
    mlog("zp = %e J_alpha_ave = %e xalpha_ave = %e Xheat_ave = %e Xion_ave = %e", MLOG_MESG, zp, J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave);
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
