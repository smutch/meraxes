#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

// DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>

#define L_FACTOR 0.620350491 // Factor relating cube length to filter radius = (4PI/3)^(-1/3)

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were
 * subsequently made by Simon Mutch & Paul Geil.
 */

double RtoM(double R)
{
    // All in internal units
    int filter = run_globals.params.ReionRtoMFilterType;
    double OmegaM = run_globals.params.OmegaM;
    double RhoCrit = run_globals.RhoCrit;

    switch (filter) {
    case 0: //top hat M = (4/3) PI <rho> R^3
        return (4.0 / 3.0) * M_PI * pow(R, 3) * (OmegaM * RhoCrit);
        break;
    case 1: //gaussian: M = (2PI)^1.5 <rho> R^3
        return pow(2 * M_PI, 1.5) * OmegaM * RhoCrit * pow(R, 3);
        break;
    default: // filter not defined
        mlog_error("Unrecognised filter (%d). Aborting...", filter);
        ABORT(EXIT_FAILURE);
        break;
    }

    return -1;
}

void find_HII_bubbles(double redshift)
{
    // TODO: TAKE A VERY VERY CLOSE LOOK AT UNITS!!!!

    double box_size = run_globals.params.BoxSize; // Mpc/h
    int ReionGridDim = run_globals.params.ReionGridDim;
    double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
    double cell_length_factor = L_FACTOR;
    double total_n_cells = pow((double)ReionGridDim, 3);
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
    int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
    int flag_ReionUVBFlag = run_globals.params.ReionUVBFlag;
    double ReionEfficiency = run_globals.params.physics.ReionEfficiency;
    double ReionNionPhotPerBary = run_globals.params.physics.ReionNionPhotPerBary;
    run_units_t* units = &(run_globals.units);
    float J_21_aux;
    double J_21_aux_constant;
    double density_over_mean;
    double sfr_density;
    double f_coll_stars;
    int i_real;
    int i_padded;

    // This parameter choice is sensitive to noise on the cell size, at least for the typical
    // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
    if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
        cell_length_factor = 1.0;

    // Init J_21
    float* J_21 = run_globals.reion_grids.J_21;
    if (flag_ReionUVBFlag)
        for (int ii = 0; ii < slab_n_real; ii++)
            J_21[ii] = 0.0;

    // Init xH
    float* xH = run_globals.reion_grids.xH;
    for (int ii = 0; ii < slab_n_real; ii++)
        xH[ii] = 1.0;

    // Init r_bubble
    float* r_bubble = run_globals.reion_grids.r_bubble;
    for (int ii = 0; ii < slab_n_real; ii++)
        r_bubble[ii] = 0.0;

    // #ifdef DEBUG
    //   {
    //     char fname_debug[STRLEN];
    //     sprintf(fname_debug, "pre_stars-%d.dat", run_globals.mpi_rank);
    //     FILE *fd_debug = fopen(fname_debug, "wb");
    //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
    //       for(int iy = 0; iy < ReionGridDim; iy++)
    //         for(int iz = 0; iz < ReionGridDim; iz++)
    //           fwrite(&(run_globals.reion_grids.stars[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
    //     fclose(fd_debug);
    //   }
    // #endif

    // #ifdef DEBUG
    //   {
    //     char fname_debug[STRLEN];
    //     sprintf(fname_debug, "pre_sfr-%d.dat", run_globals.mpi_rank);
    //     FILE *fd_debug = fopen(fname_debug, "wb");
    //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
    //       for(int iy = 0; iy < ReionGridDim; iy++)
    //         for(int iz = 0; iz < ReionGridDim; iz++)
    //           fwrite(&(run_globals.reion_grids.sfr[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
    //     fclose(fd_debug);
    //   }
    // #endif

    // Forward fourier transform to obtain k-space fields
    // TODO: Ensure that fftwf_mpi_init has been called and fftwf_mpi_cleanup will be called
    // TODO: Don't use estimate and calculate plan in code init
    float* deltax = run_globals.reion_grids.deltax;
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
    int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
    for (int ii = 0; ii < slab_n_complex; ii++) {
        deltax_unfiltered[ii] /= total_n_cells;
        stars_unfiltered[ii] /= total_n_cells;
        sfr_unfiltered[ii] /= total_n_cells;
    }

    // Loop through filter radii
    double ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMax; // Mpc/h
    double ReionRBubbleMin = run_globals.params.physics.ReionRBubbleMin; // Mpc/h
    double R = fmin(ReionRBubbleMax, L_FACTOR * box_size); // Mpc/h
    double ReionDeltaRFactor = run_globals.params.ReionDeltaRFactor;
    double ReionGammaHaloBias = run_globals.params.physics.ReionGammaHaloBias;

    bool flag_last_filter_step = false;

    while (!flag_last_filter_step) {
        // check to see if this is our last filtering step
        if (((R / ReionDeltaRFactor) <= (cell_length_factor * box_size / (double)ReionGridDim))
            || ((R / ReionDeltaRFactor) <= ReionRBubbleMin)) {
            flag_last_filter_step = true;
            R = cell_length_factor * box_size / (double)ReionGridDim;
        }

        // DEBUG
        // mlog("R = %.2e (h=0.678 -> %.2e)", MLOG_MESG, R, R/0.678);
        mlog(".", MLOG_CONT);

        // copy the k-space grids
        memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
        memcpy(stars_filtered, stars_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
        memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex);

        // do the filtering unless this is the last filter step
        int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);
        if (!flag_last_filter_step) {
            filter(deltax_filtered, local_ix_start, local_nix, ReionGridDim, (float)R);
            filter(stars_filtered, local_ix_start, local_nix, ReionGridDim, (float)R);
            filter(sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R);
        }

        // inverse fourier transform back to real space
        plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, deltax_filtered, (float*)deltax_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, stars_filtered, (float*)stars_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr_filtered, (float*)sfr_filtered, run_globals.mpi_comm, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        // Perform sanity checks to account for aliasing effects
        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                    ((float*)deltax_filtered)[i_padded] = fmaxf(((float*)deltax_filtered)[i_padded], -1 + REL_TOL);
                    ((float*)stars_filtered)[i_padded] = fmaxf(((float*)stars_filtered)[i_padded], 0.0);
                    ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);
                }

        // #ifdef DEBUG
        //   {
        //     char fname_debug[STRLEN];
        //     sprintf(fname_debug, "post_stars-%d_R%.3f.dat", run_globals.mpi_rank, R);
        //     FILE *fd_debug = fopen(fname_debug, "wb");
        //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
        //       for(int iy = 0; iy < ReionGridDim; iy++)
        //         for(int iz = 0; iz < ReionGridDim; iz++)
        //           fwrite(&(run_globals.reion_grids.stars_filtered[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
        //     fclose(fd_debug);
        //   }
        //   {
        //     char fname_debug[STRLEN];
        //     sprintf(fname_debug, "post_sfr-%d_R%.3f.dat", run_globals.mpi_rank, R);
        //     FILE *fd_debug = fopen(fname_debug, "wb");
        //     for(int ix = 0; ix < run_globals.reion_grids.slab_nix[run_globals.mpi_rank]; ix++)
        //       for(int iy = 0; iy < ReionGridDim; iy++)
        //         for(int iz = 0; iz < ReionGridDim; iz++)
        //           fwrite(&(run_globals.reion_grids.sfr_filtered[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)]), sizeof(float), 1, fd_debug);
        //     fclose(fd_debug);
        //   }
        // #endif

        /*
     * Main loop through the box...
     */

        J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI)
            * run_globals.params.physics.ReionAlphaUV * PLANCK
            * 1e21 * run_globals.params.physics.ReionEscapeFrac
            * R * units->UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS
            * units->UnitMass_in_g / pow(units->UnitLength_in_cm, 3) / units->UnitTime_in_s;

        // DEBUG
        // for (int ix=0; ix<local_nix; ix++)
        //   for (int iy=0; iy<ReionGridDim; iy++)
        //     for (int iz=0; iz<ReionGridDim; iz++)
        //       if (((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] > 0)
        //       {
        //         mlog("J_21_aux ==========", MLOG_OPEN);
        //         mlog("redshift = %.2f", MLOG_MESG, redshift);
        //         mlog("alpha = %.2e", MLOG_MESG, run_globals.params.physics.ReionAlphaUV);
        //         mlog("PLANCK = %.2e", MLOG_MESG, PLANCK);
        //         mlog("ReionEscapeFrac = %.2f", MLOG_MESG, ReionEscapeFrac);
        //         mlog("ReionNionPhotPerBary = %.1f", MLOG_MESG, ReionNionPhotPerBary);
        //         mlog("UnitMass_in_g = %.2e", MLOG_MESG, units->UnitMass_in_g);
        //         mlog("UnitLength_in_cm = %.2e", MLOG_MESG, units->UnitLength_in_cm);
        //         mlog("PROTONMASS = %.2e", MLOG_MESG, PROTONMASS);
        //         mlog("-> J_21_aux_constant = %.2e", MLOG_MESG, J_21_aux_constant);
        //         mlog("pixel_volume = %.2e", MLOG_MESG, pixel_volume);
        //         mlog("sfr_density[%d, %d, %d] = %.2e", MLOG_MESG, ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume);
        //         mlog("-> J_21_aux = %.2e", MLOG_MESG, ((float *)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] / pixel_volume * J_21_aux_constant);
        //         mlog("==========", MLOG_CLOSE);
        //         if (run_globals.mpi_rank == 0)
        //           ABORT(EXIT_SUCCESS);
        //       }

        for (int ix = 0; ix < local_nix; ix++)
            for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                    i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                    i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

                    density_over_mean = 1.0 + (double)((float*)deltax_filtered)[i_padded];

                    f_coll_stars = (double)((float*)stars_filtered)[i_padded] / (RtoM(R) * density_over_mean)
                        * (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;

                    sfr_density = (double)((float*)sfr_filtered)[i_padded] / pixel_volume; // In internal units

                    // #ifdef DEBUG
                    //           if(abs(redshift - 23.074) < 0.01)
                    //             debug("%d, %d, %d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
                    //                 run_globals.mpi_rank, ix, iy, iz,
                    //                 R, density_over_mean, f_coll_stars,
                    //                 ((float *)stars_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)],
                    //                 RtoM(R), (4.0/3.0)*M_PI*pow(R,3.0), pixel_volume, 1.0/ReionEfficiency, sfr_density, ((float*)sfr_filtered)[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)], (float)(sfr_density * J_21_aux_constant));
                    // #endif

                    if (flag_ReionUVBFlag)
                        J_21_aux = (float)(sfr_density * J_21_aux_constant);

                    // Check if ionised!
                    if (f_coll_stars > 1.0 / ReionEfficiency) // IONISED!!!!
                    {
                        // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
                        if (xH[i_real] > REL_TOL)
                            if (flag_ReionUVBFlag)
                                J_21[i_real] = J_21_aux;

                        // Mark as ionised
                        xH[i_real] = 0;

                        // Record radius
                        r_bubble[i_real] = (float)R;
                    }
                    // Check if this is the last filtering step.
                    // If so, assign partial ionisations to those cells which aren't fully ionised
                    else if (flag_last_filter_step && (xH[i_real] > REL_TOL))
                        xH[i_real] = (float)(1.0 - f_coll_stars * ReionEfficiency);

                    // Check if new ionisation
                    float* z_in = run_globals.reion_grids.z_at_ionization;
                    if ((xH[i_real] < REL_TOL) && (z_in[i_real] < 0)) // New ionisation!
                    {
                        z_in[i_real] = (float)redshift;
                        if (flag_ReionUVBFlag)
                            run_globals.reion_grids.J_21_at_ionization[i_real] = J_21_aux * (float)ReionGammaHaloBias;
                    }
                }
        // iz

        R /= ReionDeltaRFactor;
    }

    // Find the volume and mass weighted neutral fractions
    // TODO: The deltax grid will have rounding errors from forward and reverse
    //       FFT. Should cache deltax slabs prior to ffts and reuse here.
    double volume_weighted_global_xH = 0.0;
    double mass_weighted_global_xH = 0.0;
    double mass_weight = 0.0;

    for (int ix = 0; ix < local_nix; ix++)
        for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
                i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
                i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
                volume_weighted_global_xH += (double)xH[i_real];
                density_over_mean = 1.0 + (double)((float*)deltax_filtered)[i_padded];
                mass_weighted_global_xH += (double)(xH[i_real]) * density_over_mean;
                mass_weight += density_over_mean;
            }

    MPI_Allreduce(MPI_IN_PLACE, &volume_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &mass_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &mass_weight, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

    volume_weighted_global_xH /= total_n_cells;
    mass_weighted_global_xH /= mass_weight;
    run_globals.reion_grids.volume_weighted_global_xH = volume_weighted_global_xH;
    run_globals.reion_grids.mass_weighted_global_xH = mass_weighted_global_xH;
}
