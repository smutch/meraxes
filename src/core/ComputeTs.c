#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <string.h>

#include "ComputeTs.h"
#include "XRayHeatingFunctions.h"
#include "find_HII_bubbles.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "reionization.h"
#include "utils.h"

/*
 * This code is a re-write of the spin temperature calculation (Ts.c) within 21cmFAST.
 * Modified for usage within Meraxes by Bradley Greig.
 *
 * Note: 21cmFAST has little_h included, therefore below I explicitly convert (using little_h) where appropriate. Be
 * careful with units!
 */

/*
 * The minihalo feature was written by Emanuele M. Ventura, which includes an
 * amalgamation of the requisite functions for LW background using formulation of Qin2020a
 *
 * The output with II are computed accounting only for Pop II galaxies. In particular we compute Lyman-alpha, Xray, LW
 * only for Pop II. Notice that the reionization still accounts for Pop III galaxies so in this "only Pop II scenario"
 * we are still considering for the UV ionizing feedback from Pop III (xHI, and xe are computed using both PopII and Pop
 * III)
 */
void _ComputeTs(int snapshot)
{
  double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc
  int ReionGridDim = run_globals.params.ReionGridDim;
  double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc)^3
  double total_n_cells = pow((double)ReionGridDim, 3);
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
  int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
  double ReionEfficiency = run_globals.params.physics.ReionEfficiency;
  run_units_t* units = &(run_globals.units);
#if USE_MINI_HALOS
  double ReionEfficiencyIII = run_globals.params.physics.ReionEfficiencyIII;
#endif

  double redshift = run_globals.ZZ[snapshot];
  double prev_redshift;
  if (snapshot == 0) {
    prev_redshift = run_globals.ZZ[snapshot];
  } else {
    prev_redshift = run_globals.ZZ[snapshot - 1];
  }

  int i_real, i_padded, i_smoothedSFR, R_ct, x_e_ct, n_ct, m_xHII_low, m_xHII_high, NO_LIGHT;

  double prev_zpp, prev_R, zpp, zp, lower_int_limit_GAL, lower_int_limit_QSO, filling_factor_of_HI_zp, R_factor, R,
    nuprime, dzp, Luminosity_converstion_factor_GAL, Luminosity_converstion_factor_QSO;
  double collapse_fraction, density_over_mean, collapse_fraction_in_cell;

#if USE_MINI_HALOS
  double Luminosity_converstion_factor_III, collapse_fractionIII, collapse_fractionIII_in_cell;
#endif

  // TODO: Can we reduce the scope of these variables and, if not, improve the names?
  double weight = 0;
  bool first_radii = true;
  bool first_zero = true;
  int n_pts_radii = 1000;

  float curr_xalpha;
  int TsNumFilterSteps = run_globals.params.TsNumFilterSteps;

  double freq_int_heat_GAL[TsNumFilterSteps], freq_int_ion_GAL[TsNumFilterSteps], freq_int_lya_GAL[TsNumFilterSteps];
  double freq_int_heat_QSO[TsNumFilterSteps], freq_int_ion_QSO[TsNumFilterSteps], freq_int_lya_QSO[TsNumFilterSteps];

#if USE_MINI_HALOS
  double freq_int_heat_III[TsNumFilterSteps], freq_int_ion_III[TsNumFilterSteps], freq_int_lya_III[TsNumFilterSteps];
#endif

  double freq_int_heat_tbl_GAL[x_int_NXHII][TsNumFilterSteps], freq_int_ion_tbl_GAL[x_int_NXHII][TsNumFilterSteps],
    freq_int_lya_tbl_GAL[x_int_NXHII][TsNumFilterSteps];
  double freq_int_heat_tbl_QSO[x_int_NXHII][TsNumFilterSteps], freq_int_ion_tbl_QSO[x_int_NXHII][TsNumFilterSteps],
    freq_int_lya_tbl_QSO[x_int_NXHII][TsNumFilterSteps];

#if USE_MINI_HALOS
  double freq_int_heat_tbl_III[x_int_NXHII][TsNumFilterSteps], freq_int_ion_tbl_III[x_int_NXHII][TsNumFilterSteps],
    freq_int_lya_tbl_III[x_int_NXHII][TsNumFilterSteps];
#endif

  double R_values[TsNumFilterSteps];

  double dt_dzpp_list[TsNumFilterSteps];

  double ans[3], dansdz[20], xHII_call;
  double SFR_GAL[TsNumFilterSteps], SFR_QSO[TsNumFilterSteps];

#if USE_MINI_HALOS
  double SFR_III[TsNumFilterSteps];
#endif

  float* x_e_box = run_globals.reion_grids.x_e_box;
  float* x_e_box_prev = run_globals.reion_grids.x_e_box_prev;
  float* Tk_box = run_globals.reion_grids.Tk_box;
  float* TS_box = run_globals.reion_grids.TS_box;

#if USE_MINI_HALOS
  float* Tk_boxII = run_globals.reion_grids.Tk_boxII;
  float* TS_boxII = run_globals.reion_grids.TS_boxII;
  float* JLW_box = run_globals.reion_grids.JLW_box;
  float* JLW_boxII = run_globals.reion_grids.JLW_boxII;
#endif

  fftwf_complex* sfr_unfiltered = run_globals.reion_grids.sfr_unfiltered;
  fftwf_complex* sfr_filtered = run_globals.reion_grids.sfr_filtered;
  fftwf_execute(run_globals.reion_grids.sfr_forward_plan);

#if USE_MINI_HALOS
  fftwf_complex* sfrIII_unfiltered = run_globals.reion_grids.sfrIII_unfiltered;
  fftwf_complex* sfrIII_filtered = run_globals.reion_grids.sfrIII_filtered;
  fftwf_execute(run_globals.reion_grids.sfrIII_forward_plan);
#endif

  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  // TODO: Double check that looping over correct number of elements here
  for (int ii = 0; ii < slab_n_complex; ii++) {
    sfr_unfiltered[ii] /= (float)total_n_cells;
#if USE_MINI_HALOS
    sfrIII_unfiltered[ii] /= (float)total_n_cells;
#endif
  }

  double* SMOOTHED_SFR_GAL = run_globals.reion_grids.SMOOTHED_SFR_GAL;
  double* SMOOTHED_SFR_QSO;
#if USE_MINI_HALOS
  double* SMOOTHED_SFR_III = run_globals.reion_grids.SMOOTHED_SFR_III;
#endif
  if (run_globals.params.Flag_SeparateQSOXrays) {
    SMOOTHED_SFR_QSO = run_globals.reion_grids.SMOOTHED_SFR_QSO;
  }

  // Initialise the RECFAST, electron rate tables
  init_heat();

  x_e_ave = 0.0;

  double J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave, J_LW_ave;
  J_alpha_ave = xalpha_ave = Xheat_ave = Xion_ave = J_LW_ave = 0.0;

#if USE_MINI_HALOS
  double J_alpha_aveII, xalpha_aveII, Xheat_aveII, J_LW_aveII;
  J_alpha_aveII = xalpha_aveII = Xheat_aveII = J_LW_aveII = 0.0;
#endif

  // Place current redshift in 21cmFAST nomenclature (zp), delta zp (dzp) and delta z in seconds (dt_dzp)
  zp = redshift;
  dzp = zp - prev_redshift;
  dt_dzp = dtdz((float)(float)zp);

  // Check redshift against ReionMaxHeatingRedshift. If zp > ReionMaxHeatingRedshift assume the x_e (electron fraction)
  // and gas temperatures are homogenous Equivalent to the default setup of 21cmFAST.
  if ((zp - run_globals.params.physics.ReionMaxHeatingRedshift) >= -0.0001) {
    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++) {
          i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

          x_e_box_prev[i_padded] = (float)(float)xion_RECFAST((float)zp, 0);
          Tk_box[i_real] = (float)(float)T_RECFAST((float)zp, 0);

          TS_box[i_real] = get_Ts((float)zp,
                                  run_globals.reion_grids.deltax[i_padded],
                                  Tk_box[i_real],
                                  x_e_box_prev[i_padded],
                                  0,
                                  &curr_xalpha);
#if USE_MINI_HALOS
          Tk_boxII[i_real] = Tk_box[i_real];
          TS_boxII[i_real] = TS_box[i_real]; // This is true because Jalpha = 0 so curr_xalpha is the same.
#endif
        }

    // Below I calculate the collapse fraction for all sources.
    // This should be zero (especially for the default high redshift ReionMaxHeatingRedshift = 35). However, I compute
    // it anyway in case it is non-zero. In principle I think this should probably be used instead of
    // ReionMaxHeatingRedshift to switch between homogeneous/inhomogeneous. However, I do not think it'll matter too
    // much. Will look into this later.

    collapse_fraction = 0.;
#if USE_MINI_HALOS
    collapse_fractionIII = 0.;
#endif

    R = L_FACTOR * box_size / (float)ReionGridDim; // Mpc

    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++) {
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

          density_over_mean = 1.0 + run_globals.reion_grids.deltax[i_padded];

          // Multiplied by h^2 as RtoM(R) uses RhoCrit which doesn't include h factors, Geil2016
          collapse_fraction_in_cell =
            run_globals.reion_grids.stars[i_padded] /
            (RtoM(R) * run_globals.params.Hubble_h * run_globals.params.Hubble_h * density_over_mean) * (4.0 / 3.0) *
            M_PI * pow(R, 3.0) / pixel_volume;
#if USE_MINI_HALOS
          collapse_fractionIII_in_cell =
            run_globals.reion_grids.starsIII[i_padded] /
            (RtoM(R) * run_globals.params.Hubble_h * run_globals.params.Hubble_h * density_over_mean) * (4.0 / 3.0) *
            M_PI * pow(R, 3.0) / pixel_volume;
#endif

          if (collapse_fraction_in_cell > 1.0) {
            collapse_fraction_in_cell = 1.0;
          }
          collapse_fraction += collapse_fraction_in_cell;

#if USE_MINI_HALOS
          if (collapse_fractionIII_in_cell > 1.0) {
            collapse_fractionIII_in_cell = 1.0;
          }
          collapse_fractionIII += collapse_fractionIII_in_cell;
#endif
        }

    MPI_Allreduce(MPI_IN_PLACE, &collapse_fraction, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    collapse_fraction = collapse_fraction / total_n_cells;
#if USE_MINI_HALOS
    MPI_Allreduce(MPI_IN_PLACE, &collapse_fractionIII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    collapse_fractionIII = collapse_fractionIII / total_n_cells;
#endif

  } else {

    collapse_fraction = 0.;
#if USE_MINI_HALOS
    collapse_fractionIII = 0.;
#endif

    // Setup starting radius (minimum) and scaling to obtaining the maximum filtering radius for the X-ray background
    R = L_FACTOR * box_size / (float)ReionGridDim;
    R_factor = pow(R_XLy_MAX / R, 1 / (float)TsNumFilterSteps);

    // Smooth the density, stars and SFR fields over increasingly larger filtering radii (for evaluating the
    // heating/ionisation integrals)
    for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {

      R_values[R_ct] = R;

      memcpy(sfr_filtered, sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
#if USE_MINI_HALOS
      memcpy(sfrIII_filtered, sfrIII_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
#endif

      if (R_ct > 0) {
        int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

        filter(sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.TsHeatingFilterType);
#if USE_MINI_HALOS
        filter(
          sfrIII_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.TsHeatingFilterType);
#endif
      }

      // inverse fourier transform back to real space
      fftwf_execute(run_globals.reion_grids.sfr_filtered_reverse_plan);
#if USE_MINI_HALOS
      fftwf_execute(run_globals.reion_grids.sfrIII_filtered_reverse_plan);
#endif

      // Compute and store the collapse fraction and average electron fraction. Necessary for evaluating the integrals
      // back along the light-cone. Need the non-smoothed version, hence this is only done for R_ct == 0.
      if (R_ct == 0) {

        for (int ix = 0; ix < local_nix; ix++)
          for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
              i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
              i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
              i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

              ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);

              SMOOTHED_SFR_GAL[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) *
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *
                                                pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
#if USE_MINI_HALOS
              ((float*)sfrIII_filtered)[i_padded] = fmaxf(((float*)sfrIII_filtered)[i_padded], 0.0);

              SMOOTHED_SFR_III[i_smoothedSFR] = (((float*)sfrIII_filtered)[i_padded] / pixel_volume) *
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *
                                                pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
#endif

              if (run_globals.params.Flag_SeparateQSOXrays) {
                SMOOTHED_SFR_QSO[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) *
                                                  (units->UnitMass_in_g / units->UnitTime_in_s) *
                                                  pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
              }

              density_over_mean = 1.0 + run_globals.reion_grids.deltax[i_padded];

              collapse_fraction_in_cell =
                run_globals.reion_grids.stars[i_padded] /
                (RtoM(R) * run_globals.params.Hubble_h * run_globals.params.Hubble_h * density_over_mean) *
                (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;
#if USE_MINI_HALOS
              collapse_fractionIII_in_cell =
                run_globals.reion_grids.starsIII[i_padded] /
                (RtoM(R) * run_globals.params.Hubble_h * run_globals.params.Hubble_h * density_over_mean) *
                (4.0 / 3.0) * M_PI * pow(R, 3.0) / pixel_volume;
#endif

              if (collapse_fraction_in_cell > 1.0) {
                collapse_fraction_in_cell = 1.0;
              }
              collapse_fraction += collapse_fraction_in_cell;

#if USE_MINI_HALOS
              if (collapse_fractionIII_in_cell > 1.0) {
                collapse_fractionIII_in_cell = 1.0;
              }
              collapse_fractionIII += collapse_fractionIII_in_cell;
#endif
              x_e_ave += x_e_box_prev[i_padded];
            }

        MPI_Allreduce(MPI_IN_PLACE, &collapse_fraction, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &x_e_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

        collapse_fraction = collapse_fraction / total_n_cells;
        x_e_ave = x_e_ave / total_n_cells;

        stored_fcoll[snapshot] = collapse_fraction;

#if USE_MINI_HALOS
        MPI_Allreduce(MPI_IN_PLACE, &collapse_fractionIII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        collapse_fractionIII = collapse_fractionIII / total_n_cells;
        stored_fcollIII[snapshot] = collapse_fractionIII;
#endif

      } else {

        // Perform sanity checks to account for aliasing effects
        for (int ix = 0; ix < local_nix; ix++)
          for (int iy = 0; iy < ReionGridDim; iy++)
            for (int iz = 0; iz < ReionGridDim; iz++) {
              i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
              i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
              i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

              ((float*)sfr_filtered)[i_padded] = fmaxf(((float*)sfr_filtered)[i_padded], 0.0);
#if USE_MINI_HALOS
              ((float*)sfrIII_filtered)[i_padded] = fmaxf(((float*)sfrIII_filtered)[i_padded], 0.0);
#endif

              SMOOTHED_SFR_GAL[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) *
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *
                                                pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
#if USE_MINI_HALOS
              SMOOTHED_SFR_III[i_smoothedSFR] = (((float*)sfrIII_filtered)[i_padded] / pixel_volume) *
                                                (units->UnitMass_in_g / units->UnitTime_in_s) *
                                                pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
#endif
              if (run_globals.params.Flag_SeparateQSOXrays) {
                SMOOTHED_SFR_QSO[i_smoothedSFR] = (((float*)sfr_filtered)[i_padded] / pixel_volume) *
                                                  (units->UnitMass_in_g / units->UnitTime_in_s) *
                                                  pow(units->UnitLength_in_cm, -3.) / SOLAR_MASS;
              }
            }
      }

      R *= R_factor;
    }

    // A condition (defined by whether or not there are stars) for evaluating the heating/ionisation integrals
    // if (collapse_fraction > 0.0) {
#if USE_MINI_HALOS
    if ((collapse_fraction + collapse_fractionIII) > 0.0) {
#else
    if (collapse_fraction > 0.0) {
#endif
      NO_LIGHT = 0;
    } else {
      NO_LIGHT = 1;
    }

    // Populate the initial ionisation/heating tables
    for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {

      if (R_ct == 0) {
        prev_zpp = zp;
        prev_R = 0;
      } else {
        prev_zpp = zpp_edge[R_ct - 1];
        prev_R = R_values[R_ct - 1];
      }

      zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R) * MPC / (drdz((float)prev_zpp)); // cell size
      zpp = (zpp_edge[R_ct] + prev_zpp) * 0.5; // average redshift value of shell: z'' + 0.5 * dz''

      dt_dzpp_list[R_ct] = dtdz((float)zpp);

#if USE_MINI_HALOS
      filling_factor_of_HI_zp =
        1. - (ReionEfficiency * collapse_fraction + ReionEfficiencyIII * collapse_fractionIII) / (1.0 - x_e_ave);
#else
      filling_factor_of_HI_zp = 1. - ReionEfficiency * collapse_fraction / (1.0 - x_e_ave);
#endif

      lower_int_limit_GAL = fmax(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp, snapshot),
                                 run_globals.params.physics.NuXrayGalThreshold * NU_over_EV);

      if (run_globals.params.Flag_SeparateQSOXrays) {
        lower_int_limit_QSO = fmax(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp, snapshot),
                                   run_globals.params.physics.NuXrayQSOThreshold * NU_over_EV);
      }

      if (filling_factor_of_HI_zp < 0)
        filling_factor_of_HI_zp =
          0; // for global evol; nu_tau_one above treats negative (post_reionization) inferred filling factors properly

      for (x_e_ct = 0; x_e_ct < x_int_NXHII; x_e_ct++) {
        freq_int_heat_tbl_GAL[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                                x_int_XHII[x_e_ct],
                                                                lower_int_limit_GAL,
                                                                run_globals.params.physics.NuXrayGalThreshold,
                                                                run_globals.params.physics.SpecIndexXrayGal,
                                                                0);
        freq_int_ion_tbl_GAL[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                               x_int_XHII[x_e_ct],
                                                               lower_int_limit_GAL,
                                                               run_globals.params.physics.NuXrayGalThreshold,
                                                               run_globals.params.physics.SpecIndexXrayGal,
                                                               1);
        freq_int_lya_tbl_GAL[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                               x_int_XHII[x_e_ct],
                                                               lower_int_limit_GAL,
                                                               run_globals.params.physics.NuXrayGalThreshold,
                                                               run_globals.params.physics.SpecIndexXrayGal,
                                                               2);

#if USE_MINI_HALOS
        freq_int_heat_tbl_III[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                                x_int_XHII[x_e_ct],
                                                                lower_int_limit_GAL,
                                                                run_globals.params.physics.NuXrayGalThreshold,
                                                                run_globals.params.physics.SpecIndexXrayIII,
                                                                0);
        freq_int_ion_tbl_III[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                               x_int_XHII[x_e_ct],
                                                               lower_int_limit_GAL,
                                                               run_globals.params.physics.NuXrayGalThreshold,
                                                               run_globals.params.physics.SpecIndexXrayIII,
                                                               1);
        freq_int_lya_tbl_III[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                               x_int_XHII[x_e_ct],
                                                               lower_int_limit_GAL,
                                                               run_globals.params.physics.NuXrayGalThreshold,
                                                               run_globals.params.physics.SpecIndexXrayIII,
                                                               2);
#endif

        if (run_globals.params.Flag_SeparateQSOXrays) {
          freq_int_heat_tbl_QSO[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                                  x_int_XHII[x_e_ct],
                                                                  lower_int_limit_QSO,
                                                                  run_globals.params.physics.NuXrayQSOThreshold,
                                                                  run_globals.params.physics.SpecIndexXrayQSO,
                                                                  0);
          freq_int_ion_tbl_QSO[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                                 x_int_XHII[x_e_ct],
                                                                 lower_int_limit_QSO,
                                                                 run_globals.params.physics.NuXrayQSOThreshold,
                                                                 run_globals.params.physics.SpecIndexXrayQSO,
                                                                 1);
          freq_int_lya_tbl_QSO[x_e_ct][R_ct] = integrate_over_nu(zp,
                                                                 x_int_XHII[x_e_ct],
                                                                 lower_int_limit_QSO,
                                                                 run_globals.params.physics.NuXrayQSOThreshold,
                                                                 run_globals.params.physics.SpecIndexXrayQSO,
                                                                 2);
        }
      }

      // and create the sum over Lya transitions from direct Lyn flux
      sum_lyn[R_ct] = 0;
#if USE_MINI_HALOS
      sum_lyn_III[R_ct] = 0;
      if (run_globals.params.Flag_IncludeLymanWerner) {
        sum_lyn_LW[R_ct] = 0;
        sum_lyn_LW_III[R_ct] = 0;
      }
#endif

      for (n_ct = NSPEC_MAX; n_ct >= 2; n_ct--) {
        if (zpp > zmax((float)zp, n_ct))
          continue;

        nuprime = nu_n(n_ct) * (1 + zpp) / (1.0 + zp);
        sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0, 2);
#if USE_MINI_HALOS
        sum_lyn_III[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0, 3);
        if (run_globals.params.Flag_IncludeLymanWerner) {
          if (nuprime < NU_LW / NU_LL)
            nuprime = NU_LW / NU_LL;
          if (nuprime > nu_n(n_ct + 1))
            continue;
          sum_lyn_LW[R_ct] += spectral_emissivity(nuprime, 2, 2);
          sum_lyn_LW_III[R_ct] += spectral_emissivity(nuprime, 2, 3);
        }

#endif
      }

      // Find if we need to add a partial contribution to a radii to avoid kinks in the Lyman-alpha flux
      // As we look at discrete radii (light-cone redshift, zpp) we can have two radii where one has a
      // contribution and the next (larger) radii has no contribution. However, if the number of filtering
      // steps were infinitely large, we would have contributions between these two discrete radii
      // Thus, this aims to add a weighted contribution to the first radii where this occurs to smooth out
      // kinks in the average Lyman-alpha flux.

      // Note: We do not apply this correction to the LW background as it is unaffected by this. It is only
      // the Lyn contribution that experiences the kink. Applying this correction to LW introduces kinks
      // into the otherwise smooth quantity
      if (R_ct > 2 && sum_lyn[R_ct] == 0.0 && sum_lyn[R_ct - 1] > 0. && first_radii) {

        // The current zpp for which we are getting zero contribution
        double trial_zpp_max = (prev_zpp - (R_values[R_ct] - prev_R) * MPC / drdz((float)prev_zpp) + prev_zpp) * 0.5;
        // The zpp for the previous radius for which we had a non-zero contribution
        double trial_zpp_min =
          (zpp_edge[R_ct - 2] - (R_values[R_ct - 1] - R_values[R_ct - 2]) * MPC / drdz((float)zpp_edge[R_ct - 2]) +
           zpp_edge[R_ct - 2]) *
          0.5;

        // Split the previous radii and current radii into n_pts_radii smaller radii (redshift) to have fine control of
        // where it transitions from zero to non-zero This is a coarse approximation as it assumes that the linear
        // sampling is a good representation of the different volumes of the shells (from different radii).
        for (int ii = 0; ii < n_pts_radii; ii++) {
          double trial_zpp = trial_zpp_min + (trial_zpp_max - trial_zpp_min) * (float)ii / ((float)n_pts_radii - 1.);

          int counter = 0;
          for (n_ct = NSPEC_MAX; n_ct >= 2; n_ct--) {
            if (trial_zpp > zmax(zp, n_ct))
              continue;

            counter += 1;
          }
          if (counter == 0 && first_zero) {
            first_zero = false;
            weight = (float)ii / (float)n_pts_radii;
          }
        }

        // Now add a non-zero contribution to the previously zero contribution
        // The amount is the weight, multplied by the contribution from the previous radii
        sum_lyn[R_ct] = weight * sum_lyn[R_ct - 1];
#if USE_MINI_HALOS
        sum_lyn_III[R_ct] = weight * sum_lyn_III[R_ct - 1]; // I am not really sure about this line!
        if (run_globals.params.Flag_IncludeLymanWerner)
          sum_lyn_LW[R_ct] = weight * sum_lyn_LW[R_ct - 1];
#endif
        first_radii = false;
      }
    }

    growth_factor_zp = dicke(zp);
    dgrowth_factor_dzp = ddicke_dz(zp);
    dt_dzp = dtdz((float)zp);

    // Below is the converstion of the soft-band X_ray luminosity into number of X-ray photons produced. This is the
    // code taken from 21CMMC, which somewhat uses the 21cmFAST nomenclature (to ease flipping between old/new
    // parameterisation), so isn't necessarily the most intuitive way to express this.

    // Conversion of the input bolometric luminosity (new) to a ZETA_X (old) to be consistent with Ts.c from 21cmFAST
    // Conversion here means the code otherwise remains the same as the original Ts.c
    if (fabs(run_globals.params.physics.SpecIndexXrayGal - 1.0) < 0.000001) {
      Luminosity_converstion_factor_GAL =
        (run_globals.params.physics.NuXrayGalThreshold * NU_over_EV) *
        log(run_globals.params.physics.NuXraySoftCut / run_globals.params.physics.NuXrayGalThreshold);
      Luminosity_converstion_factor_GAL = 1. / Luminosity_converstion_factor_GAL;
    } else {
      Luminosity_converstion_factor_GAL =
        pow(run_globals.params.physics.NuXraySoftCut * NU_over_EV, 1. - run_globals.params.physics.SpecIndexXrayGal) -
        pow(run_globals.params.physics.NuXrayGalThreshold * NU_over_EV,
            1. - run_globals.params.physics.SpecIndexXrayGal);
      Luminosity_converstion_factor_GAL = 1. / Luminosity_converstion_factor_GAL;
      Luminosity_converstion_factor_GAL *=
        pow(run_globals.params.physics.NuXrayGalThreshold * NU_over_EV, -run_globals.params.physics.SpecIndexXrayGal) *
        (1 - run_globals.params.physics.SpecIndexXrayGal);
    }
    // Finally, convert to the correct units. NU_over_EV*hplank as only want to divide by eV -> erg (owing to the
    // definition of Luminosity)

    Luminosity_converstion_factor_GAL *= (SEC_PER_YEAR) / (PLANCK);

    // Do the same for Pop III (Once again, there are a few things you are not 100% sure)

#if USE_MINI_HALOS
    if (fabs(run_globals.params.physics.SpecIndexXrayIII - 1.0) < 0.000001) {
      Luminosity_converstion_factor_III =
        (run_globals.params.physics.NuXrayGalThreshold * NU_over_EV) *
        log(run_globals.params.physics.NuXraySoftCut / run_globals.params.physics.NuXrayGalThreshold);
      Luminosity_converstion_factor_III = 1. / Luminosity_converstion_factor_III;
    } else {
      Luminosity_converstion_factor_III =
        pow(run_globals.params.physics.NuXraySoftCut * NU_over_EV, 1. - run_globals.params.physics.SpecIndexXrayIII) -
        pow(run_globals.params.physics.NuXrayGalThreshold * NU_over_EV,
            1. - run_globals.params.physics.SpecIndexXrayIII);
      Luminosity_converstion_factor_III = 1. / Luminosity_converstion_factor_III;
      Luminosity_converstion_factor_III *=
        pow(run_globals.params.physics.NuXrayGalThreshold * NU_over_EV, -run_globals.params.physics.SpecIndexXrayIII) *
        (1 - run_globals.params.physics.SpecIndexXrayIII);
    }

    Luminosity_converstion_factor_III *= (SEC_PER_YEAR) / (PLANCK);
#endif

    // Leave the original 21cmFAST code for reference. Refer to Greig & Mesinger (2017) for the new parameterisation.
    //        const_zp_prefactor_GAL = (1.0/0.59)*( run_globals.params.physics.LXrayGal *
    //        Luminosity_converstion_factor_GAL ) / (run_globals.params.physics.NuXrayGalThreshold*NU_over_EV) *
    //        SPEED_OF_LIGHT * pow(1+zp, run_globals.params.physics.SpecIndexXrayGal+3);
    const_zp_prefactor_GAL = (run_globals.params.physics.LXrayGal * Luminosity_converstion_factor_GAL) /
                             (run_globals.params.physics.NuXrayGalThreshold * NU_over_EV) * SPEED_OF_LIGHT *
                             pow(1 + zp, run_globals.params.physics.SpecIndexXrayGal + 3);
#if USE_MINI_HALOS
    const_zp_prefactor_III = (run_globals.params.physics.LXrayGalIII * Luminosity_converstion_factor_III) /
                             (run_globals.params.physics.NuXrayGalThreshold * NU_over_EV) * SPEED_OF_LIGHT *
                             pow(1 + zp, run_globals.params.physics.SpecIndexXrayIII + 3);
#endif
    // Note the factor of 0.59 appears to be required to match 21cmFAST

    // I believe it arises from differing definitions of a stellar baryon mass
    // 21cmFAST appears to define a stellar baryon as 0.59*m_p, whereas Meraxes defines it as m_p
    // Had issues with normalisation factors comparing the codes, adding 0.59 here rectified the normalisation somewhat
    // The Lya background appeared to be a factor of ~ 2 higher than 21cmFAST for the same luminosity. Spent days
    // searching for the difference, this seems to explain it. Can either boost the luminosity conversion, or lower the
    // lya background by the same factor in XRayHeatingFunctions.c (evolveInt). Note: When comparing to 21cmFAST it is
    // important that this factor is included! Will mean the normalisation within Meraxes is (1/0.59) higher than
    // 21cmFAST, which can be trivially compensated for by reducing L_X. Ultimately the backgrounds in Meraxes will be
    // this same factor higher than 21cmFAST, but at least it is understood why and trivially accounted for.

    if (run_globals.params.Flag_SeparateQSOXrays) {

      if (fabs(run_globals.params.physics.SpecIndexXrayQSO - 1.0) < 0.000001) {
        Luminosity_converstion_factor_QSO =
          (run_globals.params.physics.NuXrayQSOThreshold * NU_over_EV) *
          log(run_globals.params.physics.NuXraySoftCut / run_globals.params.physics.NuXrayQSOThreshold);
        Luminosity_converstion_factor_QSO = 1. / Luminosity_converstion_factor_QSO;
      } else {
        Luminosity_converstion_factor_QSO =
          pow(run_globals.params.physics.NuXraySoftCut * NU_over_EV, 1. - run_globals.params.physics.SpecIndexXrayQSO) -
          pow(run_globals.params.physics.NuXrayQSOThreshold * NU_over_EV,
              1. - run_globals.params.physics.SpecIndexXrayQSO);
        Luminosity_converstion_factor_QSO = 1. / Luminosity_converstion_factor_QSO;
        Luminosity_converstion_factor_QSO *= pow(run_globals.params.physics.NuXrayQSOThreshold * NU_over_EV,
                                                 -run_globals.params.physics.SpecIndexXrayQSO) *
                                             (1 - run_globals.params.physics.SpecIndexXrayQSO);
      }
      Luminosity_converstion_factor_QSO *= (SEC_PER_YEAR) / (PLANCK);

      // Leave the original 21cmFAST code for reference. Refer to Greig & Mesinger (2017) for the new parameterisation.
      //            const_zp_prefactor_QSO = (1.0/0.59)*( run_globals.params.physics.LXrayQSO *
      //            Luminosity_converstion_factor_QSO ) / (run_globals.params.physics.NuXrayQSOThreshold*NU_over_EV) *
      //            SPEED_OF_LIGHT * pow(1+zp, run_globals.params.physics.SpecIndexXrayQSO+3);
      const_zp_prefactor_QSO = (run_globals.params.physics.LXrayQSO * Luminosity_converstion_factor_QSO) /
                               (run_globals.params.physics.NuXrayQSOThreshold * NU_over_EV) * SPEED_OF_LIGHT *
                               pow(1 + zp, run_globals.params.physics.SpecIndexXrayQSO + 3);
    }

    // interpolate to correct nu integral value based on the cell's ionization state
    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++) {
          i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

          ans[0] = x_e_box_prev[i_padded];
          ans[1] = Tk_box[i_real];
#if USE_MINI_HALOS
          ans[2] = Tk_boxII[i_real];
#endif

          for (R_ct = 0; R_ct < TsNumFilterSteps; R_ct++) {
            i_smoothedSFR = grid_index_smoothedSFR(R_ct, ix, iy, iz, TsNumFilterSteps, ReionGridDim);

            SFR_GAL[R_ct] = SMOOTHED_SFR_GAL[i_smoothedSFR];
#if USE_MINI_HALOS
            SFR_III[R_ct] = SMOOTHED_SFR_III[i_smoothedSFR];
#endif

            if (run_globals.params.Flag_SeparateQSOXrays) {
              SFR_QSO[R_ct] = SMOOTHED_SFR_QSO[i_smoothedSFR];
            }

            xHII_call = x_e_box_prev[i_padded];

            dt_dzpp = dt_dzpp_list[R_ct];

            // Check if ionized fraction is within boundaries; if not, adjust to be within
            if (xHII_call > x_int_XHII[x_int_NXHII - 1] * 0.999) {
              xHII_call = x_int_XHII[x_int_NXHII - 1] * 0.999;
            } else if (xHII_call < x_int_XHII[0]) {
              xHII_call = 1.001 * x_int_XHII[0];
            }

            m_xHII_low = locate_xHII_index((float)xHII_call);
            m_xHII_high = m_xHII_low + 1;

            // heat
            freq_int_heat_GAL[R_ct] =
              (freq_int_heat_tbl_GAL[m_xHII_high][R_ct] - freq_int_heat_tbl_GAL[m_xHII_low][R_ct]) /
              (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_heat_GAL[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_heat_GAL[R_ct] += freq_int_heat_tbl_GAL[m_xHII_low][R_ct];

            // ionization
            freq_int_ion_GAL[R_ct] =
              (freq_int_ion_tbl_GAL[m_xHII_high][R_ct] - freq_int_ion_tbl_GAL[m_xHII_low][R_ct]) /
              (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_ion_GAL[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_ion_GAL[R_ct] += freq_int_ion_tbl_GAL[m_xHII_low][R_ct];

            // lya
            freq_int_lya_GAL[R_ct] =
              (freq_int_lya_tbl_GAL[m_xHII_high][R_ct] - freq_int_lya_tbl_GAL[m_xHII_low][R_ct]) /
              (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_lya_GAL[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_lya_GAL[R_ct] += freq_int_lya_tbl_GAL[m_xHII_low][R_ct];

#if USE_MINI_HALOS
            freq_int_heat_III[R_ct] =
              (freq_int_heat_tbl_III[m_xHII_high][R_ct] - freq_int_heat_tbl_III[m_xHII_low][R_ct]) /
              (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_heat_III[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_heat_III[R_ct] += freq_int_heat_tbl_III[m_xHII_low][R_ct];

            // ionization
            freq_int_ion_III[R_ct] =
              (freq_int_ion_tbl_III[m_xHII_high][R_ct] - freq_int_ion_tbl_III[m_xHII_low][R_ct]) /
              (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_ion_III[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_ion_III[R_ct] += freq_int_ion_tbl_III[m_xHII_low][R_ct];

            // lya
            freq_int_lya_III[R_ct] =
              (freq_int_lya_tbl_III[m_xHII_high][R_ct] - freq_int_lya_tbl_III[m_xHII_low][R_ct]) /
              (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_lya_III[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_lya_III[R_ct] += freq_int_lya_tbl_III[m_xHII_low][R_ct];
#endif

            if (run_globals.params.Flag_SeparateQSOXrays) {

              // heat
              freq_int_heat_QSO[R_ct] =
                (freq_int_heat_tbl_QSO[m_xHII_high][R_ct] - freq_int_heat_tbl_QSO[m_xHII_low][R_ct]) /
                (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
              freq_int_heat_QSO[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
              freq_int_heat_QSO[R_ct] += freq_int_heat_tbl_QSO[m_xHII_low][R_ct];

              // ionization
              freq_int_ion_QSO[R_ct] =
                (freq_int_ion_tbl_QSO[m_xHII_high][R_ct] - freq_int_ion_tbl_QSO[m_xHII_low][R_ct]) /
                (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
              freq_int_ion_QSO[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
              freq_int_ion_QSO[R_ct] += freq_int_ion_tbl_QSO[m_xHII_low][R_ct];

              // lya
              freq_int_lya_QSO[R_ct] =
                (freq_int_lya_tbl_QSO[m_xHII_high][R_ct] - freq_int_lya_tbl_QSO[m_xHII_low][R_ct]) /
                (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
              freq_int_lya_QSO[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
              freq_int_lya_QSO[R_ct] += freq_int_lya_tbl_QSO[m_xHII_low][R_ct];
            }
          }

          // Perform the calculation of the heating/ionisation integrals, updating relevant quantities etc. GET BACK AT
          // THIS FOR USE_MINI_HALOS!!
#if USE_MINI_HALOS
          evolveInt((float)zp,
                    run_globals.reion_grids.deltax[i_padded],
                    SFR_GAL,
                    SFR_III,
                    SFR_QSO,
                    freq_int_heat_GAL,
                    freq_int_ion_GAL,
                    freq_int_lya_GAL,
                    freq_int_heat_III,
                    freq_int_ion_III,
                    freq_int_lya_III,
                    freq_int_heat_QSO,
                    freq_int_ion_QSO,
                    freq_int_lya_QSO,
                    NO_LIGHT,
                    ans,
                    dansdz);
#else
          evolveInt((float)zp,
                    run_globals.reion_grids.deltax[i_padded],
                    SFR_GAL,
                    0,
                    SFR_QSO,
                    freq_int_heat_GAL,
                    freq_int_ion_GAL,
                    freq_int_lya_GAL,
                    0,
                    0,
                    0,
                    freq_int_heat_QSO,
                    freq_int_ion_QSO,
                    freq_int_lya_QSO,
                    NO_LIGHT,
                    ans,
                    dansdz);
#endif

          x_e_box_prev[i_padded] += dansdz[0] * dzp; // remember dzp is negative
          if (x_e_box_prev[i_padded] > 1)            // can do this late in evolution if dzp is too large
            x_e_box_prev[i_padded] = (float)(1 - FRACT_FLOAT_ERR);
          else if (x_e_box_prev[i_padded] < 0)
            x_e_box_prev[i_padded] = 0;
          if (Tk_box[i_real] < MAX_TK)
            Tk_box[i_real] += dansdz[1] * dzp;

#if USE_MINI_HALOS
          if (Tk_boxII[i_real] < MAX_TK)
            Tk_boxII[i_real] += dansdz[6] * dzp;
          if (run_globals.params.Flag_IncludeLymanWerner) {
            JLW_box[i_real] = dansdz[5];
            JLW_boxII[i_real] = dansdz[10];
          }
#endif
          if (Tk_box[i_real] <
              0) { // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
            Tk_box[i_real] = (float)(TCMB * (1 + zp));
          }

#if USE_MINI_HALOS
          if (Tk_boxII[i_real] <
              0) { // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
            Tk_boxII[i_real] = (float)(TCMB * (1 + zp));
          }
#endif

          TS_box[i_real] = get_Ts((float)zp,
                                  run_globals.reion_grids.deltax[i_padded],
                                  Tk_box[i_real],
                                  x_e_box_prev[i_padded],
                                  (float)dansdz[2],
                                  &curr_xalpha);
#if USE_MINI_HALOS
          TS_boxII[i_real] = get_Ts((float)zp,
                                    run_globals.reion_grids.deltax[i_padded],
                                    Tk_boxII[i_real],
                                    x_e_box_prev[i_padded],
                                    (float)dansdz[7],
                                    &curr_xalpha); // It should be correct, probably I don't need a new curr_xalphaII
#endif
          J_alpha_ave += dansdz[2];
          xalpha_ave += curr_xalpha; // Double check this! It might be saving the one from PopIII!
          Xheat_ave += dansdz[3];
          Xion_ave += dansdz[4];
#if USE_MINI_HALOS
          J_alpha_aveII += dansdz[7];
          Xheat_aveII += dansdz[8];
          if (run_globals.params.Flag_IncludeLymanWerner) {
            J_LW_ave += dansdz[5];
            J_LW_aveII += dansdz[10];
          }
#endif
        }

    MPI_Allreduce(MPI_IN_PLACE, &J_alpha_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &xalpha_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Xheat_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Xion_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
#if USE_MINI_HALOS
    MPI_Allreduce(MPI_IN_PLACE, &J_alpha_aveII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &Xheat_aveII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &J_LW_ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &J_LW_aveII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
#endif

    J_alpha_ave /= total_n_cells;
    xalpha_ave /= total_n_cells;
    Xheat_ave /= total_n_cells;
    Xion_ave /= total_n_cells;
#if USE_MINI_HALOS
    J_alpha_aveII /= total_n_cells;
    Xheat_aveII /= total_n_cells;
    if (run_globals.params.Flag_IncludeLymanWerner) {
      J_LW_ave /= total_n_cells;
      J_LW_aveII /= total_n_cells;
    }
#endif

    run_globals.reion_grids.volume_ave_J_alpha = J_alpha_ave;
    run_globals.reion_grids.volume_ave_xalpha = xalpha_ave;
    run_globals.reion_grids.volume_ave_Xheat = Xheat_ave;
    run_globals.reion_grids.volume_ave_Xion = Xion_ave;
#if USE_MINI_HALOS
    run_globals.reion_grids.volume_ave_J_alphaII = J_alpha_aveII;
    run_globals.reion_grids.volume_ave_XheatII = Xheat_aveII;
    if (run_globals.params.Flag_IncludeLymanWerner) {
      run_globals.reion_grids.volume_ave_J_LW = J_LW_ave;
      run_globals.reion_grids.volume_ave_J_LWII = J_LW_aveII;
    }
#endif
  }

  memcpy(x_e_box, x_e_box_prev, sizeof(fftwf_complex) * slab_n_complex);

  double Ave_Ts = 0.0;
  double Ave_x_e = 0.0;
  double Ave_Tk = 0.0;
#if USE_MINI_HALOS
  double Ave_TkII = 0.0;
  double Ave_TsII = 0.0;
#endif

  for (int ix = 0; ix < local_nix; ix++)
    for (int iy = 0; iy < ReionGridDim; iy++)
      for (int iz = 0; iz < ReionGridDim; iz++) {
        i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

        Ave_Ts += (double)TS_box[i_real];
        Ave_Tk += (double)Tk_box[i_real];
        Ave_x_e += (double)x_e_box_prev[i_padded];
#if USE_MINI_HALOS
        Ave_TsII += (double)TS_boxII[i_real];
        Ave_TkII += (double)Tk_boxII[i_real];
#endif
      }

  MPI_Allreduce(MPI_IN_PLACE, &Ave_Ts, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &Ave_Tk, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &Ave_x_e, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

  Ave_Ts /= total_n_cells;
  Ave_Tk /= total_n_cells;
  Ave_x_e /= total_n_cells;

  run_globals.reion_grids.volume_ave_TS = Ave_Ts;
  run_globals.reion_grids.volume_ave_TK = Ave_Tk;
  run_globals.reion_grids.volume_ave_xe = Ave_x_e;

#if USE_MINI_HALOS
  MPI_Allreduce(MPI_IN_PLACE, &Ave_TsII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &Ave_TkII, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

  Ave_TsII /= total_n_cells;
  Ave_TkII /= total_n_cells;

  run_globals.reion_grids.volume_ave_TSII = Ave_TsII;
  run_globals.reion_grids.volume_ave_TKII = Ave_TkII;
#endif

  destruct_heat();

#if USE_MINI_HALOS
  mlog("zp = %e Ts_ave = %e Ts_ave (PopII) = %e Tk_ave = %e Tk_ave (PopII) = %e x_e_ave = %e",
       MLOG_MESG,
       zp,
       Ave_Ts,
       Ave_TsII,
       Ave_Tk,
       Ave_TkII,
       Ave_x_e);
  mlog("zp = %e J_alpha_ave = %e J_alpha_ave (PopII) = %e xalpha_ave = %e Xheat_ave = %e Xheat_ave (PopII) = %e "
       "Xion_ave = %e J_LW_ave = %e J_LW_ave (PopII) = %e",
       MLOG_MESG,
       zp,
       J_alpha_ave,
       J_alpha_aveII,
       xalpha_ave,
       Xheat_ave,
       Xheat_aveII,
       Xion_ave,
       J_LW_ave,
       J_LW_aveII);
#else
  mlog("zp = %e Ts_ave = %e Tk_ave = %e x_e_ave = %e", MLOG_MESG, zp, Ave_Ts, Ave_Tk, Ave_x_e);
  mlog("zp = %e J_alpha_ave = %e xalpha_ave = %e Xheat_ave = %e Xion_ave = %e",
       MLOG_MESG,
       zp,
       J_alpha_ave,
       xalpha_ave,
       Xheat_ave,
       Xion_ave);
#endif
}

// This function makes sure that the right version of ComputeTs() gets called.
// Note: Only the CPU version works for now
void ComputeTs(int snapshot, timer_info* timer_total)
{
  // Call the version of ComputeTs we've been passed (and time it)
  timer_info timer;
  // double redshift = run_globals.ZZ[snapshot];
  // Run the Meraxes version of _ComputeTs()
  mlog("Calling pure-CPU version of ComputeTs() for snap=%d/z=%.2lf...",
       MLOG_OPEN | MLOG_TIMERSTART,
       snapshot,
       run_globals.ZZ[snapshot]);
  timer_start(&timer);
  _ComputeTs(snapshot);
  timer_stop(&timer);
  timer_stop(timer_total);
  timer_gpu += timer_delta(timer);
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
  mlog("Total time spent in ComputeTs vs. total run time (snapshot %d ): %.2f of %.2f s",
       MLOG_MESG,
       snapshot,
       timer_gpu,
       timer_delta(*timer_total));
}