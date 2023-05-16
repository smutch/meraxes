#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>

#include "XRayHeatingFunctions.h"
#include "meraxes.h"
#include "meraxes_gpu.h"
#include "misc_tools.h"
#include "recombinations.h"
#include "reionization.h"
#include "utils.h"

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were
 * subsequently made by Simon Mutch & Paul Geil.
 *
 * Inclusion of electron fraction (X-ray heating) and inhomogeneous recombinations
 * by Bradley Greig. Relevant functions taken from public version of 21cmFAST.
 */

double RtoM(double R)
{
  // All in internal units
  int filter = run_globals.params.ReionRtoMFilterType;
  double OmegaM = run_globals.params.OmegaM;
  double RhoCrit = run_globals.RhoCrit;

  switch (filter) {
    case 0: // top hat M = (4/3) PI <rho> R^3
      return (4.0 / 3.0) * M_PI * pow(R, 3) * (OmegaM * RhoCrit);
    case 1: // gaussian: M = (2PI)^1.5 <rho> R^3
      return pow(2 * M_PI, 1.5) * OmegaM * RhoCrit * pow(R, 3);
    default: // filter not defined
      mlog_error("Unrecognised filter (%d). Aborting...", filter);
      ABORT(EXIT_FAILURE);
      break;
  }

  return -1;
}

void _find_HII_bubbles(const int snapshot)
{
  // TODO: TAKE A VERY VERY CLOSE LOOK AT UNITS!!!!

  const double box_size = run_globals.params.BoxSize; // Mpc/h
  const int ReionGridDim = run_globals.params.ReionGridDim;
  const double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
  double cell_length_factor = L_FACTOR;
  const double total_n_cells = pow((double)ReionGridDim, 3);
  const int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
  const int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
  const int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
  const int flag_ReionUVBFlag = run_globals.params.ReionUVBFlag;
  const double ReionEfficiency = run_globals.params.physics.ReionEfficiency;
  const double ReionNionPhotPerBary = run_globals.params.physics.ReionNionPhotPerBary;
  run_units_t* units = &(run_globals.units);
  float J_21_aux = 0;
  double J_21_aux_constant;
  double density_over_mean;
  double weighted_sfr_density;
  double f_coll_stars;
  double electron_fraction;
  double Gamma_R_prefactor;

  double dNrec, rec;
  float z_eff;

  const double redshift = run_globals.ZZ[snapshot];
  double prev_redshift;
  if (snapshot == 0) {
    prev_redshift = run_globals.ZZ[snapshot];
  } else {
    prev_redshift = run_globals.ZZ[snapshot - 1];
  }

  float zstep = (float)(prev_redshift - redshift);
  float fabs_dtdz = (float)fabs(dtdz((float)redshift) / run_globals.params.Hubble_h);

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

  // Forward fourier transform to obtain k-space fields
  // TODO: Ensure that fftwf_mpi_init has been called and fftwf_mpi_cleanup will be called

  float* deltax = run_globals.reion_grids.deltax;
  fftwf_complex* deltax_unfiltered = run_globals.reion_grids.deltax_unfiltered;
  fftwf_complex* deltax_filtered = run_globals.reion_grids.deltax_filtered;
  fftwf_execute(run_globals.reion_grids.deltax_forward_plan);

  fftwf_complex* stars_unfiltered = run_globals.reion_grids.stars_unfiltered;
  fftwf_complex* stars_filtered = run_globals.reion_grids.stars_filtered;
  fftwf_execute(run_globals.reion_grids.stars_forward_plan);

  fftwf_complex* weighted_sfr_unfiltered = run_globals.reion_grids.weighted_sfr_unfiltered;
  fftwf_complex* weighted_sfr_filtered = run_globals.reion_grids.weighted_sfr_filtered;
  fftwf_execute(run_globals.reion_grids.weighted_sfr_forward_plan);

  // The free electron fraction from X-rays
  // TODO: Only necessary if we aren't using the GPU (not implemented there yet)
  float* x_e_box = NULL;
  fftwf_complex* x_e_unfiltered = NULL;
  fftwf_complex* x_e_filtered = NULL;
  if (run_globals.params.Flag_IncludeSpinTemp) {
    x_e_box = run_globals.reion_grids.x_e_box;
    x_e_unfiltered = run_globals.reion_grids.x_e_unfiltered;
    x_e_filtered = run_globals.reion_grids.x_e_filtered;
    fftwf_execute(run_globals.reion_grids.x_e_box_forward_plan);
  }

  // Fields relevant for computing the inhomogeneous recombinations
  float* Gamma12 = run_globals.reion_grids.Gamma12;
  float* N_rec = run_globals.reion_grids.N_rec;
  fftwf_complex* N_rec_unfiltered = NULL;
  fftwf_complex* N_rec_filtered = NULL;
  if (run_globals.params.Flag_IncludeRecombinations) {
    N_rec_unfiltered = run_globals.reion_grids.N_rec_unfiltered;
    N_rec_filtered = run_globals.reion_grids.N_rec_filtered;
    fftwf_execute(run_globals.reion_grids.N_rec_forward_plan);
  }

  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  // TODO: Double check that looping over correct number of elements here
  for (int ii = 0; ii < slab_n_complex; ii++) {
    deltax_unfiltered[ii] /= total_n_cells;
    stars_unfiltered[ii] /= total_n_cells;
    weighted_sfr_unfiltered[ii] /= total_n_cells;
    if (run_globals.params.Flag_IncludeRecombinations) {
      N_rec_unfiltered[ii] /= total_n_cells;
    }
    if (run_globals.params.Flag_IncludeSpinTemp) {
      x_e_unfiltered[ii] /= total_n_cells;
    }
  }

  // Loop through filter radii
  double ReionRBubbleMax;
  if (run_globals.params.Flag_IncludeRecombinations) {
    ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMaxRecomb; // Mpc/h
  } else {
    ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMax; // Mpc/h
  }
  double ReionRBubbleMin = run_globals.params.physics.ReionRBubbleMin; // Mpc/h
  double R = fmin(ReionRBubbleMax, L_FACTOR * box_size);               // Mpc/h
  double ReionDeltaRFactor = run_globals.params.ReionDeltaRFactor;
  double ReionGammaHaloBias = run_globals.params.physics.ReionGammaHaloBias;

  bool flag_last_filter_step = false;

  // set recombinations to zero (for case when recombinations are not used)
  rec = 0.0;

  while (!flag_last_filter_step) {
    // check to see if this is our last filtering step
    if (((R / ReionDeltaRFactor) <= (cell_length_factor * box_size / (double)ReionGridDim)) ||
        ((R / ReionDeltaRFactor) <= ReionRBubbleMin)) {
      flag_last_filter_step = true;
      R = cell_length_factor * box_size / (double)ReionGridDim;
    }

    // DEBUG
    // mlog("R = %.2e (h=0.678 -> %.2e)", MLOG_MESG, R, R/0.678);
    mlog(".", MLOG_CONT);

    // copy the k-space grids
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
    memcpy(stars_filtered, stars_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
    memcpy(weighted_sfr_filtered, weighted_sfr_unfiltered, sizeof(fftwf_complex) * slab_n_complex);

    if (run_globals.params.Flag_IncludeRecombinations) {
      memcpy(N_rec_filtered, N_rec_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
    }
    if (run_globals.params.Flag_IncludeSpinTemp) {
      memcpy(x_e_filtered, x_e_unfiltered, sizeof(fftwf_complex) * slab_n_complex);
    }

    // do the filtering unless this is the last filter step
    int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);
    if (!flag_last_filter_step) {
      filter(deltax_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
      filter(stars_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
      filter(
        weighted_sfr_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);

      if (run_globals.params.Flag_IncludeRecombinations) {
        filter(N_rec_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
      }
      if (run_globals.params.Flag_IncludeSpinTemp) {
        filter(x_e_filtered, local_ix_start, local_nix, ReionGridDim, (float)R, run_globals.params.ReionFilterType);
      }
    }

    // inverse fourier transform back to real space
    fftwf_execute(run_globals.reion_grids.deltax_filtered_reverse_plan);
    fftwf_execute(run_globals.reion_grids.stars_filtered_reverse_plan);
    fftwf_execute(run_globals.reion_grids.weighted_sfr_filtered_reverse_plan);

    if (run_globals.params.Flag_IncludeRecombinations) {
      fftwf_execute(run_globals.reion_grids.N_rec_filtered_reverse_plan);
    }

    if (run_globals.params.Flag_IncludeSpinTemp) {
      fftwf_execute(run_globals.reion_grids.x_e_filtered_reverse_plan);
    }

    // Perform sanity checks to account for aliasing effects
    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++) {
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
          ((float*)deltax_filtered)[i_padded] = fmaxf(((float*)deltax_filtered)[i_padded], -1 + REL_TOL);
          ((float*)stars_filtered)[i_padded] = fmaxf(((float*)stars_filtered)[i_padded], 0.0);
          if (((float*)stars_filtered)[i_padded] < ABS_TOL) {
            ((float*)stars_filtered)[i_padded] = 0;
          }
          ((float*)weighted_sfr_filtered)[i_padded] = fmaxf(((float*)weighted_sfr_filtered)[i_padded], 0.0);
          if (((float*)weighted_sfr_filtered)[i_padded] < ABS_TOL) {
            ((float*)weighted_sfr_filtered)[i_padded] = 0;
          }

          if (run_globals.params.Flag_IncludeRecombinations) {
            ((float*)N_rec_filtered)[i_padded] = fmaxf(((float*)N_rec_filtered)[i_padded], 0.0);
            if (((float*)N_rec_filtered)[i_padded] < ABS_TOL) {
              ((float*)N_rec_filtered)[i_padded] = 0;
            }
          }
          if (run_globals.params.Flag_IncludeSpinTemp) {
            ((float*)x_e_filtered)[i_padded] = fmaxf(((float*)x_e_filtered)[i_padded], 0.0);
            if (((float*)x_e_filtered)[i_padded] < ABS_TOL) {
              ((float*)x_e_filtered)[i_padded] = 0;
            }
            ((float*)x_e_filtered)[i_padded] = fminf(((float*)x_e_filtered)[i_padded], 0.999);
          }
        }

    /*
     * Main loop through the box...
     */

    J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI) * run_globals.params.physics.ReionAlphaUV *
                        PLANCK * 1e21 // * run_globals.params.physics.ReionEscapeFrac
                        * R * units->UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS * units->UnitMass_in_g /
                        pow(units->UnitLength_in_cm, 3) / units->UnitTime_in_s;

    if (run_globals.params.Flag_IncludeRecombinations) {
      Gamma_R_prefactor = (1.0 + redshift) * (1.0 + redshift) * R *
                          (units->UnitLength_in_cm / run_globals.params.Hubble_h) * SIGMA_HI *
                          run_globals.params.physics.ReionAlphaUV / (run_globals.params.physics.ReionAlphaUV + 2.75) /
                          1.0e-12; // Converting R h^-1 to R.
      Gamma_R_prefactor *= (units->UnitMass_in_g / units->UnitTime_in_s) *
                           pow(units->UnitLength_in_cm / run_globals.params.Hubble_h, -3.) * ReionNionPhotPerBary /
                           PROTONMASS; // Convert pixel volume (Mpc/h)^3 -> (cm)^3
    }

    double M_mean = RtoM(R);
    double R_cubed = R * R * R;

    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++) {
          i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
          i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);

          density_over_mean = 1.0 + (double)((float*)deltax_filtered)[i_padded];

          f_coll_stars = (double)((float*)stars_filtered)[i_padded] / (M_mean * density_over_mean) * (4.0 / 3.0) *
                         M_PI * R_cubed / pixel_volume;

          weighted_sfr_density = (double)((float*)weighted_sfr_filtered)[i_padded] / pixel_volume; // In internal units

          // Calculate the recombinations within the cell
          if (run_globals.params.Flag_IncludeRecombinations) {
            rec = (double)((float*)N_rec_filtered)[i_padded] / density_over_mean;
          }

          // Account for the partial ionisation of the cell from X-rays
          if (run_globals.params.Flag_IncludeSpinTemp) {
            electron_fraction = 1.0 - ((float*)x_e_filtered)[i_padded];
          } else {
            electron_fraction = 1.0;
          }

          if (flag_ReionUVBFlag)
            J_21_aux = (float)(weighted_sfr_density * J_21_aux_constant);

          // Modified reionisation condition, including recombinations and partial ionisations from X-rays
          // Check if ionised!
          if (f_coll_stars > (electron_fraction / ReionEfficiency) * (1. + rec)) // IONISED!!!!
          {
            // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
            if (xH[i_real] > REL_TOL) {
              if (flag_ReionUVBFlag)
                J_21[i_real] = J_21_aux;

              // Store the ionisation background and the reionisation redshift for each cell
              if (run_globals.params.Flag_IncludeRecombinations) {
                Gamma12[i_real] = (float)(Gamma_R_prefactor * weighted_sfr_density);
              }
            }

            // Mark as ionised
            xH[i_real] = 0;

            // Record radius
            r_bubble[i_real] = (float)R;
          }
          // Check if this is the last filtering step.
          // If so, assign partial ionisations to those cells which aren't fully ionised
          else if (flag_last_filter_step && (xH[i_real] > REL_TOL)) {
            xH[i_real] = (float)(electron_fraction - f_coll_stars * ReionEfficiency);
            if (xH[i_real] < 0.) {
              xH[i_real] = (float)0.;
            } else if (xH[i_real] > 1.0) {
              xH[i_real] = (float)1.;
            }
          }

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
  double volume_weighted_global_J_21 = 0.0;

  for (int ix = 0; ix < local_nix; ix++)
    for (int iy = 0; iy < ReionGridDim; iy++)
      for (int iz = 0; iz < ReionGridDim; iz++) {
        i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
        i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
        double cell_xH = (double)(xH[i_real]);
        volume_weighted_global_xH += cell_xH;

        if (flag_ReionUVBFlag) {
          volume_weighted_global_J_21 += (double)J_21[i_real];
        }

        density_over_mean = 1.0 + (double)((float*)deltax)[i_padded];
        mass_weighted_global_xH += cell_xH * density_over_mean;
        mass_weight += density_over_mean;

        if (run_globals.params.Flag_IncludeRecombinations) {
          // Store the resultant recombination cell
          z_eff = (float)((1. + redshift) * pow(density_over_mean, 1.0 / 3.0) - 1);
          dNrec = splined_recombination_rate(z_eff, Gamma12[i_real]) * fabs_dtdz * zstep * (1. - cell_xH);
          N_rec[i_padded] += dNrec;
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, &volume_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  if (flag_ReionUVBFlag) {
    MPI_Allreduce(MPI_IN_PLACE, &volume_weighted_global_J_21, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  }
  MPI_Allreduce(MPI_IN_PLACE, &mass_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &mass_weight, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

  volume_weighted_global_xH /= total_n_cells;
  if (flag_ReionUVBFlag) {
    volume_weighted_global_J_21 /= total_n_cells;
  }
  mass_weighted_global_xH /= mass_weight;
  run_globals.reion_grids.volume_weighted_global_xH = volume_weighted_global_xH;
  run_globals.reion_grids.volume_weighted_global_J_21 = volume_weighted_global_J_21;
  run_globals.reion_grids.mass_weighted_global_xH = mass_weighted_global_xH;
}

// This function makes sure that the right version of find_HII_bubbles() gets called.
void find_HII_bubbles(int snapshot, timer_info* timer_total)
{
  // Call the version of find_HII_bubbles we've been passed (and time it)
  double redshift = run_globals.ZZ[snapshot];
  timer_info timer;
#ifdef USE_CUDA
  mlog("Calling hybrid-GPU/FFTW version of find_HII_bubbles() for snap=%d/z=%.2lf...",
       MLOG_OPEN | MLOG_TIMERSTART,
       snapshot,
       redshift);

  // Run the GPU version of _find_HII_bubbles()
  timer_start(&timer);

  int flag_write_validation_data = false;
  _find_HII_bubbles_gpu(snapshot, flag_write_validation_data);
#else
  // Run the Meraxes version of _find_HII_bubbles()
  mlog("Calling pure-CPU version of find_HII_bubbles() for snap=%d/z=%.2lf...",
       MLOG_OPEN | MLOG_TIMERSTART,
       snapshot,
       redshift);
  timer_start(&timer);
  _find_HII_bubbles(snapshot);
#endif
  timer_stop(&timer);
  timer_stop(timer_total);
  timer_gpu += timer_delta(timer);
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
  mlog("Total time spent in find_HII_bubbles vs. total run time (snapshot %d ): %.2f of %.2f s",
       MLOG_MESG,
       snapshot,
       timer_gpu,
       timer_delta(*timer_total));
}
