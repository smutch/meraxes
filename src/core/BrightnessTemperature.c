#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

#include "BrightnessTemperature.h"
#include "XRayHeatingFunctions.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "reionization.h"

/*
 * This code computes the 21cm brightness temperature, which is somewhat
 * of a re-write of the first part of delta_T.c from 21cmFAST. Also includes
 * a prescription for line-of-sight redshift-space distortions, taken from
 * 21CMMC (Greig & Mesinger 2018).
 * Written by Bradley Greig.
 */

void ComputeBrightnessTemperatureBox(int snapshot)
{

  int ii, jj, kk, i_real, i_padded, iii;

  double T_rad, pixel_Ts_factor, pixel_deltax, pixel_x_HI;

  float* xH = run_globals.reion_grids.xH;
  float* deltax = run_globals.reion_grids.deltax;

  float* delta_T = run_globals.reion_grids.delta_T;

  int ReionGridDim = run_globals.params.ReionGridDim;
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
  int slab_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
  double total_n_cells = pow((double)ReionGridDim, 3);

  double BoxSize = run_globals.params.BoxSize / run_globals.params.Hubble_h;

  // Set some redshift dependant values
  float redshift = (float)(float)run_globals.ZZ[snapshot];
  float Hubble_h = (float)(float)run_globals.params.Hubble_h;
  float OmegaM = (float)(float)run_globals.params.OmegaM;
  float OmegaB = (float)(OmegaM * run_globals.params.BaryonFrac);
  float const_factor = (float)(27.0 * (OmegaB * Hubble_h * Hubble_h / 0.023) *
                               sqrt((0.15 / OmegaM / Hubble_h / Hubble_h) * (1.0 + redshift) / 10.0));

  float H_z = (float)(float)hubble(redshift);

  double max_v_deriv, gradient_component, dvdx, subcell_width, d1_low, d2_low, d1_high, d2_high, subcell_displacement;
  double x_val1, x_val2, RSD_pos_new, RSD_pos_new_boundary_low, RSD_pos_new_boundary_high, cell_distance,
    fraction_outside, fraction_within;

  x_val1 = 0.;
  x_val2 = 1.;

  float* x_pos = calloc(N_RSD_STEPS, sizeof(float));
  float* x_pos_offset = calloc(N_RSD_STEPS, sizeof(float));
  float* delta_T_RSD_LOS = calloc((size_t)ReionGridDim, sizeof(float));

  // Initialise the 21cm brightness temperature box
  for (ii = 0; ii < local_nix; ii++) {
    for (jj = 0; jj < ReionGridDim; jj++)
      for (kk = 0; kk < ReionGridDim; kk++) {
        i_real = grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL);

        delta_T[i_real] = 0.0;
      }
  }

  T_rad = TCMB * (1. + redshift);

  for (ii = 0; ii < local_nix; ii++) {
    for (jj = 0; jj < ReionGridDim; jj++) {
      for (kk = 0; kk < ReionGridDim; kk++) {

        i_real = grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL);
        i_padded = grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED);

        pixel_deltax = deltax[i_padded];
        pixel_x_HI = xH[i_real];

        delta_T[i_real] = (float)(const_factor * pixel_x_HI * (1.0 + pixel_deltax));

        // Whether or not to include the spin temperature effects
        if (run_globals.params.Flag_IncludeSpinTemp) {

          // We are not using the optically thin limit approximation, compute full 21cm brightness temperature
          if (run_globals.params.Flag_IncludePecVelsFor21cm > 1) {

            // Converting the prefactors into the optical depth, tau. Factor of 1000 is the conversion of spin
            // temperature from K to mK
            delta_T[i_real] *= (1. + redshift) / (1000. * run_globals.reion_grids.TS_box[i_real]);

          } else {

            pixel_Ts_factor = (1 - T_rad / run_globals.reion_grids.TS_box[i_real]);
            delta_T[i_real] *= pixel_Ts_factor;
          }
        }
      }
    }
  }

  // Whether or not we are going to include the effects of peculiar velocities
  // There are two effects of peculiar velocities:
  // (i) modifications to the brightness temperature in real space owing to peculiar velocity gradients
  // (ii) observations are performed in redshift-space, and thus the signal is expanded/contracted by dense structures
  // (analagous to the fingers of god).

  // The input flag Flag_IncludePecVelsFor21cm controls how peculiar velocities will be treated within Meraxes.
  // Flag_IncludePecVelsFor21cm = 0: No velocity effects included
  // Flag_IncludePecVelsFor21cm = 1: Only (i) included, and done so by capping the maximum allowed peculiar velocity
  // gradient (this is the 21cmFAST treatment) Flag_IncludePecVelsFor21cm = 2: Only (i) but allowing for any peculiar
  // velocity gradient (i.e. uncapped) Flag_IncludePecVelsFor21cm = 3: Both (i) and (ii). Uncapped peculiar velocity
  // gradients and a treatment for line-of-sight peculiar velocity gradients from Greig & Mesinger (2018),
  //                                 based on a modified version of Jensen et al. (2013) & Mao et al. (2012)

  // **NOTE** : line-of-sight RSDs can be set in the saturated spin temperature limit (TS >> TCMB), however it defaults
  // to capping the maximum allowed peculiar velocity gradient.
  //          : This arises owing to the optical thin approximation implicit when TS >> TCMB.

  subcell_width = (BoxSize / ((double)ReionGridDim)) / (double)N_RSD_STEPS;

  float* vel;
  fftwf_complex* vel_gradient;

  double min_vel, max_vel, min_vel_grad, max_vel_grad;

  min_vel = 0.0;
  max_vel = 0.0;

  min_vel_grad = 0.0;
  max_vel_grad = 0.0;

  // If we aren't using velocities, then we are already done!
  if ((run_globals.params.Flag_IncludeSpinTemp) && (run_globals.params.Flag_IncludePecVelsFor21cm > 0)) {

    // Compute the velocity gradient, given the velocity field
    vel = run_globals.reion_grids.vel;

    // Temporary fix to the potential units issue with the velocity field
    // Multiply by sqrt(a) to convert Gadget internal units to proper velocities.
    // Dividing by 1000. because I believe the units are actually m/s not km/s (too large otherwise)
    // I am just going to divide by 1000. until I recieve confirmation
    for (ii = 0; ii < local_nix; ii++) {
      for (jj = 0; jj < ReionGridDim; jj++) {
        for (kk = 0; kk < ReionGridDim; kk++) {
          i_padded = grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED);

          vel[i_padded] = (float)(sqrt(1. / (1. + redshift)) * vel[i_padded] / 1000.);

          vel[i_padded] *= 1000. * 100. / MPC;

          // The algorithm used in 21cmFAST for dealing with velocities uses the **comoving** velocity. Therefore,
          // convert from peculiar to comoving velocity
          vel[i_padded] = (float)(vel[i_padded] / (1. / (1. + redshift)));

          if (vel[i_padded] < min_vel) {
            min_vel = vel[i_padded];
          }
          if (vel[i_padded] > max_vel) {
            max_vel = vel[i_padded];
          }
        }
      }
    }

    vel_gradient = run_globals.reion_grids.vel_gradient;
    fftwf_execute(run_globals.reion_grids.vel_forward_plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    // TODO: Double check that looping over correct number of elements here
    for (ii = 0; ii < slab_n_complex; ii++) {
      vel_gradient[ii] /= total_n_cells;

      // Include h_factor
      vel_gradient[ii] *= run_globals.params.Hubble_h;
    }

    int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);
    velocity_gradient(vel_gradient, local_ix_start, local_nix, ReionGridDim);

    fftwf_execute(run_globals.reion_grids.vel_gradient_reverse_plan);

    if (run_globals.params.Flag_IncludePecVelsFor21cm == 1) {

      // now add the velocity correction to the delta_T maps (only used for T_S >> T_CMB case).
      max_v_deriv = fabs(MAX_DVDR * H_z);

      for (ii = 0; ii < local_nix; ii++) {
        for (jj = 0; jj < ReionGridDim; jj++) {
          for (kk = 0; kk < ReionGridDim; kk++) {

            i_real = grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL);
            i_padded = grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED);

            dvdx = (double)((float*)vel_gradient)[i_padded];

            // set maximum allowed gradient for this linear approximation
            if (fabs(dvdx) > max_v_deriv) {
              if (dvdx < 0)
                dvdx = -max_v_deriv;
              else
                dvdx = max_v_deriv;
            }

            delta_T[i_real] /= (dvdx / H_z + 1.0);
          }
        }
      }
    }

    // now add the line-of-sight velocity correction to the 21cm brightness temperature maps (line-of-sight redshift
    // space distortions)
    if (run_globals.params.Flag_IncludePecVelsFor21cm > 1) {

      for (ii = 0; ii < local_nix; ii++) {
        for (jj = 0; jj < ReionGridDim; jj++) {
          for (kk = 0; kk < ReionGridDim; kk++) {

            i_real = grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL);
            i_padded = grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED);

            gradient_component = fabs((double)((float*)vel_gradient)[i_padded] / H_z + 1.0);

            if ((double)((float*)vel_gradient)[i_padded] < min_vel_grad) {
              min_vel_grad = (double)((float*)vel_gradient)[i_padded];
            }
            if ((double)((float*)vel_gradient)[i_padded] > max_vel_grad) {
              max_vel_grad = (double)((float*)vel_gradient)[i_padded];
            }

            if (run_globals.params.Flag_IncludeSpinTemp) {

              // Calculate the brightness temperature, using the optical depth
              if (gradient_component < FRACT_FLOAT_ERR) {

                // Gradient component goes to zero, optical depth diverges. But, since we take exp(-tau), this goes to
                // zero and (1 - exp(-tau)) goes to unity. Again, factors of 1000. are conversions from K to mK

                delta_T[i_real] = (float)(1000. * (run_globals.reion_grids.TS_box[i_real] - T_rad) / (1. + redshift));
              } else {
                delta_T[i_real] = (float)((1. - exp(-delta_T[i_real] / gradient_component)) * 1000. *
                                          (run_globals.reion_grids.TS_box[i_real] - T_rad) / (1. + redshift));
              }

            } else {

              dvdx = (double)((float*)vel_gradient)[i_padded];

              // set maximum allowed gradient for this linear approximation
              if (fabs(dvdx) > max_v_deriv) {
                if (dvdx < 0)
                  dvdx = -max_v_deriv;
                else
                  dvdx = max_v_deriv;
              }

              delta_T[i_real] /= (dvdx / H_z + 1.0);
            }
          }
        }
      }

      if (run_globals.params.Flag_IncludePecVelsFor21cm == 3) {

        // normalised units of cell length. 0 equals beginning of cell, 1 equals end of cell
        // These are the sub-cell central positions (x_pos_offset), and the corresponding normalised value (x_pos)
        // between 0 and 1
        for (iii = 0; iii < N_RSD_STEPS; iii++) {

          x_pos_offset[iii] = (float)(subcell_width * (float)iii + subcell_width / 2.);
          x_pos[iii] = (float)(x_pos_offset[iii] / (BoxSize / (float)ReionGridDim));
        }

        for (ii = 0; ii < local_nix; ii++) {
          for (jj = 0; jj < ReionGridDim; jj++) {

            for (kk = 0; kk < ReionGridDim; kk++) {
              delta_T_RSD_LOS[kk] = 0.0;
            }

            for (kk = 0; kk < ReionGridDim; kk++) {

              i_real = grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL);

              if ((fabs(delta_T[i_real]) >= FRACT_FLOAT_ERR) && (xH[i_real] >= FRACT_FLOAT_ERR)) {

                if (kk == 0) {
                  d1_low = vel[grid_index(ii, jj, ReionGridDim - 1, ReionGridDim, INDEX_PADDED)] / H_z;
                  d2_low = vel[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] / H_z;
                } else {
                  d1_low = vel[grid_index(ii, jj, kk - 1, ReionGridDim, INDEX_PADDED)] / H_z;
                  d2_low = vel[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] / H_z;
                }

                // Displacements (converted from velocity) for the original cell centres straddling half of the
                // sub-cells (cell after)
                if (kk == (ReionGridDim - 1)) {
                  d1_high = vel[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] / H_z;
                  d2_high = vel[grid_index(ii, jj, 0, ReionGridDim, INDEX_PADDED)] / H_z;
                } else {
                  d1_high = vel[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] / H_z;
                  d2_high = vel[grid_index(ii, jj, kk + 1, ReionGridDim, INDEX_PADDED)] / H_z;
                }

                for (iii = 0; iii < N_RSD_STEPS; iii++) {

                  // linearly interpolate the displacements to determine the corresponding displacements of the
                  // sub-cells Checking of 0.5 is for determining if we are left or right of the mid-point of the
                  // original cell (for the linear interpolation of the displacement) to use the appropriate cell

                  if (x_pos[iii] <= 0.5) {
                    subcell_displacement =
                      d1_low + ((x_pos[iii] + 0.5) - x_val1) * (d2_low - d1_low) / (x_val2 - x_val1);
                  } else {
                    subcell_displacement =
                      d1_high + ((x_pos[iii] - 0.5) - x_val1) * (d2_high - d1_high) / (x_val2 - x_val1);
                  }

                  // The new centre of the sub-cell post R.S.D displacement. Normalised to units of cell width for
                  // determining it's displacement
                  RSD_pos_new = (x_pos_offset[iii] + subcell_displacement) / (BoxSize / (float)ReionGridDim);

                  // The sub-cell boundaries of the sub-cell, for determining the fractional contribution of the
                  // sub-cell to neighbouring cells when the sub-cell straddles two cell positions
                  RSD_pos_new_boundary_low = RSD_pos_new - (subcell_width / 2.) / (BoxSize / (float)ReionGridDim);
                  RSD_pos_new_boundary_high = RSD_pos_new + (subcell_width / 2.) / (BoxSize / (float)ReionGridDim);

                  if (RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_high < 1.0) {
                    // sub-cell has remained in the original cell (just add it back to the original cell)

                    delta_T_RSD_LOS[kk] += delta_T[i_real] / (float)N_RSD_STEPS;
                  } else if (RSD_pos_new_boundary_low < 0.0 && RSD_pos_new_boundary_high < 0.0) {
                    // sub-cell has moved completely into a new cell (toward the observer)

                    // determine how far the sub-cell has moved in units of original cell boundary
                    cell_distance = ceil(fabs(RSD_pos_new_boundary_low)) - 1.;

                    // Determine the location of the sub-cell relative to the original cell binning
                    if (fabs(RSD_pos_new_boundary_high) > cell_distance) {
                      // sub-cell is entirely contained within the new cell (just add it to the new cell)

                      // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                      if (kk < ((int)cell_distance + 1)) {
                        delta_T_RSD_LOS[kk - ((int)cell_distance + 1) + ReionGridDim] +=
                          delta_T[i_real] / (float)N_RSD_STEPS;
                      } else {
                        delta_T_RSD_LOS[kk - ((int)cell_distance + 1)] += delta_T[i_real] / (float)N_RSD_STEPS;
                      }
                    } else {
                      // sub-cell is partially contained within the cell

                      // Determine the fraction of the sub-cell which is in either of the two original cells
                      fraction_outside = (fabs(RSD_pos_new_boundary_low) - cell_distance) /
                                         (subcell_width / (BoxSize / (float)ReionGridDim));
                      fraction_within = 1. - fraction_outside;

                      // Check if the first part of the sub-cell is at the box edge
                      if (kk < (((int)cell_distance))) {
                        delta_T_RSD_LOS[kk - ((int)cell_distance) + ReionGridDim] +=
                          fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                      } else {
                        delta_T_RSD_LOS[kk - ((int)cell_distance)] +=
                          fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                      }
                      // Check if the second part of the sub-cell is at the box edge
                      if (kk < (((int)cell_distance + 1))) {
                        delta_T_RSD_LOS[kk - ((int)cell_distance + 1) + ReionGridDim] +=
                          fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                      } else {
                        delta_T_RSD_LOS[kk - ((int)cell_distance + 1)] +=
                          fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                      }
                    }
                  } else if (RSD_pos_new_boundary_low < 0.0 &&
                             (RSD_pos_new_boundary_high > 0.0 && RSD_pos_new_boundary_high < 1.0)) {
                    // sub-cell has moved partially into a new cell (toward the observer)

                    // Determine the fraction of the sub-cell which is in either of the two original cells
                    fraction_within = RSD_pos_new_boundary_high / (subcell_width / (BoxSize / (float)ReionGridDim));
                    fraction_outside = 1. - fraction_within;

                    // Check the periodic boundaries conditions and move the fraction of each sub-cell to the
                    // appropriate new cell
                    if (kk == 0) {
                      delta_T_RSD_LOS[ReionGridDim - 1] += fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                      delta_T_RSD_LOS[kk] += fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                    } else {
                      delta_T_RSD_LOS[kk - 1] += fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                      delta_T_RSD_LOS[kk] += fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                    }
                  } else if ((RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_low < 1.0) &&
                             (RSD_pos_new_boundary_high >= 1.0)) {
                    // sub-cell has moved partially into a new cell (away from the observer)

                    // Determine the fraction of the sub-cell which is in either of the two original cells
                    fraction_outside =
                      (RSD_pos_new_boundary_high - 1.) / (subcell_width / (BoxSize / (float)ReionGridDim));
                    fraction_within = 1. - fraction_outside;

                    // Check the periodic boundaries conditions and move the fraction of each sub-cell to the
                    // appropriate new cell
                    if (kk == (ReionGridDim - 1)) {
                      delta_T_RSD_LOS[kk] += fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                      delta_T_RSD_LOS[0] += fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                    } else {
                      delta_T_RSD_LOS[kk] += fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                      delta_T_RSD_LOS[kk + 1] += fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                    }
                  } else {
                    // sub-cell has moved completely into a new cell (away from the observer)

                    // determine how far the sub-cell has moved in units of original cell boundary
                    cell_distance = floor(fabs(RSD_pos_new_boundary_high));

                    if (RSD_pos_new_boundary_low >= cell_distance) {
                      // sub-cell is entirely contained within the new cell (just add it to the new cell)

                      // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                      if (kk > (ReionGridDim - 1 - (int)cell_distance)) {
                        delta_T_RSD_LOS[kk + (int)cell_distance - ReionGridDim] += delta_T[i_real] / (float)N_RSD_STEPS;
                      } else {
                        delta_T_RSD_LOS[kk + (int)cell_distance] += delta_T[i_real] / (float)N_RSD_STEPS;
                      }
                    } else {
                      // sub-cell is partially contained within the cell

                      // Determine the fraction of the sub-cell which is in either of the two original cells
                      fraction_outside =
                        (RSD_pos_new_boundary_high - cell_distance) / (subcell_width / (BoxSize / (float)ReionGridDim));
                      fraction_within = 1. - fraction_outside;

                      // Check if the first part of the sub-cell is at the box edge
                      if (kk > (ReionGridDim - 1 - ((int)cell_distance - 1))) {
                        delta_T_RSD_LOS[kk + (int)cell_distance - 1 - ReionGridDim] +=
                          fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                      } else {
                        delta_T_RSD_LOS[kk + (int)cell_distance - 1] +=
                          fraction_within * delta_T[i_real] / (float)N_RSD_STEPS;
                      }
                      // Check if the second part of the sub-cell is at the box edge
                      if (kk > (ReionGridDim - 1 - ((int)cell_distance))) {
                        delta_T_RSD_LOS[kk + (int)cell_distance - ReionGridDim] +=
                          fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                      } else {
                        delta_T_RSD_LOS[kk + (int)cell_distance] +=
                          fraction_outside * delta_T[i_real] / (float)N_RSD_STEPS;
                      }
                    }
                  }
                }
              }
            }

            for (kk = 0; kk < ReionGridDim; kk++) {
              delta_T[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] = delta_T_RSD_LOS[kk];
            }
          }
        }
      }
      // End of line-of-sight redshift space distortions
    }
  }

  double Ave_Tb = 0.0;

  for (int ix = 0; ix < local_nix; ix++)
    for (int iy = 0; iy < ReionGridDim; iy++)
      for (int iz = 0; iz < ReionGridDim; iz++) {
        i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

        Ave_Tb += (double)delta_T[i_real];
      }

  Ave_Tb /= total_n_cells;

  MPI_Allreduce(MPI_IN_PLACE, &Ave_Tb, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
  mlog("zp = %e Tb_ave = %e", MLOG_MESG, redshift, Ave_Tb);

  free(delta_T_RSD_LOS);
  free(x_pos_offset);
  free(x_pos);

  run_globals.reion_grids.volume_ave_Tb = Ave_Tb;
}
