#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

#include "compute_power_spectrum.h"
#include "meraxes.h"
#include "misc_tools.h"

/*
 * A generic function to compute the 21cm PS of any field. Algorithm taken from delta_T.c
 * from 21cmFAST, but generalised to take any input cubic field.
 * Written by Bradley Greig
 *
 */

void Initialise_PowerSpectrum()
{

  double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc
  int ReionGridDim = run_globals.params.ReionGridDim;

  float k_factor = 1.35;
  float delta_k = (float)(2. * M_PI / box_size);
  float k_first_bin_ceil = delta_k;
  float k_max = delta_k * ReionGridDim;

  // Initialise arrays
  // ghetto counting (lookup how to do logs of arbitrary bases in c...)
  int num_bins = 0;
  float k_ceil = k_first_bin_ceil;
  float k_floor = k_ceil;
  while (k_ceil < k_max) {
    num_bins++;
    k_floor = k_ceil;
    k_ceil *= k_factor;
  }

  run_globals.params.PS_Length = num_bins;
  mlog("Initialise_PowerSpectrum set PS_Length to %d.", MLOG_MESG, run_globals.params.PS_Length);
}

void Compute_PS(int snapshot)
{

  float* delta_T = run_globals.reion_grids.delta_T;

  double box_size = run_globals.params.BoxSize / run_globals.params.Hubble_h; // Mpc

  fftwf_complex* deldel_ps = fftwf_alloc_complex((size_t)run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);

  float volume = powf((float)(float)box_size, 3);

  int ReionGridDim = run_globals.params.ReionGridDim;
  double total_n_cells = pow((double)ReionGridDim, 3);
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);

  int ii, jj, kk, n_x, n_y, n_z;

  double ave;

  mlog("Calculating the 21cm power spectrum (dimensional, i.e mK^2)", MLOG_MESG);

  ave = 0.0;
  for (ii = 0; ii < local_nix; ii++) {
    for (jj = 0; jj < ReionGridDim; jj++) {
      for (kk = 0; kk < ReionGridDim; kk++) {
        ave += delta_T[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)];
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &ave, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);

  ave /= total_n_cells;

  for (ii = 0; ii < local_nix; ii++) {
    for (jj = 0; jj < ReionGridDim; jj++) {
      for (kk = 0; kk < ReionGridDim; kk++) {
        ((float*)deldel_ps)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] =
          (float)((delta_T[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] / ave - 1) * volume /
                  (float)total_n_cells);
        ((float*)deldel_ps)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] *= ave;
      }
    }
  }

  fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(
    ReionGridDim, ReionGridDim, ReionGridDim, (float*)deldel_ps, deldel_ps, run_globals.mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  // Calculate power spectrum
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  unsigned long long ct;

  float k_mag;

  float k_factor = 1.35;
  float delta_k = (float)(2. * M_PI / box_size);
  float k_first_bin_ceil = delta_k;
  float k_max = delta_k * ReionGridDim;

  float k_floor = 0.0;
  float k_ceil = k_first_bin_ceil;

  double* p_box = malloc(sizeof(double) * run_globals.params.PS_Length);
  double* k_ave = malloc(sizeof(double) * run_globals.params.PS_Length);
  unsigned long long* in_bin_ct = malloc(sizeof(unsigned long long) * run_globals.params.PS_Length);

  for (ii = 0; ii < run_globals.params.PS_Length; ii++) {
    p_box[ii] = 0.0;
    k_ave[ii] = 0.0;
    in_bin_ct[ii] = 0;
  }

  // Co-eval box, so should sample the entire cube

  // now construct the power spectrum file
  int HII_middle = ReionGridDim / 2;

  int local_ix_start = (int)(run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank]);

  for (n_x = 0; n_x < local_nix; n_x++) {
    float k_x = (n_x + local_ix_start) * delta_k;
    if ((n_x + local_ix_start) > HII_middle) {
      k_x = ((n_x + local_ix_start) - ReionGridDim) * delta_k; // wrap around for FFT convention
    }

    for (n_y = 0; n_y < ReionGridDim; n_y++) {
      float k_y = n_y * delta_k;
      if (n_y > HII_middle)
        k_y = (n_y - ReionGridDim) * delta_k;

      for (n_z = 0; n_z <= HII_middle; n_z++) {
        float k_z = n_z * delta_k;

        k_mag = (float)sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

        // now go through the k bins and update
        ct = 0;
        k_floor = 0;
        k_ceil = k_first_bin_ceil;

        while (k_ceil < k_max) {
          // check if we fall in this bin
          if ((k_mag >= k_floor) && (k_mag < k_ceil)) {
            in_bin_ct[ct]++;
            p_box[ct] += pow(k_mag, 3) *
                         pow(cabs(deldel_ps[grid_index(n_x, n_y, n_z, ReionGridDim, INDEX_COMPLEX_HERM)]), 2.) /
                         (2.0 * M_PI * M_PI * volume);
            // note the 1/VOLUME factor, which turns this into a power density in k-space

            k_ave[ct] += k_mag;
            break;
          }

          ct++;
          k_floor = k_ceil;
          k_ceil *= k_factor;
        }
      }
    }
  } // end looping through k box

  float* PS_k = run_globals.reion_grids.PS_k;
  float* PS_data = run_globals.reion_grids.PS_data;
  float* PS_error = run_globals.reion_grids.PS_error;

  // NOTE - previous ct ran from 1 (not zero) to NUM_BINS
  for (ii = 0; ii < run_globals.params.PS_Length; ii++) {

    MPI_Allreduce(MPI_IN_PLACE, &k_ave[ii], 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &p_box[ii], 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &in_bin_ct[ii], 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, run_globals.mpi_comm);

    PS_k[ii] = (float)(k_ave[ii] / (double)in_bin_ct[ii]);
    PS_data[ii] = (float)(p_box[ii] / (double)in_bin_ct[ii]);
    PS_error[ii] = (float)((p_box[ii] / (double)in_bin_ct[ii]) / sqrt((double)in_bin_ct[ii]));
  }

  // Deallocate
  // ------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------

  free(p_box);
  free(k_ave);
  free(in_bin_ct);
  fftwf_free(deldel_ps);
}
