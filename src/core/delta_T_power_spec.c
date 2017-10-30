#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were
 * subsequently made by Simon Mutch & Paul Geil.
 */

int delta_T_ps(
    int snapshot,
    float* average_T,
    float* delta_T,
    float** ps,
    int* ps_nbins)
{
    unsigned long long ct, temp_ct;

    float pixel_x_HI, pixel_deltax;
    float k_mag;

    float* xH = run_globals.reion_grids.xH;
    float* deltax = run_globals.reion_grids.deltax;

    // Begin initialisation
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------

    float min = 1e3;
    float max = -1e3;
    double ave = 0.0;
    int ReionGridDim = run_globals.params.ReionGridDim;
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);

    // Set some redshift dependant values
    float redshift = run_globals.ZZ[snapshot];
    float Hubble_h = run_globals.params.Hubble_h;
    float OmegaM = run_globals.params.OmegaM;
    float OmegaB = OmegaM * run_globals.params.BaryonFrac;
    float const_factor = 27.0 * (OmegaB * Hubble_h * Hubble_h / 0.023) * sqrt((0.15 / OmegaM / Hubble_h / Hubble_h) * (1.0 + redshift) / 10.0);

    // delta_T grid (same size as the bubble box)
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------

    temp_ct = 0;

    for (int ii = 0; ii < local_nix; ii++) {
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++) {
                pixel_deltax = deltax[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];

                int index = grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL);
                pixel_x_HI = xH[index];

                if (pixel_x_HI > ABS_TOL)
                    temp_ct++;

                delta_T[index] = const_factor * pixel_x_HI * (1.0 + pixel_deltax);

                if (max < delta_T[index])
                    max = delta_T[index];

                if (min > delta_T[index])
                    min = delta_T[index];

                ave += delta_T[index];
            }
    }

    int tot_num_pixels = (int)pow(ReionGridDim, 3);

    MPI_Allreduce(MPI_IN_PLACE, &ave, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    ave /= (double)tot_num_pixels;

    // Calculate power spectrum
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------

    float k_factor = 1.5;
    float delta_k = run_globals.params.ReionPowerSpecDeltaK;
    float k_first_bin_ceil = delta_k;
    float k_max = delta_k * ReionGridDim;

    // Initialise arrays
    // ghetto counting (lookup how to do logs of arbitrary bases in c...)
    int num_bins = 0;
    float k_floor = 0.0;
    float k_ceil = k_first_bin_ceil;
    while (k_ceil < k_max) {
        num_bins++;
        k_floor = k_ceil;
        k_ceil *= k_factor;
    }

    //    printf("    NUM_BINS       = %d\n", NUM_BINS);

    double* p_box = malloc(sizeof(double) * num_bins);
    double* k_ave = malloc(sizeof(double) * num_bins);
    unsigned long long* in_bin_ct = malloc(sizeof(unsigned long long) * num_bins);

    for (int ii = 0; ii < num_bins; ii++) {
        p_box[ii] = 0.0;
        k_ave[ii] = 0.0;
        in_bin_ct[ii] = 0;
    }

    fftwf_complex* deldel_T = fftwf_alloc_complex(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);

    // Fill-up the real-space of the deldel box
    // Note: we include the V/N factor for the scaling after the fft
    float volume = powf(run_globals.params.BoxSize, 3);
    for (int ii = 0; ii < local_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++)
                ((float*)deldel_T)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] = (delta_T[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] / ave - 1) * volume / (float)tot_num_pixels;

    // Transform to k-space
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, (float*)deldel_T,
        deldel_T, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Now construct the power spectrum
    int HII_middle = ReionGridDim / 2;
    for (int n_x = 0; n_x < ReionGridDim; n_x++) {
        float k_x = n_x * delta_k;
        if (n_x > HII_middle)
            k_x = (n_x - ReionGridDim) * delta_k; // Wrap around for FFT convention

        for (int n_y = 0; n_y < ReionGridDim; n_y++) {
            float k_y = n_y * delta_k;
            if (n_y > HII_middle)
                k_y = (n_y - ReionGridDim) * delta_k;

            for (int n_z = 0; n_z <= HII_middle; n_z++) {
                float k_z = n_z * delta_k;

                k_mag = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

                // Now go through the k bins and update
                ct = 0;
                k_floor = 0.0;
                k_ceil = k_first_bin_ceil;
                while (k_ceil < k_max) {
                    // Check if we fal in this bin
                    if ((k_mag >= k_floor) && (k_mag < k_ceil)) {
                        in_bin_ct[ct]++;
                        p_box[ct] += pow(k_mag, 3) * pow(cabs(deldel_T[grid_index(n_x, n_y, n_z, ReionGridDim, INDEX_PADDED)]), 2) / (2.0 * M_PI * M_PI * volume);
                        // Note the 1/VOLUME factor, which turns this into a power density in k-space

                        k_ave[ct] += k_mag;
                        break;
                    }

                    ct++;
                    k_floor = k_ceil;
                    k_ceil *= k_factor;
                }
            }
        }
    } // End looping through k box

    // Malloc and store the result
    *ps = (float*)calloc(3 * num_bins, sizeof(float));

    // NOTE - previous ct ran from 1 (not zero) to NUM_BINS
    for (int ii = 0; ii < num_bins; ii++)
        if (in_bin_ct[ii] > 0) {
            (*ps)[0 + 3 * ii] = k_ave[ii] / (float)in_bin_ct[ii]; // Wavenumber
            (*ps)[1 + 3 * ii] = p_box[ii] / (float)in_bin_ct[ii]; // Power
            (*ps)[2 + 3 * ii] = p_box[ii] / (float)in_bin_ct[ii] / sqrt((float)in_bin_ct[ii]); // Error in power?
        }

    *ps_nbins = num_bins;
    *average_T = (float)ave;

    // Deallocate
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------

    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    fftwf_free(deldel_T);

    return 0;
}
