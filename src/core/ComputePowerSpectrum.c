#include "meraxes.h"
#include <assert.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>

/*
 * A generic function to compute the 21cm PS of any field. Algorithm taken from delta_T.c
 * from 21cmFAST, but generalised to take any input cubic field.
 * Written by Bradley Greig
 * 
 */

void Compute_PS(int snapshot, int field)
{
    // Fields of interest:
    // 1. Density field (already a fluctuation)
    // 2. 21cm brightness temperature (not a fluctuating quantity)

    float* deltax = run_globals.reion_grids.deltax;
    float* delta_T = run_globals.reion_grids.delta_T;

    fftwf_complex* deldel_ps = fftwf_alloc_complex(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);
    
    float volume = powf(run_globals.params.BoxSize, 3);
    
    double box_size = run_globals.params.BoxSize; // Mpc/h
    int ReionGridDim = run_globals.params.ReionGridDim;
    double total_n_cells = pow((double)ReionGridDim, 3);
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);

    int ii, jj, kk, n_x, n_y, n_z;

    double ave;

    if(field==1) {

        mlog("Calculating the 21cm power spectrum (dimensional, i.e mK^2)",MLOG_MESG);

        ave = 0.0;
        for (ii=0; ii<local_nix; ii++){
            for (jj=0; jj<ReionGridDim; jj++){
                for (kk=0; kk<ReionGridDim; kk++){
                    ave += delta_T[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)];
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &ave, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);

        ave /= total_n_cells;

        for (int ii = 0; ii < local_nix; ii++) {
            for (int jj = 0; jj < ReionGridDim; jj++) {
                for (int kk = 0; kk < ReionGridDim; kk++) {
                    ((float*)deldel_ps)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] = (delta_T[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] / ave - 1) * volume / (float)total_n_cells;
                    ((float*)deldel_ps)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] *= ave;
                }
            }
        }        

    }
    else if (field==0) {

        mlog("Calculating the matter power spectrum",MLOG_MESG);

        for (int ii = 0; ii < local_nix; ii++) {
            for (int jj = 0; jj < ReionGridDim; jj++) {
                for (int kk = 0; kk < ReionGridDim; kk++) {
       	            ((float*)deldel_ps)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] = deltax[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] * volume / (float)total_n_cells;
       	       	}
       	    }
       	}

    }

    else {
        mlog("Not a valid argument for computing the power spectrum",MLOG_MESG);
        mlog("Options are: 0: Matter power spectrum, 1: 21cm power spectrum",MLOG_MESG);

        exit(0);
    }


    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, (float *)deldel_ps, deldel_ps, run_globals.mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Calculate power spectrum
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------

    unsigned long long ct;

    float k_mag;

    float k_factor = 1.35;
    float delta_k = 2. * M_PI / box_size;
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

    double* p_box = malloc(sizeof(double) * num_bins);
    double* k_ave = malloc(sizeof(double) * num_bins);
    unsigned long long* in_bin_ct = malloc(sizeof(unsigned long long) * num_bins);

    for (int ii = 0; ii < num_bins; ii++) {
        p_box[ii] = 0.0;
        k_ave[ii] = 0.0;
        in_bin_ct[ii] = 0;
    }

    // Co-eval box, so should sample the entire cube
        
    // now construct the power spectrum file
    int HII_middle = ReionGridDim / 2;

    for (n_x=0; n_x<ReionGridDim; n_x++) {
        float k_x = n_x * delta_k;
        if (n_x > HII_middle) {
            k_x = ( n_x - ReionGridDim ) * delta_k;  // wrap around for FFT convention
        }
        
        for (n_y=0; n_y < ReionGridDim; n_y++){
            float k_y = n_y * delta_k;        
            if (n_y > HII_middle)
                k_y = ( n_y - ReionGridDim ) * delta_k;
                    
            for (n_z=0; n_z <= HII_middle; n_z++){
                float k_z = n_z * delta_k;
                        
                k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                        
                // now go through the k bins and update
                ct = 0;
                k_floor = 0;
                k_ceil = k_first_bin_ceil;
                
                while (k_ceil < k_max){
                    // check if we fall in this bin
                    if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                        in_bin_ct[ct]++;
                        p_box[ct] += pow(k_mag,3) * pow( cabs(deldel_ps[grid_index(n_x, n_y, n_z, ReionGridDim, INDEX_PADDED)]), 2.)/(2.0 * M_PI * M_PI * volume);
                        // note the 1/VOLUME factor, which turns this into a power density in k-space
                            
                        k_ave[ct] += k_mag;
                        break;
                    }
                            
                    ct++;
                    k_floor=k_ceil;
                    k_ceil*=k_factor;
                }
            }
        }
    } // end looping through k box
 
    // NOTE - previous ct ran from 1 (not zero) to NUM_BINS
    for (int ii = 0; ii < num_bins; ii++)
        if (in_bin_ct[ii] > 0) {
            k_ave[ii] = k_ave[ii] / (float)in_bin_ct[ii]; // Wavenumber
            p_box[ii] = p_box[ii] / (float)in_bin_ct[ii]; // Power
//           = p_box[ii] / (float)in_bin_ct[ii] / sqrt((float)in_bin_ct[ii]); // Error in power?
            mlog("PS Data: ii = %d k_ave = %e p_box = %e",MLOG_MESG,ii,k_ave[ii],p_box[ii]);            
        }

    // Deallocate
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------

    free(p_box);
    free(k_ave);
    free(in_bin_ct);
    fftwf_free(deldel_ps);

}
