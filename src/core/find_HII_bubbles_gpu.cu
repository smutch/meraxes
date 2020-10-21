//==============================================================================
//
// This code was developed as part of the Astronomy Data and Computing Services
// (ADACS; https://adacs.org.au) 2017B Software Support program.
//
// Written by: Gregory B. Poole
// Date:       September 2017
//
// It is distributed under the MIT (Expat) License (see https://opensource.org/):
//
// Copyright (c) 2017 Astronomy Data and Computing Services (ADACS)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//==============================================================================

#include <complex.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <math.h>

#include "XRayHeatingFunctions.h"
#include "find_HII_bubbles.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "recombinations.h"
#include "reionization.h"

/*
 * ⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻
 * ------------------------------------------------------------
 * TODO LIST - @smutch (See also TODO items throughout the code)
 * ------------------------------------------------------------
 * - [X] Why are the unfiltered grids being copied to the GPU and then copied back each filtering step if using FFTW?
 * - [ ] Investigate the use of streams.
 * - [ ] The error handling is horrendously complicated. Is this really necessary?
 *
 * ⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻
 */

/*
 * This code is a re-write of the modified version of 21cmFAST used in Mutch et
 * al. (2016; Meraxes paper).  The original code was written by Andrei Mesinger
 * with additions as detailed in Sobacchi & Mesinger (2013abc).  Updates were
 * subsequently made by Simon Mutch & Paul Geil.
 */

// This is the CUDA-enabled version of find_HII_bubbles().  It uses cuFFT for
//    all FFTs if the USE_CUFFT compiler flag has been set.  Otherwise, it uses
//    fftw.  In both cases, everything else is done with the GPU.
void _find_HII_bubbles_gpu(const int snapshot, const bool flag_write_validation_output)
{
  // Fetch needed things from run_globals
  const MPI_Comm mpi_comm = run_globals.mpi_comm;
  const int mpi_rank = run_globals.mpi_rank;
  const double box_size = run_globals.params.BoxSize;
  const int ReionGridDim = run_globals.params.ReionGridDim;
  const int ReionUVBFlag = run_globals.params.ReionUVBFlag;
  const int Flag_IncludeRecombinations = run_globals.params.Flag_IncludeRecombinations;
  const double ReionEfficiency = run_globals.params.physics.ReionEfficiency;
  const double ReionNionPhotPerBary = run_globals.params.physics.ReionNionPhotPerBary;
  const double UnitLength_in_cm = run_globals.units.UnitLength_in_cm;
  const double UnitMass_in_g = run_globals.units.UnitMass_in_g;
  const double UnitTime_in_s = run_globals.units.UnitTime_in_s;
  const double ReionRBubbleMax = run_globals.params.physics.ReionRBubbleMax;
  const double ReionRBubbleMin = run_globals.params.physics.ReionRBubbleMin;
  const double ReionDeltaRFactor = run_globals.params.ReionDeltaRFactor;
  const double ReionGammaHaloBias = run_globals.params.physics.ReionGammaHaloBias;
  const double ReionAlphaUV = run_globals.params.physics.ReionAlphaUV;
  const double Hubble_h = run_globals.params.Hubble_h;
  // const double ReionEscapeFrac = run_globals.params.physics.ReionEscapeFrac;
  // grid parameters
  const ptrdiff_t* slabs_nix = run_globals.reion_grids.slab_nix;
  const ptrdiff_t* slabs_n_complex = run_globals.reion_grids.slab_n_complex;
  const ptrdiff_t* slabs_ix_start = run_globals.reion_grids.slab_ix_start;
  const int local_nix = (int)slabs_nix[mpi_rank];
  const int local_ix_start = (int)slabs_ix_start[mpi_rank];
  const int slab_n_complex = (int)(slabs_n_complex[mpi_rank]);
  const int slab_n_real = local_nix * ReionGridDim * ReionGridDim;
  // input grids
  float* deltax = run_globals.reion_grids.deltax; // real & padded
  float* stars = run_globals.reion_grids.stars;   // real & padded
  float* sfr = run_globals.reion_grids.sfr;       // real & padded
  // preallocated grids
  float* J_21 = run_globals.reion_grids.J_21;         // real
  float* r_bubble = run_globals.reion_grids.r_bubble; // real
  // output grids
  float* xH = run_globals.reion_grids.xH;                                 // real
  float* z_at_ionization = run_globals.reion_grids.z_at_ionization;       // real
  float* J_21_at_ionization = run_globals.reion_grids.J_21_at_ionization; // real
  // output values
  double* volume_weighted_global_xH = &(run_globals.reion_grids.volume_weighted_global_xH);
  double* mass_weighted_global_xH = &(run_globals.reion_grids.mass_weighted_global_xH);
  // a few needed constants
  const double pixel_volume = pow(box_size / (double)ReionGridDim, 3); // (Mpc/h)^3
  const double total_n_cells = pow((double)ReionGridDim, 3);
  const double inv_total_n_cells = 1. / total_n_cells;

  // Inhomogeneous recombinations {{
  float* Gamma12 = run_globals.reion_grids.Gamma12;

  float* N_rec = run_globals.reion_grids.N_rec;
  float* N_rec_prev = run_globals.reion_grids.N_rec_prev;

  const double redshift = run_globals.ZZ[snapshot];
  double prev_redshift;
  if (snapshot == 0) {
    prev_redshift = run_globals.ZZ[snapshot];
  } else {
    prev_redshift = run_globals.ZZ[snapshot - 1];
  }

  float zstep = (float)(prev_redshift - redshift);
  float fabs_dtdz = (float)fabs(dtdz((float)redshift) / Hubble_h);

  // Check that a valid filter option has been specified
  try {
    throw_on_generic_error(run_globals.params.ReionFilterType < 0 || run_globals.params.ReionFilterType > 2,
                           meraxes_cuda_exception::INVALID_FILTER);
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

  // Initialize device arrays
  cufftComplex* deltax_unfiltered_device = NULL;
  cufftComplex* stars_unfiltered_device = NULL;
  cufftComplex* sfr_unfiltered_device = NULL;
  cufftComplex* deltax_filtered_device = NULL;
  cufftComplex* stars_filtered_device = NULL;
  cufftComplex* sfr_filtered_device = NULL;
  cufftComplex* N_rec_unfiltered_device = NULL;
  cufftComplex* N_rec_filtered_device = NULL;
  float* xH_device = NULL;
  float* r_bubble_device = NULL;
  float* z_at_ionization_device = NULL;
  float* J_21_at_ionization_device = NULL;
  float* J_21_device = NULL;
  float* Gamma12_device = NULL;

  try {
    throw_on_cuda_error(cudaMalloc((void**)&deltax_unfiltered_device, sizeof(cufftComplex) * slab_n_complex),
                        meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&stars_unfiltered_device, sizeof(cufftComplex) * slab_n_complex),
                        meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&sfr_unfiltered_device, sizeof(cufftComplex) * slab_n_complex),
                        meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&deltax_filtered_device, sizeof(cufftComplex) * slab_n_complex),
                        meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&stars_filtered_device, sizeof(cufftComplex) * slab_n_complex),
                        meraxes_cuda_exception::MALLOC);
    throw_on_cuda_error(cudaMalloc((void**)&sfr_filtered_device, sizeof(cufftComplex) * slab_n_complex),
                        meraxes_cuda_exception::MALLOC);

    if (Flag_IncludeRecombinations) {
      throw_on_cuda_error(cudaMalloc((void**)&N_rec_unfiltered_device, sizeof(cufftComplex) * slab_n_complex),
                          meraxes_cuda_exception::MALLOC);
      throw_on_cuda_error(cudaMalloc((void**)&N_rec_filtered_device, sizeof(cufftComplex) * slab_n_complex),
                          meraxes_cuda_exception::MALLOC);
    }

    if (slab_n_real > 0) {
      throw_on_cuda_error(cudaMalloc((void**)&xH_device, sizeof(float) * slab_n_real), meraxes_cuda_exception::MALLOC);
      throw_on_cuda_error(cudaMalloc((void**)&r_bubble_device, sizeof(float) * slab_n_real),
                          meraxes_cuda_exception::MALLOC);
      throw_on_cuda_error(cudaMalloc((void**)&z_at_ionization_device, sizeof(float) * slab_n_real),
                          meraxes_cuda_exception::MALLOC);
      throw_on_cuda_error(cudaMalloc((void**)&J_21_at_ionization_device, sizeof(float) * slab_n_real),
                          meraxes_cuda_exception::MALLOC);
      if (ReionUVBFlag)
        throw_on_cuda_error(cudaMalloc((void**)&J_21_device, sizeof(float) * slab_n_real),
                            meraxes_cuda_exception::MALLOC);
      if (Flag_IncludeRecombinations) {
        throw_on_cuda_error(cudaMalloc((void**)&Gamma12_device, sizeof(float) * slab_n_real),
                            meraxes_cuda_exception::MALLOC);
      }
    }

    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

// If we're not using CUFFT, do the forward FFT first, before sending it to the device
#ifndef USE_CUFFT
  // The following are only needed if we are using FFTW
  Complex* deltax_filtered =
    (Complex*)run_globals.reion_grids.deltax_filtered; // complex TODO: Check the consistancy of Complex instances
  Complex* stars_filtered =
    (Complex*)run_globals.reion_grids.stars_filtered; // complex TODO: Check the consistancy of Complex instances
  Complex* sfr_filtered =
    (Complex*)run_globals.reion_grids.sfr_filtered; // complex TODO: Check the consistancy of Complex instances

  Complex* N_rec_filtered = NULL;
  if (Flag_IncludeRecombinations)
    N_rec_filtered =
      (Complex*)run_globals.reion_grids.N_rec_filtered; // complex TODO: Check the consistancy of Complex instances

  // Forward fourier transform to obtain k-space fields
  fftwf_complex* deltax_unfiltered = (fftwf_complex*)deltax; // WATCH OUT!
  fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(
    ReionGridDim, ReionGridDim, ReionGridDim, deltax, deltax_unfiltered, mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  fftwf_complex* stars_unfiltered = (fftwf_complex*)stars; // WATCH OUT!
  plan = fftwf_mpi_plan_dft_r2c_3d(
    ReionGridDim, ReionGridDim, ReionGridDim, stars, stars_unfiltered, mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  fftwf_complex* sfr_unfiltered = (fftwf_complex*)sfr; // WATCH OUT!
  plan =
    fftwf_mpi_plan_dft_r2c_3d(ReionGridDim, ReionGridDim, ReionGridDim, sfr, sfr_unfiltered, mpi_comm, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  if (Flag_IncludeRecombinations) {
    fftwf_complex* N_rec_unfiltered = (fftwf_complex*)N_rec_prev; // WATCH OUT!
    plan = fftwf_mpi_plan_dft_r2c_3d(
      ReionGridDim, ReionGridDim, ReionGridDim, N_rec_prev, N_rec_unfiltered, mpi_comm, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
  }
#endif

  // Perform host -> device transfer of input grids (note that these grids are k-space if we are using FFTW
  //    but are real-space if we are using CUFFT.  They will be transformed once on the device if the latter.
  try {
    throw_on_cuda_error(
      cudaMemcpy(deltax_unfiltered_device, deltax, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
      meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(
      cudaMemcpy(stars_unfiltered_device, stars, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
      meraxes_cuda_exception::MEMCPY);
    throw_on_cuda_error(
      cudaMemcpy(sfr_unfiltered_device, sfr, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
      meraxes_cuda_exception::MEMCPY);

    if (Flag_IncludeRecombinations)
      throw_on_cuda_error(
        cudaMemcpy(N_rec_unfiltered_device, N_rec_prev, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
        meraxes_cuda_exception::MEMCPY);

    if (slab_n_real) {
      throw_on_cuda_error(
        cudaMemcpy(z_at_ionization_device, z_at_ionization, sizeof(float) * slab_n_real, cudaMemcpyHostToDevice),
        meraxes_cuda_exception::MEMCPY);
      throw_on_cuda_error(
        cudaMemcpy(J_21_at_ionization_device, J_21_at_ionization, sizeof(float) * slab_n_real, cudaMemcpyHostToDevice),
        meraxes_cuda_exception::MEMCPY);

      if (Flag_IncludeRecombinations) {
        throw_on_cuda_error(cudaMemcpy(Gamma12_device, Gamma12, sizeof(float) * slab_n_real, cudaMemcpyHostToDevice),
                            meraxes_cuda_exception::MEMCPY);
      }
    }

    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

// If we're using CUFFT, perform the forward FFT now that the data is on the device
#ifdef USE_CUFFT
  // Initialize cuFFT
  cufftHandle plan;
  try {
    throw_on_cuFFT_error(cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_R2C),
                         meraxes_cuda_exception::CUFFT_CREATE_PLAN);

    // depreciated in CUFFT > v9.1
    // throw_on_cuFFT_error(cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL),
    // meraxes_cuda_exception::CUFFT_SET_COMPATIBILITY);

    // Perform FFTs
    throw_on_cuFFT_error(cufftExecR2C(plan, (cufftReal*)deltax_unfiltered_device, deltax_unfiltered_device),
                         meraxes_cuda_exception::CUFFT_R2C);
    throw_on_cuFFT_error(cufftExecR2C(plan, (cufftReal*)stars_unfiltered_device, stars_unfiltered_device),
                         meraxes_cuda_exception::CUFFT_R2C);
    throw_on_cuFFT_error(cufftExecR2C(plan, (cufftReal*)sfr_unfiltered_device, sfr_unfiltered_device),
                         meraxes_cuda_exception::CUFFT_R2C);

    if (Flag_IncludeRecombinations)
      throw_on_cuFFT_error(cufftExecR2C(plan, (cufftReal*)N_rec_unfiltered_device, N_rec_unfiltered_device),
                           meraxes_cuda_exception::CUFFT_R2C);

    // Clean-up
    throw_on_cuFFT_error(cufftDestroy(plan), meraxes_cuda_exception::CUFFT_PLAN_DESTROY);

    // Make sure that the device has synchronized
    throw_on_cuda_error(cudaThreadSynchronize(), meraxes_cuda_exception::SYNC);

    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }
#endif

  // Initialize GPU block and thread count
  int threads = run_globals.gpu->n_threads;
  int grid_complex = (slab_n_complex + (threads - 1)) / threads;
  int grid_real = (slab_n_real + (threads - 1)) / threads;

  // mlog("threads = %d", MLOG_ALLRANKS|MLOG_MESG, threads);
  // mlog("slab_n_complex = %d", MLOG_ALLRANKS|MLOG_MESG, slab_n_complex);
  // mlog("grid_complex = %d", MLOG_ALLRANKS|MLOG_MESG, grid_complex);
  // mlog("slab_n_real = %d", MLOG_ALLRANKS|MLOG_MESG, slab_n_real);
  // mlog("grid_real = %d", MLOG_ALLRANKS|MLOG_MESG|MLOG_FLUSH, grid_real);

  MPI_Barrier(run_globals.mpi_comm);

  // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  try {
    throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(
                            deltax_unfiltered_device, inv_total_n_cells, slab_n_complex)),
                          meraxes_cuda_exception::KERNEL_CMPLX_AX);
    throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(
                            stars_unfiltered_device, inv_total_n_cells, slab_n_complex)),
                          meraxes_cuda_exception::KERNEL_CMPLX_AX);
    throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(
                            sfr_unfiltered_device, inv_total_n_cells, slab_n_complex)),
                          meraxes_cuda_exception::KERNEL_CMPLX_AX);

    if (Flag_IncludeRecombinations)
      throw_on_kernel_error((complex_vector_times_scalar<<<grid_complex, threads>>>(
                              N_rec_unfiltered_device, inv_total_n_cells, slab_n_complex)),
                            meraxes_cuda_exception::KERNEL_CMPLX_AX);

    check_thread_sync(meraxes_cuda_exception::KERNEL_CMPLX_AX);
    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

  // Initialize a few of the output grids
  try {
    if (slab_n_real > 0) {
      throw_on_kernel_error((set_array_gpu<<<grid_real, threads>>>(xH_device, slab_n_real, 1.f)),
                            meraxes_cuda_exception::KERNEL_SET_ARRAY);
      throw_on_kernel_error((set_array_gpu<<<grid_real, threads>>>(r_bubble_device, slab_n_real, 0.f)),
                            meraxes_cuda_exception::KERNEL_SET_ARRAY);
      if (ReionUVBFlag)
        throw_on_kernel_error((set_array_gpu<<<grid_real, threads>>>(J_21_device, slab_n_real, 0.f)),
                              meraxes_cuda_exception::KERNEL_SET_ARRAY);
    }
    check_thread_sync(meraxes_cuda_exception::KERNEL_SET_ARRAY);
    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

  // This parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations. It probably doesn't matter for larger cell sizes.
  double cell_length_factor = L_FACTOR;
  if ((box_size / (double)ReionGridDim) < 1.0) // Fairly arbitrary length based on 2 runs Sobacchi did
    cell_length_factor = 1.0;

// Initialize inverse FFTs
#ifdef USE_CUFFT
  try {
    throw_on_cuFFT_error(cufftPlan3d(&plan, ReionGridDim, ReionGridDim, ReionGridDim, CUFFT_C2R),
                         meraxes_cuda_exception::CUFFT_C2R);
    // depreciated in CUFFT > v9.1
    // throw_on_cuFFT_error(cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL),
    // meraxes_cuda_exception::CUFFT_SET_COMPATIBILITY); Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }
#else
  fftwf_plan plan_deltax = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                     ReionGridDim,
                                                     ReionGridDim,
                                                     (fftwf_complex*)deltax_filtered,
                                                     (float*)deltax_filtered,
                                                     mpi_comm,
                                                     FFTW_ESTIMATE);
  fftwf_plan plan_stars = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                    ReionGridDim,
                                                    ReionGridDim,
                                                    (fftwf_complex*)stars_filtered,
                                                    (float*)stars_filtered,
                                                    mpi_comm,
                                                    FFTW_ESTIMATE);
  fftwf_plan plan_sfr = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                  ReionGridDim,
                                                  ReionGridDim,
                                                  (fftwf_complex*)sfr_filtered,
                                                  (float*)sfr_filtered,
                                                  mpi_comm,
                                                  FFTW_ESTIMATE);

  fftwf_plan plan_N_rec;
  if (Flag_IncludeRecombinations)
    plan_N_rec = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                           ReionGridDim,
                                           ReionGridDim,
                                           (fftwf_complex*)N_rec_filtered,
                                           (float*)N_rec_filtered,
                                           mpi_comm,
                                           FFTW_ESTIMATE);
#endif

  // Loop through filter radii
  double R = fmin(ReionRBubbleMax, L_FACTOR * box_size); // Mpc/h
  bool flag_last_filter_step = false;
  int i_R = 0;
  while (!flag_last_filter_step) {
    i_R++;

    // check to see if this is our last filtering step
    if (((R / ReionDeltaRFactor) <= (cell_length_factor * box_size / (double)ReionGridDim)) ||
        ((R / ReionDeltaRFactor) <= ReionRBubbleMin)) {
      flag_last_filter_step = true;
      R = cell_length_factor * box_size / (double)ReionGridDim;
    }

    mlog(".", MLOG_CONT);

    // Create working copies of the k-space grids
    try {
      throw_on_cuda_error(
        cudaMemcpy(
          deltax_filtered_device, deltax_unfiltered_device, sizeof(Complex) * slab_n_complex, cudaMemcpyDeviceToDevice),
        meraxes_cuda_exception::MEMCPY);
      throw_on_cuda_error(
        cudaMemcpy(
          stars_filtered_device, stars_unfiltered_device, sizeof(Complex) * slab_n_complex, cudaMemcpyDeviceToDevice),
        meraxes_cuda_exception::MEMCPY);
      throw_on_cuda_error(
        cudaMemcpy(
          sfr_filtered_device, sfr_unfiltered_device, sizeof(Complex) * slab_n_complex, cudaMemcpyDeviceToDevice),
        meraxes_cuda_exception::MEMCPY);

      if (Flag_IncludeRecombinations)
        throw_on_cuda_error(
          cudaMemcpy(
            N_rec_filtered_device, N_rec_unfiltered_device, sizeof(Complex) * slab_n_complex, cudaMemcpyDeviceToDevice),
          meraxes_cuda_exception::MEMCPY);

      // Throw an exception if another rank has thrown one
      throw_on_global_error();
    } catch (const meraxes_cuda_exception e) {
      e.process_exception();
    }

    // Perform convolution
    if (!flag_last_filter_step) {
      try {
        throw_on_kernel_error((filter_gpu<<<grid_complex, threads>>>(deltax_filtered_device,
                                                                     ReionGridDim,
                                                                     local_ix_start,
                                                                     slab_n_complex,
                                                                     R,
                                                                     box_size,
                                                                     run_globals.params.ReionFilterType)),
                              meraxes_cuda_exception::KERNEL_FILTER);
        throw_on_kernel_error((filter_gpu<<<grid_complex, threads>>>(stars_filtered_device,
                                                                     ReionGridDim,
                                                                     local_ix_start,
                                                                     slab_n_complex,
                                                                     R,
                                                                     box_size,
                                                                     run_globals.params.ReionFilterType)),
                              meraxes_cuda_exception::KERNEL_FILTER);
        throw_on_kernel_error((filter_gpu<<<grid_complex, threads>>>(sfr_filtered_device,
                                                                     ReionGridDim,
                                                                     local_ix_start,
                                                                     slab_n_complex,
                                                                     R,
                                                                     box_size,
                                                                     run_globals.params.ReionFilterType)),
                              meraxes_cuda_exception::KERNEL_FILTER);

        if (Flag_IncludeRecombinations)
          throw_on_kernel_error((filter_gpu<<<grid_complex, threads>>>(N_rec_filtered_device,
                                                                       ReionGridDim,
                                                                       local_ix_start,
                                                                       slab_n_complex,
                                                                       R,
                                                                       box_size,
                                                                       run_globals.params.ReionFilterType)),
                                meraxes_cuda_exception::KERNEL_FILTER);

        check_thread_sync(meraxes_cuda_exception::KERNEL_FILTER);
        // Throw an exception if another rank has thrown one
        throw_on_global_error();
      } catch (const meraxes_cuda_exception e) {
        e.process_exception();
      }
    }

    // inverse fourier transform back to real space
    try {
#ifdef USE_CUFFT
      throw_on_cuFFT_error(
        cufftExecC2R(plan, (cufftComplex*)deltax_filtered_device, (cufftReal*)deltax_filtered_device),
        meraxes_cuda_exception::CUFFT_C2R);
      throw_on_cuFFT_error(cufftExecC2R(plan, (cufftComplex*)stars_filtered_device, (cufftReal*)stars_filtered_device),
                           meraxes_cuda_exception::CUFFT_C2R);
      throw_on_cuFFT_error(cufftExecC2R(plan, (cufftComplex*)sfr_filtered_device, (cufftReal*)sfr_filtered_device),
                           meraxes_cuda_exception::CUFFT_C2R);

      if (Flag_IncludeRecombinations)
        throw_on_cuFFT_error(
          cufftExecC2R(plan, (cufftComplex*)N_rec_filtered_device, (cufftReal*)N_rec_filtered_device),
          meraxes_cuda_exception::CUFFT_C2R);

#else
      throw_on_cuda_error(
        cudaMemcpy(deltax_filtered, deltax_filtered_device, sizeof(float) * 2 * slab_n_complex, cudaMemcpyDeviceToHost),
        meraxes_cuda_exception::MEMCPY);
      fftwf_execute(plan_deltax);
      throw_on_cuda_error(
        cudaMemcpy(deltax_filtered_device, deltax_filtered, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
        meraxes_cuda_exception::MEMCPY);

      throw_on_cuda_error(
        cudaMemcpy(stars_filtered, stars_filtered_device, sizeof(float) * 2 * slab_n_complex, cudaMemcpyDeviceToHost),
        meraxes_cuda_exception::MEMCPY);
      fftwf_execute(plan_stars);
      throw_on_cuda_error(
        cudaMemcpy(stars_filtered_device, stars_filtered, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
        meraxes_cuda_exception::MEMCPY);

      throw_on_cuda_error(
        cudaMemcpy(sfr_filtered, sfr_filtered_device, sizeof(float) * 2 * slab_n_complex, cudaMemcpyDeviceToHost),
        meraxes_cuda_exception::MEMCPY);
      fftwf_execute(plan_sfr);
      throw_on_cuda_error(
        cudaMemcpy(sfr_filtered_device, sfr_filtered, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
        meraxes_cuda_exception::MEMCPY);

      if (Flag_IncludeRecombinations) {
        throw_on_cuda_error(
          cudaMemcpy(N_rec_filtered, N_rec_filtered_device, sizeof(float) * 2 * slab_n_complex, cudaMemcpyDeviceToHost),
          meraxes_cuda_exception::MEMCPY);
        fftwf_execute(plan_N_rec);
        throw_on_cuda_error(
          cudaMemcpy(N_rec_filtered_device, N_rec_filtered, sizeof(float) * 2 * slab_n_complex, cudaMemcpyHostToDevice),
          meraxes_cuda_exception::MEMCPY);
      }
#endif
      // Throw an exception if another rank has thrown one
      throw_on_global_error();
    } catch (const meraxes_cuda_exception e) {
      e.process_exception();
    }

    // Perform sanity checks to account for aliasing effects
    try {
      if (slab_n_real > 0) {
        throw_on_kernel_error((sanity_check_aliasing<<<grid_real, threads>>>(
                                deltax_filtered_device, ReionGridDim, slab_n_real, -1.f + REL_TOL)),
                              meraxes_cuda_exception::KERNEL_CHECK);
        throw_on_kernel_error(
          (sanity_check_aliasing<<<grid_real, threads>>>(stars_filtered_device, ReionGridDim, slab_n_real, 0.f)),
          meraxes_cuda_exception::KERNEL_CHECK);
        throw_on_kernel_error(
          (sanity_check_aliasing<<<grid_real, threads>>>(sfr_filtered_device, ReionGridDim, slab_n_real, 0.f)),
          meraxes_cuda_exception::KERNEL_CHECK);

        if (Flag_IncludeRecombinations)
          throw_on_kernel_error(
            (sanity_check_aliasing<<<grid_real, threads>>>(N_rec_filtered_device, ReionGridDim, slab_n_real, 0.f)),
            meraxes_cuda_exception::KERNEL_CHECK);
      }
      check_thread_sync(meraxes_cuda_exception::KERNEL_CHECK);

      // Throw an exception if another rank has thrown one
      throw_on_global_error();
    } catch (const meraxes_cuda_exception e) {
      e.process_exception();
    }

    // Main loop through the box...
    const double J_21_aux_constant = (1.0 + redshift) * (1.0 + redshift) / (4.0 * M_PI) * ReionAlphaUV * PLANCK *
                                     1e21 // * ReionEscapeFrac
                                     * R * UnitLength_in_cm * ReionNionPhotPerBary / PROTONMASS * UnitMass_in_g /
                                     pow(UnitLength_in_cm, 3) / UnitTime_in_s;
    const double inv_pixel_volume = 1. / pixel_volume;

    double Gamma_R_prefactor = 1.0;
    if (Flag_IncludeRecombinations) {
      Gamma_R_prefactor = (1.0 + redshift) * (1.0 + redshift) * R * (UnitLength_in_cm / Hubble_h) * SIGMA_HI *
                          ReionAlphaUV / (ReionAlphaUV + 2.75) / 1.0e-12; // Converting R h^-1 to R.
    }

    try {
      if (slab_n_real > 0) {
        throw_on_kernel_error((find_HII_bubbles_gpu_main_loop<<<grid_real, threads>>>(redshift,
                                                                                      slab_n_real,
                                                                                      flag_last_filter_step,
                                                                                      ReionUVBFlag,
                                                                                      Flag_IncludeRecombinations,
                                                                                      ReionGridDim,
                                                                                      R,
                                                                                      RtoM(R),
                                                                                      ReionEfficiency,
                                                                                      inv_pixel_volume,
                                                                                      J_21_aux_constant,
                                                                                      ReionGammaHaloBias,
                                                                                      UnitMass_in_g,
                                                                                      UnitTime_in_s,
                                                                                      UnitLength_in_cm,
                                                                                      Hubble_h,
                                                                                      ReionNionPhotPerBary,
                                                                                      Gamma_R_prefactor,
                                                                                      xH_device,
                                                                                      J_21_device,
                                                                                      r_bubble_device,
                                                                                      J_21_at_ionization_device,
                                                                                      z_at_ionization_device,
                                                                                      Gamma12_device,
                                                                                      deltax_filtered_device,
                                                                                      stars_filtered_device,
                                                                                      sfr_filtered_device,
                                                                                      N_rec_unfiltered_device)),
                              meraxes_cuda_exception::KERNEL_MAIN_LOOP);
      }
      check_thread_sync(meraxes_cuda_exception::KERNEL_MAIN_LOOP);
      // Throw an exception if another rank has thrown one
      throw_on_global_error();
    } catch (const meraxes_cuda_exception e) {
      e.process_exception();
    }

    R /= ReionDeltaRFactor;
  }

// Clean-up FFT plan(s)
#ifdef USE_CUFFT
  try {
    throw_on_cuFFT_error(cufftDestroy(plan), meraxes_cuda_exception::CUFFT_PLAN_DESTROY);
    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }
#else
  fftwf_destroy_plan(plan_deltax);
  fftwf_destroy_plan(plan_stars);
  fftwf_destroy_plan(plan_sfr);

  if (Flag_IncludeRecombinations)
    fftwf_destroy_plan(plan_N_rec);
#endif

  // Perform device -> host transfer
  try {
    if (slab_n_real > 0) {
      throw_on_cuda_error(cudaMemcpy((void*)xH, (void*)xH_device, sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),
                          meraxes_cuda_exception::MEMCPY);
      throw_on_cuda_error(
        cudaMemcpy((void*)r_bubble, (void*)r_bubble_device, sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),
        meraxes_cuda_exception::MEMCPY);
      if (ReionUVBFlag)
        throw_on_cuda_error(
          cudaMemcpy((void*)J_21, (void*)J_21_device, sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),
          meraxes_cuda_exception::MEMCPY);
      throw_on_cuda_error(
        cudaMemcpy(
          (void*)z_at_ionization, (void*)z_at_ionization_device, sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),
        meraxes_cuda_exception::MEMCPY);
      throw_on_cuda_error(cudaMemcpy((void*)J_21_at_ionization,
                                     (void*)J_21_at_ionization_device,
                                     sizeof(float) * slab_n_real,
                                     cudaMemcpyDeviceToHost),
                          meraxes_cuda_exception::MEMCPY);

      if (Flag_IncludeRecombinations) {
        throw_on_cuda_error(
          cudaMemcpy((void*)Gamma12, (void*)Gamma12_device, sizeof(float) * slab_n_real, cudaMemcpyDeviceToHost),
          meraxes_cuda_exception::MEMCPY);
      }
    }
    throw_on_cuda_error(
      cudaMemcpy(
        (void*)deltax, (void*)deltax_filtered_device, sizeof(float) * 2 * slab_n_complex, cudaMemcpyDeviceToHost),
      meraxes_cuda_exception::MEMCPY);
    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

  // Clean-up device
  try {
    throw_on_cuda_error(cudaFree(deltax_unfiltered_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(stars_unfiltered_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(sfr_unfiltered_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(deltax_filtered_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(stars_filtered_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(sfr_filtered_device), meraxes_cuda_exception::FREE);

    if (Flag_IncludeRecombinations)
      throw_on_cuda_error(cudaFree(N_rec_filtered_device), meraxes_cuda_exception::FREE);

    throw_on_cuda_error(cudaFree(xH_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(r_bubble_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(z_at_ionization_device), meraxes_cuda_exception::FREE);
    throw_on_cuda_error(cudaFree(J_21_at_ionization_device), meraxes_cuda_exception::FREE);

    if (ReionUVBFlag)
      throw_on_cuda_error(cudaFree(J_21_device), meraxes_cuda_exception::FREE);

    if (Flag_IncludeRecombinations) {
      throw_on_cuda_error(cudaFree(Gamma12_device), meraxes_cuda_exception::FREE);
    }

    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception e) {
    e.process_exception();
  }

  // Find the volume and mass weighted neutral fractions
  // TODO: The deltax grid will have rounding errors from forward and reverse
  //       FFT. Should cache deltax slabs prior to ffts and reuse here.
  *volume_weighted_global_xH = 0.0;
  *mass_weighted_global_xH = 0.0;
  double mass_weight = 0.0;

  // Calculate neutral fractions.
  // TODO: A parallel reduction could be done for this before results are off-loaded
  //       from the GPU.
  int ix, iy, iz;
  for (ix = 0; ix < local_nix; ix++)
    for (iy = 0; iy < ReionGridDim; iy++)
      for (iz = 0; iz < ReionGridDim; iz++) {
        const int i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);
        const int i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
        double density_over_mean = 1.0 + (double)((float*)deltax)[i_padded];
        *volume_weighted_global_xH += (double)xH[i_real];
        *mass_weighted_global_xH += (double)(xH[i_real]) * density_over_mean;
        mass_weight += density_over_mean;
      }
  MPI_Allreduce(MPI_IN_PLACE, volume_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, mass_weighted_global_xH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, &mass_weight, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  *volume_weighted_global_xH *= inv_total_n_cells;
  *mass_weighted_global_xH /= mass_weight;

  if (Flag_IncludeRecombinations) {
    // Store the resultant recombination grid
    for (int ix = 0; ix < local_nix; ix++)
      for (int iy = 0; iy < ReionGridDim; iy++)
        for (int iz = 0; iz < ReionGridDim; iz++) {
          const int i_padded = grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED);
          const int i_real = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

          const double density_over_mean = 1.0 + deltax[i_padded];
          const float z_eff = (float)((1. + redshift) * pow(density_over_mean, 1.0 / 3.0) - 1);
          const float dNrec =
            splined_recombination_rate(z_eff, Gamma12[i_real]) * fabs_dtdz * zstep * (1. - xH[i_real]);
          N_rec[i_padded] += dNrec;
        }
  }
}

// vim:set et sw=2 ts=2:
