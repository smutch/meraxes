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

#include <assert.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "meraxes.h"
#include "misc_tools.h"

#include <cuda_runtime.h>
#include <cufft.h>

// These functions deal with any GPU exceptions, but should be called with the macros defined in the corresponding .hh
// file
__host__ void _throw_on_generic_error(bool check_failure,
                                      int implementation_code,
                                      const std::string& file,
                                      const std::string& func,
                                      int line)
{
  if (check_failure)
    throw(meraxes_cuda_exception(GENERIC_CUDA_ERROR_CODE, implementation_code, file, func, line));
}
__host__ void _throw_on_cuda_error(cudaError_t cuda_code,
                                   int implementation_code,
                                   const std::string& file,
                                   const std::string& func,
                                   int line)
{
  if (cuda_code != cudaSuccess)
    throw(meraxes_cuda_exception((int)cuda_code, implementation_code, file, func, line));
}
__host__ void _throw_on_cuFFT_error(cufftResult cufft_code,
                                    int implementation_code,
                                    const std::string& file,
                                    const std::string& func,
                                    int line)
{
  if (cufft_code != CUFFT_SUCCESS)
    throw(meraxes_cuda_exception((int)cufft_code, implementation_code, file, func, line));
}
__host__ void _check_for_cuda_error(int implementation_code, const std::string& file, const std::string& func, int line)
{
  try {
    cudaError_t cuda_code = cudaPeekAtLastError();
    if (cuda_code != cudaSuccess)
      throw(
        meraxes_cuda_exception((int)cuda_code, implementation_code, "CUDA error detected after ", file, func, line));
  } catch (const meraxes_cuda_exception& e) {
    e.process_exception();
  }
}
__host__ void _check_thread_sync(int implementation_code, const std::string& file, const std::string& func, int line)
{
  try {
    cudaError_t cuda_code = cudaDeviceSynchronize();
    if (cuda_code != cudaSuccess)
      throw(meraxes_cuda_exception(
        (int)cuda_code, implementation_code, "Threads not synchronised after ", file, func, line));
  } catch (const meraxes_cuda_exception& e) {
    e.process_exception();
  }
}
__host__ void _throw_on_global_error(const std::string& file, const std::string& func, int line)
{
  int error_code = 0;
  MPI_Allreduce(MPI_IN_PLACE, &error_code, 1, MPI_INT, MPI_MAX, run_globals.mpi_comm);
  if (error_code != 0)
    throw(meraxes_cuda_exception(0, meraxes_cuda_exception::GLOBAL, file, func, line));
}
__host__ void notify_of_global_error(int error_code)
{
  int result = (int)error_code;
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MAX, run_globals.mpi_comm);
}

static void dump_gpu_properties()
{
  // This should be called after the logic of init_CUDA() is complete.
  gpu_info info = *run_globals.gpu;

  mlog("\nGPU properties", MLOG_MESG);
  mlog("==============", MLOG_MESG);
  mlog("device = %s", MLOG_MESG, info.device);
  mlog("name = %s", MLOG_MESG, info.properties.name);
  mlog("flag_use_cuFFT = %d", MLOG_MESG, (int)info.flag_use_cuFFT);
  mlog("n_threads = %d", MLOG_MESG, info.n_threads);
  mlog("maxThreadsPerBlock = %d", MLOG_MESG, info.properties.maxThreadsPerBlock);
  mlog("maxThreadsPerMultiProcessor = %d", MLOG_MESG, info.properties.maxThreadsPerMultiProcessor);
  mlog("concurrentKernels = %d", MLOG_MESG, info.properties.concurrentKernels);
  mlog("==============", MLOG_MESG);
}

// Initialize device.  Called by init_gpu().
void init_CUDA()
{

  // TODO: Ask Greg about `throw_on_global_error`
  // TODO: Can I just use the NVIDIA sample code for error handling?

  try {
    // get the number of devices available to this rank
    int num_devices = 0;
    throw_on_cuda_error(cudaGetDeviceCount(&num_devices), meraxes_cuda_exception::INIT);

    int world_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // do your thang to assign this rank to a device
    throw_on_cuda_error(cudaSetDevice(world_rank % num_devices),
                        meraxes_cuda_exception::INIT); // alternate assignment between ranks

    // do a check to make sure that we have a working assigned device
    throw_on_cuda_error(cudaFree(0), meraxes_cuda_exception::INIT);
    throw_on_global_error();

    // Get the device assigned to this context
    throw_on_cuda_error(cudaGetDevice(&(run_globals.gpu->device)), meraxes_cuda_exception::INIT);

    // Get the properties of the device assigned to this context
    throw_on_cuda_error(cudaGetDeviceProperties(&(run_globals.gpu->properties), run_globals.gpu->device),
                        meraxes_cuda_exception::INIT);

    // Throw an exception if another rank has thrown one
    throw_on_global_error();
  } catch (const meraxes_cuda_exception& e) {
    e.process_exception();
  }

  // Set the number of threads to use.  Perhaps something
  //    more sophisticated should be done here eventually ...
  run_globals.gpu->n_threads = 256;

  // Send a status message to mlog
  if (run_globals.mpi_size == 1) {
    mlog("GPU context established on device %d (%s; %.1fGBs of global memory).",
         MLOG_MESG,
         run_globals.gpu->device,
         run_globals.gpu->properties.name,
         (float)(run_globals.gpu->properties.totalGlobalMem / (1024 * 1024 * 1024)));
  } else {
    // Report device details for each rank
    for (int i_rank = 0; i_rank < run_globals.mpi_size; i_rank++) {
      if (i_rank == run_globals.mpi_rank) {
        char* mps_pipe = NULL;
        if ((mps_pipe = secure_getenv("CUDA_MPS_PIPE_DIRECTORY")) != NULL) {
          mlog(
            "Context established on GPU device with MPS server: %s", MLOG_MESG | MLOG_ALLRANKS | MLOG_FLUSH, mps_pipe);
        } else {
          mlog("Context established on GPU device %d (%s; %.1fGBs of global memory).",
               MLOG_MESG | MLOG_ALLRANKS | MLOG_FLUSH,
               run_globals.gpu->device,
               run_globals.gpu->properties.name,
               (float)(run_globals.gpu->properties.totalGlobalMem / (1024 * 1024 * 1024)));
        }
      }
      MPI_Barrier(run_globals.mpi_comm);
    }
  }

  dump_gpu_properties();
}

// Call this function in kernels to put the GPU in an error state that can be caught after as an exception
//    This is not necessarily the best way, but it will do the job for now.  This is based on:
// https://devtalk.nvidia.com/default/topic/418479/how-to-trigger-a-cuda-error-from-inside-a-kernel/
__device__ void inline cause_cuda_error()
{
  int* adr = (int*)0xffffffff;
  *adr = 12;
}

// Convert a grid vector index to an (i,j,k)-triplet
// Pass i_x_start=0 for local indices, rather than global indices
__device__ void inline grid_index2indices(const int idx,
                                          const int dim,
                                          const int i_x_start,
                                          const int mode,
                                          int* ix,
                                          int* iy,
                                          int* iz)
{
  int n;
  int remainder = idx;
  switch (mode) {
    case INDEX_PADDED:
      n = (2 * (dim / 2 + 1));
      break;
    case INDEX_REAL:
      n = dim;
      break;
    case INDEX_COMPLEX_HERM:
      n = (dim / 2 + 1);
      break;
    default:
      cause_cuda_error();
      break;
  }
  (*iz) = remainder % n;
  remainder = (remainder - (*iz)) / n;
  (*iy) = remainder % dim;
  remainder = (remainder - (*iy)) / dim;
  (*ix) = remainder + i_x_start;
}

__device__ int inline grid_index_gpu(int i, int j, int k, int dim, int type)
{
  int ind;

  switch (type) {
    case INDEX_PADDED:
      ind = k + (2 * (dim / 2 + 1)) * (j + dim * i);
      break;
    case INDEX_REAL:
      ind = k + dim * (j + dim * i);
      break;
    case INDEX_COMPLEX_HERM:
      ind = k + (dim / 2 + 1) * (j + dim * i);
      break;
    default:
      cause_cuda_error();
      break;
  }

  return ind;
}

__device__ float inline k_mag_gpu(const int n_x, const int n_y, const int n_z, const int grid_dim, const float box_size)
{
  float delta_k = 2.0 * M_PI / box_size;
  int middle = grid_dim / 2;
  // Compute k_x
  float k_x;
  if (n_x > middle)
    k_x = (n_x - grid_dim) * delta_k;
  else
    k_x = n_x * delta_k;
  // Compute k_y
  float k_y;
  if (n_y > middle)
    k_y = (n_y - grid_dim) * delta_k;
  else
    k_y = n_y * delta_k;
  // Compute k_z
  float k_z;
  k_z = n_z * delta_k;
  return (sqrtf(k_x * k_x + k_y * k_y + k_z * k_z));
}

// Kernel to set all the values of a vector to a constant
__global__ void set_array_gpu(float* array, int n_real, float val)
{
  int i_real = blockIdx.x * blockDim.x + threadIdx.x;
  if (i_real < n_real)
    array[i_real] = val;
}

// Kernel to multiply a complex vector by a scalar
__global__ void complex_vector_times_scalar(Complex* vector, double scalar, int n_complex)
{
  int i_complex = blockIdx.x * blockDim.x + threadIdx.x;
  if (i_complex < n_complex) {
    vector[i_complex].x *= scalar;
    vector[i_complex].y *= scalar;
  }
}

// Kernel to perform aliasing sanity checks
__global__ void sanity_check_aliasing(Complex* grid, int grid_dim, int n_real, float val)
{
  int i_real = blockIdx.x * blockDim.x + threadIdx.x;
  if (i_real < n_real) {
    int i_x, i_y, i_z;
    grid_index2indices(
      i_real, grid_dim, 0, INDEX_REAL, &i_x, &i_y, &i_z); // pass i_x_start=0 'cause we want the local indices
    const int i_padded = grid_index_gpu(i_x, i_y, i_z, grid_dim, INDEX_PADDED);
    ((float*)grid)[i_padded] = fmaxf(((float*)grid)[i_padded], val);

    if ((val == 0.f) && (((float*)grid)[i_padded] < ABS_TOL))
      ((float*)grid)[i_padded] = 0.f;
  }
}

// Kernel to perform filtering convolution
__global__ void
filter_gpu(Complex* box, int grid_dim, int local_ix_start, int n_complex, float R, double box_size_in, int filter_type)
{
  int i_complex_herm = blockIdx.x * blockDim.x + threadIdx.x;
  if (i_complex_herm < n_complex) {

    // Map i_complex to global grid indices
    int n_x, n_y, n_z;
    grid_index2indices(i_complex_herm,
                       grid_dim,
                       local_ix_start,
                       INDEX_COMPLEX_HERM,
                       &n_x,
                       &n_y,
                       &n_z); // pass i_x_start becuase we need the global indices for k_mag

    // Compute k_mag for given grid indices
    float box_size = box_size_in;
    float k_mag = k_mag_gpu(n_x, n_y, n_z, grid_dim, box_size);

    // Calculate filter
    float kR = k_mag * R;
    float scalar = 0.f;
    int support = false;
    switch (filter_type) {
      case 0: // Real space top-hat
        scalar = (3.0 * (sinf(kR) / powf(kR, 3) - cosf(kR) / powf(kR, 2)));
        support = (kR > 1e-4);
        break;

      case 1:               // k-space top hat
        kR *= 0.413566994f; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
        scalar = 0.f;
        support = (kR > 1);
        break;

      case 2:         // Gaussian
        kR *= 0.643f; // Equates integrated volume to the real space top-hat
        scalar = powf(M_E, -kR * kR / 2.0f);
        support = true;
        break;
    }

    // Apply filter
    if (support) {
      box[i_complex_herm].x *= scalar;
      box[i_complex_herm].y *= scalar;
    }
  }
}

// Kernel to apply the logic of the main loop of find_HII_bubbles
__global__ void find_HII_bubbles_gpu_main_loop(const float redshift,
                                               const int n_real,
                                               const int flag_last_filter_step,
                                               const int ReionUVBFlag,
                                               const int Flag_IncludeRecombinations,
                                               const int ReionGridDim,
                                               const float R,
                                               const float M,
                                               const float ReionEfficiency,
                                               const float inv_pixel_volume,
                                               const float J_21_aux_constant,
                                               const double ReionGammaHaloBias,
                                               const double UnitMass_in_g,
                                               const double UnitTime_in_s,
                                               const double UnitLength_in_cm,
                                               const double Hubble_h,
                                               const double ReionNionPhotPerBary,
                                               const double Gamma_R_prefactor,
                                               float* xH,
                                               float* J_21,
                                               float* r_bubble,
                                               float* J_21_at_ionization,
                                               float* z_at_ionization,
                                               float* Gamma12,
                                               Complex* deltax_filtered_device,
                                               Complex* stars_filtered_device,
                                               Complex* sfr_filtered_device,
                                               Complex* N_rec_filtered_device)
{
  int i_real = blockIdx.x * blockDim.x + threadIdx.x;
  double Gamma_R_conversion = Gamma_R_prefactor * (UnitMass_in_g / UnitTime_in_s) *
                              pow(UnitLength_in_cm / Hubble_h, -3.) * ReionNionPhotPerBary /
                              PROTONMASS; // Convert pixel volume (Mpc/h)^3 -> (cm)^3

  // This will be reset below if Flag_IncludeRecombinations.
  double rec = 0.0;

  if (i_real < n_real) {
    int ix, iy, iz;
    grid_index2indices(
      i_real, ReionGridDim, 0, INDEX_REAL, &ix, &iy, &iz); // pass i_x_start=0 'cause we want the local indices
    const int i_padded = grid_index_gpu(ix, iy, iz, ReionGridDim, INDEX_PADDED);

    double density_over_mean = 1.0 + (double)((float*)deltax_filtered_device)[i_padded];

    double f_coll_stars = (double)((float*)stars_filtered_device)[i_padded] / (M * density_over_mean) * (4.0 / 3.0) *
                          M_PI * (R * R * R) * inv_pixel_volume;

    double sfr_density = (double)((float*)sfr_filtered_device)[i_padded] * inv_pixel_volume; // In internal units

    // Calculate the recombinations within the cell
    if (Flag_IncludeRecombinations)
      rec = (double)((float*)N_rec_filtered_device)[i_padded] / density_over_mean;

    float J_21_aux;
    if (ReionUVBFlag)
      J_21_aux = (float)(sfr_density * J_21_aux_constant);

    // Modified reionisation condition, including recombinations.
    if (f_coll_stars > (1.0 / ReionEfficiency) * (1. + rec)) // IONISED!!!!
    {
      // If it is the first crossing of the ionisation barrier for this cell (largest R), let's record J_21
      if (xH[i_real] > REL_TOL) {
        if (ReionUVBFlag)
          J_21[i_real] = J_21_aux;

        // Store the ionisation background and the reionisation redshift for each cell
        if (Flag_IncludeRecombinations) {
          Gamma12[i_real] = (float)(Gamma_R_prefactor * sfr_density * (UnitMass_in_g / UnitTime_in_s) *
                                    pow(UnitLength_in_cm / Hubble_h, -3.) * ReionNionPhotPerBary /
                                    PROTONMASS); // Convert pixel volume (Mpc/h)^3 -> (cm)^3
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
      xH[i_real] = (float)(1.0 - f_coll_stars * ReionEfficiency);
      if (xH[i_real] < 0.) {
        xH[i_real] = (float)0.;
      } else if (xH[i_real] > 1.0) {
        xH[i_real] = (float)1.;
      }
    }

    // Check if new ionisation
    float* z_in = z_at_ionization;
    if ((xH[i_real] < REL_TOL) && (z_in[i_real] < 0)) // New ionisation!
    {
      z_in[i_real] = (float)redshift;
      if (ReionUVBFlag)
        J_21_at_ionization[i_real] = J_21_aux * (float)ReionGammaHaloBias;
    }
  }
}

// vim:set et sw=2 ts=2:
