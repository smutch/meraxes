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

#ifndef _MERAXES_GPU_HH
#define _MERAXES_GPU_HH

#ifdef __cplusplus
#ifdef __NVCC__
#include <cuda_runtime.h>
#include <cufft.h>

#define cuda_check _cuda_check(cudaError_t error, __FILE__, __FUNCTION__, __LINE__)
#define cuda_check_last _cuda_check(cudaGetLastError(), __FILE__, __FUNCTION__, __LINE__)
__global__ void _cuda_check(cudaError error, const char* file, const char* function, const int line_number);

#define cufft_check _cufft_check(cufftResult_t result, __FILE__, __FUNCTION__, __LINE__)
__global__ void _cufft_check(cufftResult_t result, const char* file, const char* function, const int line_number);

// These routines are needed by the kernels
__device__ void inline grid_index2indices(const int idx,
                                          const int dim,
                                          const int local_ix_start,
                                          const int mode,
                                          int* ix,
                                          int* iy,
                                          int* iz);
__device__ int inline grid_index_gpu(int i, int j, int k, int dim, int type);
__device__ float inline k_mag_gpu(const int n_x, const int n_y, const int n_z, const int dim, const float box_size);

// Kernel definitions
__global__ void set_array_gpu(float* array, int n_real, float val);
__global__ void sanity_check_aliasing(Complex* grid, int grid_dim, int n_real, float val);
__global__ void complex_vector_times_scalar(Complex* vector, double scalar, int n_complex);
__global__ void
filter_gpu(Complex* box, int grid_dim, int local_ix_start, int n_complex, float R, double box_size_in, int filter_type);
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
                                               Complex* N_rec_filtered_device);

#endif // __NVCC__
#endif // __cplusplus
