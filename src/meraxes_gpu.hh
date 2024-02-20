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

// Host-side exception-handling routines
__host__ void _throw_on_generic_error(bool check_failure,
                                      int implementation_code,
                                      const std::string& file,
                                      const std::string& func,
                                      int line);
__host__ void _throw_on_cuda_error(cudaError_t cuda_code,
                                   int implementation_code,
                                   const std::string& file,
                                   const std::string& func,
                                   int line);
__host__ void _throw_on_kernel_error(int implementation_code,
                                     const std::string& file,
                                     const std::string& func,
                                     int line);
__host__ void _check_for_cuda_error(int implementation_code,
                                    const std::string& file,
                                    const std::string& func,
                                    int line);
__host__ void _check_thread_sync(int implementation_code, const std::string& file, const std::string& func, int line);
__host__ void _throw_on_global_error(const std::string& file, const std::string& func, int line);
__host__ void notify_of_global_error(int error_code);

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
                                               Complex* weighted_sfr_filtered_device,
                                               Complex* N_rec_filtered_device);

// Wrap exception-handling calls with these macros to add exception location information to the error messages
// N.B.: The ',' in the execution configuration part of a CUDA kernel call confuses the pre-processor ... so
//       make sure to add ()'s around the kernel call when using throw_on_kernel_error()
#define throw_on_generic_error(check_failure, implementation_code)                                                     \
  {                                                                                                                    \
    _throw_on_generic_error(check_failure, implementation_code, __FILE__, __func__, __LINE__);                         \
  }
#define throw_on_cuda_error(cuda_code, implementation_code)                                                            \
  {                                                                                                                    \
    _throw_on_cuda_error((cuda_code), implementation_code, __FILE__, __func__, __LINE__);                              \
  }
#define throw_on_kernel_error(kernel_call, implementation_code)                                                        \
  {                                                                                                                    \
    (kernel_call);                                                                                                     \
    _check_for_cuda_error(implementation_code, __FILE__, __func__, __LINE__);                                          \
  }
#define check_thread_sync(implementation_code)                                                                         \
  {                                                                                                                    \
    _check_thread_sync(implementation_code, __FILE__, __func__, __LINE__);                                             \
  }
#define throw_on_global_error()                                                                                        \
  {                                                                                                                    \
    _throw_on_global_error(__FILE__, __func__, __LINE__);                                                              \
  }
#endif

// Define base exception classes
#include <sstream>
#include <string>
#define GENERIC_CUDA_ERROR_CODE 160614
class meraxes_exception_base : public std::exception
{
protected:
  int _error_code;
  std::string _message;
  std::string _file;
  std::string _func;
  int _line;
  std::string _composition;

public:
  // Constructor (C++ STL strings).
  explicit meraxes_exception_base(int code,
                                  const std::string& message,
                                  const std::string& file,
                                  const std::string& func,
                                  const int line)
    : _error_code(code)
    , _message(message)
    , _file(file)
    , _func(func)
    , _line(line)
  {
    // Create the error message ...
    std::stringstream s;
    // ... add the error description to the exception message
    s << _message << " (code=";
    // ... add the error location to the exception message
    s << _error_code << ")";
    // Convert stream to string
    this->_composition = s.str();
  }

  // Destructor.  Virtual to allow for subclassing.
  virtual ~meraxes_exception_base() noexcept {}

  // The following functions return pointers.
  //    The underlying memory is possessed by the
  //    meraxes_exception object. Callers must
  //    not attempt to free the memory.
  virtual const char* what() const noexcept { return (_composition.c_str()); }
  virtual const char* file() const noexcept { return (_file.c_str()); }
  virtual const char* func() const noexcept { return (_func.c_str()); }

  // Returns an integer expressing the error code.
  virtual int error_code() const noexcept { return (this->_error_code); }

  // Returns an integer expressing the line number of the error.
  virtual int line() const noexcept { return (this->_line); }

  // Call this method inside catch blocks to process the exception
  // THIS IS THE CODE TO MODIFY IF YOU WANT TO ADJUST THE WAY
  //    MERAXES RESPONDS TO GPU ERRORS
  virtual void process_exception() const
  {
    // This loop aids in making the error reporting a bit cleaner in the log
    for (int i_rank = 0; i_rank < run_globals.mpi_size; i_rank++) {
      if (i_rank == run_globals.mpi_rank) {
        if (this->error_code() != 0)
          _mlog_error(this->what(), this->file(), this->func(), this->line());
      }
      MPI_Barrier(run_globals.mpi_comm);
    }
    if (this->error_code() != 0)
      ABORT(EXIT_FAILURE);
    else
      ABORT(EXIT_SUCCESS);
  }
};

// Derived meraxes exception class for CUDA exceptions
//    Define all implementation error codes and associated error messages here
class meraxes_cuda_exception : public meraxes_exception_base
{
public:
  // List all the implementation exception codes here
  enum _code_list
  {
    INVALID_FILTER,
    GLOBAL,
    INIT,
    INIT_PBS_GPUFILE,
    MALLOC,
    FREE,
    MEMCPY,
    SYNC,
    KERNEL_CMPLX_AX,
    KERNEL_SET_ARRAY,
    KERNEL_FILTER,
    KERNEL_CHECK,
    KERNEL_MAIN_LOOP
  };

private:
  std::string _set_message(int code)
  {
    // List the exception message information for each implementation exception code here
    switch (code) {
      case INVALID_FILTER:
        return ("Invalid filter specified in run_globals.");
      case GLOBAL:
        return ("Notified of error on differing rank.");
      case INIT:
        return ("CUDA error while initializing device.");
      case INIT_PBS_GPUFILE:
        return ("CUDA error while using PBS_GPUFILE to set device number");
      case MALLOC:
        return ("CUDA error while calling cudaMalloc()");
      case FREE:
        return ("CUDA error while calling cudaFree()");
      case MEMCPY:
        return ("CUDA error while calling cudaMemcpy()");
      case SYNC:
        return ("CUDA error while syncing device");
      // The following kernel error messages are meant to have
      //   modifiers placed in front of them to descriminate
      //   between CUDA errors and thread-sync errors.
      case KERNEL_CMPLX_AX:
        return ("scalar-times-complex-vector kernel execution");
      case KERNEL_SET_ARRAY:
        return ("initialize-vector-to-constant kernel execution");
      case KERNEL_FILTER:
        return ("filter/convolution kernel execution");
      case KERNEL_CHECK:
        return ("post-convolution-sanity-check kernel execution");
      case KERNEL_MAIN_LOOP:
        return ("main find_HII_bubbles() kernel execution");
      // This should never happen.
      default:
        return ("Undocumented CUDA error.");
    }
  }

public:
  // This is the constructor used for most standard exception handling
  explicit meraxes_cuda_exception(int cuda_code,
                                  int implementation_code,
                                  const std::string file,
                                  const std::string func,
                                  const int line)
    : meraxes_exception_base((cuda_code), _set_message(implementation_code), file, func, line)
  {
    if (implementation_code != GLOBAL)
      notify_of_global_error(cuda_code);
  }
  // This constructor deals with the case when we want to modify the _set_message() string.  This is
  //    used for specifying whether kernel errors are CUDA errors or thread-sync errors, for example.
  explicit meraxes_cuda_exception(int cuda_code,
                                  int implementation_code,
                                  const std::string& modifier,
                                  const std::string file,
                                  const std::string func,
                                  const int line)
    : meraxes_exception_base((cuda_code), modifier + _set_message(implementation_code), file, func, line)
  {
    if (implementation_code != GLOBAL)
      notify_of_global_error(cuda_code);
  }
};

#endif
#endif
