#ifndef _MERAXES_GPU_HH
#ifdef __cplusplus
#define _MERAXES_GPU_HH

#ifdef __NVCC__
#include <cuda_runtime.h>
#include <cufft.h>
// These routines are needed by the kernel
__device__ void  inline grid_index2indices(const int idx,const int dim,const int local_ix_start,const int mode,int *ix,int *iy,int *iz);
__device__ int   inline grid_index_gpu(int i, int j, int k, int dim, int type);
__device__ float inline k_mag_gpu(const int n_x,const int n_y,const int n_z,const int dim,const float box_size);

// Kernel definitions
__global__ void set_array_gpu(float *array,int n_real,float val);
__global__ void sanity_check_aliasing(Complex *grid,int grid_dim,int local_ix_start,int n_real,float val);
__global__ void complex_vector_times_scalar(Complex *vector,double scalar,int n_complex);
__global__ void filter_gpu(Complex *box,int slab_nx,int grid_dim,int local_ix_start,int n_complex,float R,double box_size_in,int filter_type);
__global__ void find_HII_bubbles_gpu_main_loop(
                   float    redshift,
                   int      n_real,
                   int      flag_last_filter_step,
                   int      flag_ReionUVBFlag,
                   int      ReionGridDim,
                   int      local_ix_start,
                   float    R,
                   float    M,
                   float    ReionEfficiency,
                   float    inv_pixel_volume,
                   float    J_21_aux_constant,
                   double   ReionGammaHaloBias,
                   float   *xH,
                   float   *J_21,
                   float   *r_bubble,
                   float   *J_21_at_ionization,
                   float   *z_at_ionization,
                   Complex *deltax_filtered_device,
                   Complex *stars_filtered_device,
                   Complex *sfr_filtered_device);

// Host-side exception-handling routines
__host__ void  _throw_on_cuda_error  (cudaError_t cuda_code, int implementation_code, const std::string file, int line);
__host__ void  _throw_on_cuFFT_error (cufftResult cuda_code, int implementation_code, const std::string file, int line);
__host__ void  _throw_on_kernel_error(int implementation_code, const std::string file, int line);
__host__ void  _check_for_cuda_error (int implementation_code,const std::string file, int line);
__host__ void  _check_thread_sync    (int implementation_code,const std::string file, int line);

// Wrap exception-handling calls with these macros to add exception location information to the error messages
// N.B.: The ',' in the execution configuration part of a CUDA kernel call confuses the pre-processor ... so 
//       make sure to add ()'s around the kernel call when using throw_on_kernel_error()
#define throw_on_cuda_error(cuda_code,implementation_code)    { _throw_on_cuda_error ((cuda_code),implementation_code, __FILE__, __LINE__); }
#define throw_on_cuFFT_error(cufft_code,implementation_code)  { _throw_on_cuFFT_error((cufft_code),implementation_code, __FILE__, __LINE__); }
#define throw_on_kernel_error(kernel_call,implementation_code){ (kernel_call);_check_for_cuda_error(implementation_code, __FILE__, __LINE__); }
#define check_thread_sync(implementation_code)                { _check_thread_sync(implementation_code,__FILE__, __LINE__); }
#endif

// Define exception classes
#include <sstream>
#include <string>
class meraxes_exception_base : public std::exception {
    protected:
        int         _error_code;
        std::string _message;
        std::string _file;
        int         _line;
        std::string _composition;
    public:
        // Constructor (C++ STL strings).
        explicit meraxes_exception_base(int code,const std::string& message,const std::string& file,const int line):
          _error_code(code),
          _message(message),
          _file(file),
          _line(line)
          {
            // Create the error message ...
            std::stringstream s;
            // ... add the error description to the exception message
            s << _message << ": ";
            // ... add the error location to the exception message
            s << _file << "(" << std::to_string(_line) << ")";
            // Convert stream to string
            this->_composition = s.str();
          }

        // Destructor.  Virtual to allow for subclassing.
        virtual ~meraxes_exception_base() noexcept {}

        // Returns a pointer to the error description.
        //    The underlying memory is possessed by the
        //    meraxes_exception object. Callers must
        //    not attempt to free the memory.
        virtual const char* what() const noexcept {
            return(_composition.c_str());
        }

        // Returns an integer expressing the error code.
        virtual int error_code() const noexcept {
           return(this->_error_code);
        }

        // Call this method inside catch blocks to process the exception
        // THIS IS THE CODE TO MODIFY IF YOU WANT TO ADJUST THE WAY
        //    MERAXES RESPONDS TO GPU ERRORS
        virtual void process_exception() const{
            mlog(this->what(),MLOG_MESG);
            ABORT(this->error_code());
        }
};

// Derived meraxes exception class for CUDA exceptions
//    Defineall implementation error codes and associated error messages here
class meraxes_cuda_exception : public meraxes_exception_base {
    public:
        // List all the implementation exception codes here
        enum _code_list{
            MALLOC,
            FREE,
            MEMCPY,
            SYNC,
            CUFFT_CREATE_PLAN,
            CUFFT_SET_COMPATIBILITY,
            CUFFT_R2C,
            CUFFT_C2R,
            CUFFT_PLAN_DESTROY,
            KERNEL_CMPLX_AX,
            KERNEL_SET_ARRAY,
            KERNEL_FILTER,
            KERNEL_CHECK,
            KERNEL_MAIN_LOOP
            }; 
    private:
        std::string _set_message(int code){
            // List the exception message information for each implementation exception code here
            switch (code){
            case MALLOC:
                return("CUDA error while calling cudaMalloc()");
            case FREE:
                return("CUDA error while calling cudaFree()");
            case MEMCPY:
                return("CUDA error while calling cudaMemcpy()");
            case SYNC:
                return("CUDA error while syncing device");
            case CUFFT_CREATE_PLAN:
                return("cuFFT error while creating cuFFT plan");
            case CUFFT_SET_COMPATIBILITY:
                return("cuFFT error while setting cuFFT compatibility");
            case CUFFT_R2C:
                return("cuFFT error while performing real->complex cuFFT transform");
            case CUFFT_C2R:
                return("cuFFT error while performing complex->real cuFFT transform");
            case CUFFT_PLAN_DESTROY:
                return("cuFFT error while destroying cuFFT plan");
            // The following kernel error messages are meant to have
            //   modifiers placed in front of them to descriminate
            //   between CUDA errors and thread-sync errors.
            case KERNEL_CMPLX_AX:
                return("scalar-times-complex-vector kernel execution");
            case KERNEL_SET_ARRAY:
                return("initialize-vector-to-constant kernel execution");
            case KERNEL_FILTER:
                return("filter/convolution kernel execution");
            case KERNEL_CHECK:
                return("post-convolution-sanity-check kernel execution");
            case KERNEL_MAIN_LOOP:
                return("main find_HII_bubbles() kernel execution");
            // This should never happen.
            default:
                return("Undocumented CUDA error.");
            }
        }
    public:
        // This is the constructor used for most standard exception handling
        explicit meraxes_cuda_exception(int cuda_code,int implementation_code,const std::string file,const int line) : 
                 meraxes_exception_base((cuda_code),_set_message(implementation_code),file,line) {
        }
        // This constructor deals with the case when we want to modify the _set_message() string.  This is
        //    used for specifying whether kernel errors are CUDA errors or thread-sync errors, for example.
        explicit meraxes_cuda_exception(int cuda_code,int implementation_code,const std::string& modifier,const std::string file,const int line) : 
                 meraxes_exception_base((cuda_code),modifier+_set_message(implementation_code),file,line) {
        }
};

#endif
#endif
