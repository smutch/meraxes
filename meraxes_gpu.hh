#ifndef _MERAXES_GPU_HH
#ifdef __cplusplus
#define _MERAXES_GPU_HH

#include <iostream>
#include <sstream>
#include <string>

// Wrap all cuda calls with these macros inside try blocks
#define throw_on_cuda_error(cuda_code,implementation_code)   { _throw_on_cuda_error ((cuda_code), implementation_code, __FILE__, __LINE__); }
#define throw_on_cuFFT_error(cufft_code,implementation_code) { _throw_on_cuFFT_error((cufft_code),implementation_code, __FILE__, __LINE__); }

// Define exception classes
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
        //    The underlying memory is in posession of the
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
            default:
                return("Undocumented CUDA error.");
            }
        }
    public:
        explicit meraxes_cuda_exception(int cuda_code,int implementation_code,const std::string file,const int line) : 
                 meraxes_exception_base((cuda_code),_set_message(implementation_code),file,line) {
        }
};

#endif
#endif
