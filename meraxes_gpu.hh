#ifndef _MERAXES_GPU_H_H
#ifdef __cplusplus
#define _MERAXES_GPU_H_H

#include <iostream>
#include <sstream>
#include <string>

// Stuff for exception handling
#define ERROR_CUDA_MALLOC 1
#define ERROR_CUDA_MEMCPY_H2D 2
#define ERROR_CUFFT_CREATE_PLAN 3
#define ERROR_CUFFT_SET_COMPATIBILITY 4
#define ERROR_CUFFT_R2C 5
#define ERROR_CUFFT_C2R 6
#define ERROR_CUFFT_PLAN_DESTROY 7
#define ERROR_CUDA_SYNC 8
class local_exception : public std::exception {
    protected:
        int         _error_code;
        std::string _message;

    public:

        // Constructor (C strings).
        //     The string contents are copied upon construction.
        //     Hence, responsibility for deleting the char* lies
        //     with the caller.
        explicit local_exception(int code,const char* message):
          _error_code(code),
          _message(message)
          {
          }

        // Constructor (C++ STL strings).
        explicit local_exception(int code,const std::string& message):
          _error_code(code),
          _message(message)
          {}

        // Destructor.
        // Virtual to allow for subclassing.
        virtual ~local_exception() noexcept {}

        // Returns a pointer to the (constant) error description.
        //    The underlying memory is in posession of the
        //    local_exception object. Callers must
        //    not attempt to free the memory.
        virtual const char* what() const noexcept {
           return this->_message.c_str();
        }

        // Returns an integer expressing the error code.
        virtual int error_code() const noexcept {
           return this->_error_code;
        }

        // Compose an exception message.
        virtual const std::string compose_message() const {
            std::ostringstream s;
            s << "code=" << this->error_code() << "; " << this->what();
            return s.str();
        }
};

#define throw_on_cuda_error(cuda_code,implementation_code)   { _throw_on_cuda_error ((cuda_code), implementation_code, __FILE__, __LINE__); }
#define throw_on_cuFFT_error(cufft_code,implementation_code) { _throw_on_cuFFT_error((cufft_code),implementation_code, __FILE__, __LINE__); }

#endif
#endif
