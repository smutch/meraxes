#ifndef _MERAXES_GPU_H
#define _MERAXES_GPU_H

// Stuff used for testing timer
#include "utils.h"

#ifndef __cplusplus
typedef fftwf_complex Complex;
#else
#include <cuda_runtime.h>
typedef float2 Complex;
extern "C" {
#endif
void _find_HII_bubbles_gpu(double redshift,const bool flag_write_validation_output);
void init_CUDA();

#include "utils.h"
void find_HII_bubbles_driver(
    int        snapshot,
    const char *reference_directory,
    const bool flag_write_validation_data,
    timer_info *timer);

#ifdef __cplusplus
    }
#endif

// All the c++11-specific stuff 
#include "meraxes_gpu.hh"

#define _MERAXES_GPU_H
#endif
