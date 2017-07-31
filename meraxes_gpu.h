#ifndef _MERAXES_GPU_H
#define _MERAXES_GPU_H
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "utils.h"
#include "meraxes.h"

#ifdef __NVCC__
typedef float2 Complex;
#else
typedef fftwf_complex Complex;
#endif
//static __device__ __host__ Complex ComplexAdd    (Complex, Complex);
//static __device__ __host__ Complex ComplexScale  (Complex, float);
//static __device__ __host__ Complex ComplexMul    (Complex, Complex);
//static __device__ __host__ Complex ComplexScalMul(Complex, double);

__global__ void complex_vector_times_scalar(Complex *vector,double scalar,int n);

#ifdef __cplusplus
extern "C" {
#endif
void  _find_HII_bubbles_gpu(
    // input
    double redshift,
    MPI_Comm mpi_comm,
    int mpi_rank,
    double box_size,
    int ReionGridDim,
    int local_nix,
    int flag_ReionUVBFlag,
    double ReionEfficiency,
    double ReionNionPhotPerBary,
    double UnitLength_in_cm,
    double UnitMass_in_g,
    double UnitTime_in_s,
    double ReionRBubbleMax,
    double ReionRBubbleMin,
    double ReionDeltaRFactor,
    double ReionGammaHaloBias,
    double ReionAlphaUV,
    double ReionEscapeFrac,

    bool validation_output,

    // preallocated 1D grids (local_nix * ReionGridDim * ReionGridDim)
    float *J_21,  // real
    float *r_bubble, // real

    // input grids
    float *deltax,  // real & padded
    float *stars,  // real & padded
    float *sfr,  // real & padded

    // preallocated
    Complex *deltax_filtered,  // complex
    Complex *stars_filtered,  // complex
    Complex *sfr_filtered,  // complex

    // length = mpi.size
    ptrdiff_t *slabs_n_complex,
    ptrdiff_t *slabs_ix_start,

    // output - preallocated real grids (local_nix * ReionGridDim * ReionGridDim)
    float *xH, // real
    float *z_at_ionization,
    float *J_21_at_ionization,

    // output - single values
    double *volume_weighted_global_xH,
    double *mass_weighted_global_xH
    );

void find_HII_bubbles_driver(
    double redshift,
    void  (*find_HII_bubbles_passed)(
        // input
        double redshift,
        MPI_Comm mpi_comm,
        int mpi_rank,
        double box_size,
        int ReionGridDim,
        int local_nix,
        int flag_ReionUVBFlag,
        double ReionEfficiency,
        double ReionNionPhotPerBary,
        double UnitLength_in_cm,
        double UnitMass_in_g,
        double UnitTime_in_s,
        double ReionRBubbleMax,
        double ReionRBubbleMin,
        double ReionDeltaRFactor,
        double ReionGammaHaloBias,
        double ReionAlphaUV,
        double ReionEscapeFrac,

        bool validation_output,

        // preallocated 1D grids (local_nix * ReionGridDim * ReionGridDim)
        float *J_21,  // real
        float *r_bubble, // real

        // input grids
        float *deltax,  // real & padded
        float *stars,  // real & padded
        float *sfr,  // real & padded

        // preallocated
        Complex *deltax_filtered,  // complex
        Complex *stars_filtered,  // complex
        Complex *sfr_filtered,  // complex

        // length = mpi.size
        ptrdiff_t *slabs_n_complex,
        ptrdiff_t *slabs_ix_start,

        // output - preallocated real grids (local_nix * ReionGridDim * ReionGridDim)
        float *xH, // real
        float *z_at_ionization,
        float *J_21_at_ionization,

        // output - single values
        double *volume_weighted_global_xH,
        double *mass_weighted_global_xH
        ),
    const char *reference_directory,
    timer_info *timer);
#ifdef __cplusplus
}
#endif
#define _MERAXES_GPU_H
#endif
