#ifndef _MERAXES_GPU_H

double RtoM(double R);

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
    fftwf_complex *deltax_filtered,  // complex
    fftwf_complex *stars_filtered,  // complex
    fftwf_complex *sfr_filtered,  // complex

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
        fftwf_complex *deltax_filtered,  // complex
        fftwf_complex *stars_filtered,  // complex
        fftwf_complex *sfr_filtered,  // complex

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
        ));

#define _MERAXES_GPU_H
#endif
