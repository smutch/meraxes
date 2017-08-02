#include "meraxes.h"
#include "meraxes_gpu.h"
#include "utils.h"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <limits.h>

#include <hdf5.h>
#include <hdf5_hl.h>

void find_HII_bubbles_driver(
    double redshift,
    void  (*_find_HII_bubbles_passed)(
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
        Complex *deltax_filtered_in,  // complex
        Complex *stars_filtered_in,  // complex
        Complex *sfr_filtered_in,  // complex
    
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
    timer_info *timer)
{
    int mpi_rank=run_globals.mpi_rank;
    int n_rank  =run_globals.mpi_size;
    if(n_rank!=1){
        if(mpi_rank==0) fprintf(stderr,"n_rank=%d but only n_rank==1 supported at this point.\n",n_rank);
        exit(1);
    }
    
    // Open inputs file
    char fname[STRLEN];
    sprintf(fname, "%s/validation_input-core%03d-z%.2f.h5",reference_directory,mpi_rank, redshift);
    hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    // Read input attributes
    double box_size;
    int    ReionGridDim;
    int    local_nix;
    int    flag_ReionUVBFlag;
    double ReionEfficiency;
    double ReionNionPhotPerBary;
    double UnitLength_in_cm;
    double UnitMass_in_g;
    double UnitTime_in_s;
    double ReionRBubbleMax;
    double ReionRBubbleMin;
    double ReionDeltaRFactor;
    double ReionGammaHaloBias;
    double ReionAlphaUV;
    double ReionEscapeFrac;
    H5LTget_attribute_double(file_id, "/", "redshift", &redshift);
    H5LTget_attribute_int   (file_id, "/", "mpi_rank", &mpi_rank);
    H5LTget_attribute_double(file_id, "/", "box_size", &box_size);
    H5LTget_attribute_int   (file_id, "/", "ReionGridDim", &ReionGridDim);
    H5LTget_attribute_int   (file_id, "/", "local_nix", &local_nix);
    H5LTget_attribute_int   (file_id, "/", "flag_ReionUVBFlag", &flag_ReionUVBFlag);
    H5LTget_attribute_double(file_id, "/", "ReionEfficiency", &ReionEfficiency);
    H5LTget_attribute_double(file_id, "/", "ReionNionPhotPerBary", &ReionNionPhotPerBary);
    H5LTget_attribute_double(file_id, "/", "UnitLength_in_cm", &UnitLength_in_cm);
    H5LTget_attribute_double(file_id, "/", "UnitMass_in_g", &UnitMass_in_g);
    H5LTget_attribute_double(file_id, "/", "UnitTime_in_s", &UnitTime_in_s);
    H5LTget_attribute_double(file_id, "/", "ReionRBubbleMax", &ReionRBubbleMax);
    H5LTget_attribute_double(file_id, "/", "ReionRBubbleMin", &ReionRBubbleMin);
    H5LTget_attribute_double(file_id, "/", "ReionDeltaRFactor", &ReionDeltaRFactor);
    H5LTget_attribute_double(file_id, "/", "ReionGammaHaloBias", &ReionGammaHaloBias);
    H5LTget_attribute_double(file_id, "/", "ReionAlphaUV", &ReionAlphaUV);
    H5LTget_attribute_double(file_id, "/", "ReionEscapeFrac", &ReionEscapeFrac);

    // Initialize fftw here (we need to do this because we need 'slab_n_complex' to set the size of the input datasets
    ptrdiff_t local_nix_2;
    ptrdiff_t local_ix_start;
    int       local_n_complex = fftwf_mpi_local_size_3d(ReionGridDim, ReionGridDim, ReionGridDim / 2 + 1, run_globals.mpi_comm, &local_nix_2, &local_ix_start);
    if(local_nix!=local_nix_2){
        printf("Error: local_nix!=local_nix_2 (ie. %d!=%d)\n",local_nix,local_nix_2);
        exit(1);
    }
    int       slab_n_complex    =  local_n_complex;
    int       slab_n_real       =  local_nix*ReionGridDim*ReionGridDim;
    ptrdiff_t slabs_n_complex[] = {slab_n_complex}; // works for 1 core only
    ptrdiff_t slabs_ix_start[]  = {0};              // works for 1 core only

    // Read input datasets
    float *deltax             = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    float *stars              = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    float *sfr                = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    float *z_at_ionization    = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    float *J_21_at_ionization = (float *)malloc(sizeof(float)*(slab_n_complex * 2));
    H5LTread_dataset_float(file_id, "deltax", deltax);
    H5LTread_dataset_float(file_id, "sfr",  sfr);
    H5LTread_dataset_float(file_id, "stars", stars);
    H5LTread_dataset_float(file_id, "z_at_ionization",    z_at_ionization);
    H5LTread_dataset_float(file_id, "J_21_at_ionization", J_21_at_ionization);
    H5Fclose(file_id);

    // Initialize outputs    
    float   *J_21                     =(float   *)malloc(sizeof(float)  *slab_n_real);
    float   *r_bubble                 =(float   *)malloc(sizeof(float)  *slab_n_real); 
    Complex *deltax_filtered          =(Complex *)malloc(sizeof(Complex)*slab_n_complex);
    Complex *stars_filtered           =(Complex *)malloc(sizeof(Complex)*slab_n_complex);  
    Complex *sfr_filtered             =(Complex *)malloc(sizeof(Complex)*slab_n_complex);
    float   *xH                       =(float   *)malloc(sizeof(float)  *slab_n_real);
    double   volume_weighted_global_xH;
    double   mass_weighted_global_xH;
    
    // Call the version of find_HII_bubbles we've been passed
    timer_start(timer);
    _find_HII_bubbles_passed(
        redshift, //
        run_globals.mpi_comm,
        mpi_rank,
        box_size, //
        ReionGridDim, //
        local_nix, //
        flag_ReionUVBFlag, //
        ReionEfficiency, //
        ReionNionPhotPerBary, //
        UnitLength_in_cm, //
        UnitMass_in_g, //
        UnitTime_in_s, //
        ReionRBubbleMax, //
        ReionRBubbleMin, //
        ReionDeltaRFactor, //
        ReionGammaHaloBias, //
        ReionAlphaUV, //
        ReionEscapeFrac, //
        true, //
        J_21,  
        r_bubble, 
        deltax, //
        stars, //
        sfr, //
        deltax_filtered,  
        stars_filtered,  
        sfr_filtered,  
        slabs_n_complex,
        slabs_ix_start,
        xH, 
        z_at_ionization,
        J_21_at_ionization,
        &volume_weighted_global_xH,
        &mass_weighted_global_xH);
    timer_stop(timer);
    
    // Clean-up
    free(J_21);
    free(r_bubble);
    free(deltax);
    free(stars);
    free(sfr);
    free(deltax_filtered);
    free(stars_filtered);
    free(sfr_filtered);
    free(xH);
    free(z_at_ionization);
    free(J_21_at_ionization);
}

