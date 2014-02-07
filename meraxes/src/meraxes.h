#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <hdf5.h>

#ifndef _INIT_MERAXES
#define _INIT_MERAXES

#ifdef USE_TOCF
#include <21cmfast.h>
#endif

/*
 * Definitions
 */

#define STRLEN  256  //!< Default string length

// TODO: This should not be hard coded if at all possible...
#ifndef MAXSNAPS
#define MAXSNAPS 235  //!< Maximum number of snapshots
#endif

#ifndef NOUT
#define NOUT 1
#endif

#ifndef MAX_PHOTO_NBANDS
#define MAX_PHOTO_NBANDS 5
#endif
#ifndef N_PHOTO_JUMPS
#define N_PHOTO_JUMPS 1000
#endif

#define MVIR_PROP 1
#define VMAX_PROP 2

#define ABORT(sigterm)                                                                 \
do {                                                                                   \
  SID_log_error("in file: %s\tfunc: %s\tline: %i", __FILE__, __FUNCTION__, __LINE__);  \
  myexit(sigterm);                                                                     \
} while(0)

// Units (cgs):
#define GRAVITY          6.672e-8
#define SOLAR_MASS       1.989e33
#define SOLAR_LUM        3.826e33
#define RAD_CONST        7.565e-15
#define AVOGADRO         6.0222e23
#define BOLTZMANN        1.3806e-16
#define GAS_CONST        8.31425e7
#define C                2.9979e10
#define PLANCK           6.6262e-27
#define PROTONMASS       1.6726e-24
#define HUBBLE           3.2407789e-18 //! [h/sec]
#define SEC_PER_MEGAYEAR 3.155e13
#define SEC_PER_YEAR     3.155e7


/*
 * Structures
 */

//! Physics parameter values
struct physics_params_t{
  int    funcprop;
  double peak;
  double sigma;
  double stellarfrac;
  double peak_evo;
  double sigma_evo;
  double stellarfrac_evo;
  double bhgrowthfactor;
  double reion_z_re;
  double reion_delta_z_re;
  double reion_delta_z_sc;
  double reion_T0;
  double reion_Tcool;
  // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value
  double reion_Nion_phot_per_bary;
  double reion_escape_frac;
  double reion_mean_n_rec; // Mean number of recombinations per baryon
};
typedef struct physics_params_t physics_params_t;

//! Run params
//! Everything in this structure is supplied by the user...
struct run_params_t{
  char                  filename[STRLEN];
  char                  OutputDir[STRLEN];
  char                  FileNameGalaxies[STRLEN];
  char                  SimName[STRLEN];
  char                  SimulationDir[STRLEN];
  char                  CatalogFilePrefix[STRLEN];
  char                  FileWithOutputSnaps[STRLEN];
  char                  PhotometricTablesDir[STRLEN];
  char                  SSPModel[STRLEN];
  char                  IMF[STRLEN];
  char                  MagSystem[STRLEN];
  char                  MagBands[STRLEN];
  int                   NEverySnap;
  int                   NScanSnap;
  int                   TotalSimSnaps;
  int                   LastSnapshotNr;
  int                   FirstFile;
  int                   LastFile;
  double                BoxSize;
  double                VolumeFactor;
  double                ThreshMajorMerger;
  double                RecycleFraction;
  double                Hubble_h;
  double                BaryonFrac;
  double                OmegaM;
  double                OmegaK;
  double                OmegaR;
  double                OmegaLambda;
  double                Sigma8;
  double                wLambda;
  double                SpectralIndex;
  double                PartMass;
  double                MergerTimeFactor;
  int                   SnaplistLength;
  int                   RandomSeed;
  physics_params_t      physics;
  int                   TOCF_Flag;
};
typedef struct run_params_t run_params_t;

struct run_units_t{
  double UnitTime_in_s;
  double UnitLength_in_cm;
  double UnitVelocity_in_cm_per_s;
  double UnitTime_in_Megayears;
  double UnitMass_in_g;
  double UnitDensity_in_cgs;
  double UnitPressure_in_cgs;
  double UnitCoolingRate_in_cgs;
  double UnitEnergy_in_cgs;
};
typedef struct run_units_t run_units_t;

struct hdf5_output_t
{
  size_t         dst_size;
  hid_t          array3f_tid;
  hid_t          array_nmag_f_tid;
  size_t        *dst_offsets;
  size_t        *dst_field_sizes;
  const char   **field_names;
  hid_t         *field_types;
  int            n_props;
};
typedef struct hdf5_output_t hdf5_output_t;

struct phototabs_t{
  int   JumpTable[N_PHOTO_JUMPS];
  int   NAges;
  int   NBands;
  int   NMetals;
  float JumpFactor;
  float *Table;
  float *Ages;
  float *Metals;
  char  (*MagBands)[5];
};
typedef struct phototabs_t phototabs_t;

#ifdef USE_TOCF
struct tocf_grids_t
{
  float *xH;
  fftwf_complex *stars;
  fftwf_complex *stars_filtered;
  fftwf_complex *stars_copy;
  fftwf_complex *deltax;
  fftwf_complex *deltax_filtered;
  fftwf_complex *deltax_copy;
  fftwf_complex *sfr;
  fftwf_complex *sfr_filtered;
  fftwf_complex *sfr_copy;
  fftwf_complex *N_rec;
  fftwf_complex *N_rec_filtered;
  float *z_at_ionization;
  float *J_21_at_ionization;
  float *J_21;
  float *Mvir_crit;
  float *mfp;
  float  global_xH;
};
typedef struct tocf_grids_t tocf_grids_t;
#endif

//! Global variables which will will be passed around
struct run_globals_t{
  int                        LastOutputSnap;
  int                        ListOutputSnaps[NOUT];
  int                        NGhosts;
  double                     AA[MAXSNAPS];
  double                     ZZ[MAXSNAPS];
  double                     LTTime[MAXSNAPS];
  double                     Hubble;
  double                     RhoCrit;
  double                     G;
  char                       FNameOut[STRLEN];
  struct galaxy_t      *FirstGal;
  struct galaxy_t      *LastGal;
  gsl_rng                   *random_generator;
  struct run_params_t   params;
  struct run_units_t    units;
  hdf5_output_t         hdf5props;
  phototabs_t           photo;
#ifdef USE_TOCF
  tocf_grids_t          tocf_grids;
#endif
};
typedef struct run_globals_t run_globals_t;

//! Tree info struct
typedef struct trees_info_t{
  int n_step;
  int n_search;
  int n_halos;
  int n_halos_max;
  int max_tree_id;
  int n_fof_groups;
} trees_info_t;

//! Tree entry struct
typedef struct tree_entry_t{
  int id;
  int flags;
  int desc_id;
  int tree_id;
  int file_offset;
  int desc_index;
  int central_index;
  int forest_id;
  double fof_mvir;
} tree_entry_t;

//! This is the structure for a halo in the catalog files
struct catalog_halo_t{
  long long id_MBP;                    //!< ID of most bound particle in structure
  double    M_vir;                     //!< Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int       n_particles;               //!< Number of particles in the structure
  float     position_COM[3];           //!< Centre-of-mass position      [Mpc/h]
  float     position_MBP[3];           //!< Most bound particle position [Mpc/h]
  float     velocity_COM[3];           //!< Centre-of-mass velocity      [km/s]
  float     velocity_MBP[3];           //!< Most bound particle velocity [km/s]
  float     R_vir;                     //!< Virial radius [Mpc/h]
  float     R_halo;                    //!< Distance of last halo particle from MBP [Mpc/h]
  float     R_max;                     //!< Radius of maximum circular velocity     [Mpc/h]
  float     V_max;                     //!< Maximum circular velocity               [km/s]
  float     sigma_v;                   //!< Total 3D velocity dispersion            [km/s]
  float     spin[3];                   //!< Specific angular momentum vector        [Mpc/h*km/s]
  float     q_triaxial;                //!< Triaxial shape parameter q=b/a
  float     s_triaxial;                //!< Triaxial shape parameter s=c/a
  float     shape_eigen_vectors[3][3]; //!< Normalized triaxial shape eigenvectors
  char      padding[8];                //!< Alignment padding
};
typedef struct catalog_halo_t catalog_halo_t;


//! The meraxis halo structure
struct halo_t{
  long long id_MBP;      //!< ID of most bound particle
  int    ID;             //!< Halo ID
  int    Type;           //!< Type (0 for central, 1 for satellite)
  int    SnapOffset;     //!< Number of snapshots this halo skips before reappearing
  int    DescIndex;      //!< Index of descendant in next relevant snapshot
  int    TreeFlags;      //!< Bitwise flag indicating the type of match in the trees
  struct fof_group_t *FOFGroup;
  struct halo_t      *NextHaloInFOFGroup;
  struct galaxy_t    *Galaxy;
  double Mvir;           //!< virial mass [M_sol/h]
  int    Len;            //!< Number of particles in the structure
  float  Pos[3];         //!< Most bound particle position [Mpc/h]
  float  Vel[3];         //!< Centre-of-mass velocity [km/s]
  float  Rvir;           //!< Virial radius [Mpc/h]
  float  Rhalo;          //!< Distance of last halo particle from MBP [Mpc/h]
  float  Rmax;           //!< Radius of maximum circular velocity [Mpc/h]
  float  Vvir;           //!< Virial velocity [km/s]
  float  Vmax;           //!< Maximum circular velocity [km/s]
  float  VelDisp;        //!< Total 3D velocity dispersion [km/s]
  float  Spin[3];        //!< Specific angular momentum vector [Mpc/h *km/s]
};
typedef struct halo_t halo_t;

struct fof_group_t{
  halo_t *FirstHalo;
};
typedef struct fof_group_t fof_group_t;

struct galaxy_t
{
  long long id_MBP;
  int    ID;
  int    Type;
  int    OldType;
  int    SnapSkipCounter;
  int    HaloDescIndex;
  int    TreeFlags;
  bool   ghost_flag;
  struct halo_t         *Halo;
  struct galaxy_t       *FirstGalInHalo;
  struct galaxy_t       *NextGalInHalo;
  struct galaxy_t       *Next;
  struct galaxy_t       *MergerTarget;
  int    Len;
  double dt;      //!< Time between current snapshot and last identification
  double LTTime;  //!< Lookback time at the last time this galaxy was identified

  // properties of subhalo at the last time this galaxy was a central galaxy
  double Pos[3];
  double Vel[3];
  double Mvir;
  double dM;
  double Rvir;
  double Vvir;
  double Vmax;

  // baryonic reservoirs
  double Gas;
  double StellarMass;
  double Sfr;

  // misc
  double Cos_Inc;
  double MergTime;

  // write index
  int output_index;

#ifdef CALC_MAGS
  double Lum[MAX_PHOTO_NBANDS][NOUT];
#endif

};
typedef struct galaxy_t galaxy_t;

struct galaxy_output_t
{
  long long id_MBP;
  int   ID;
  int   Type;
  int   CentralGal;
  int   GhostFlag;

  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;
  float Mvir;
  float dM;
  float dMdt;
  float Rvir;
  float Vvir;
  float Vmax;

  // baryonic reservoirs
  float StellarMass;
  float Sfr;

  // misc
  float Cos_Inc;
  float MergTime;
  float LTTime;

#ifdef CALC_MAGS
  float Mag[MAX_PHOTO_NBANDS];
  float MagDust[MAX_PHOTO_NBANDS];
#endif
};
typedef struct galaxy_output_t galaxy_output_t;


/*
 * Functions
 */

void    myexit(int signum);
void    read_parameter_file(run_globals_t *run_globals, char *fname);
void    init_meraxes(run_globals_t *run_globals);
void    dracarys(run_globals_t *run_globals);
int     evolve_galaxies(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot, int NGal, int NFof);
trees_info_t read_halos(run_globals_t *run_globals, int snapshot, halo_t **halo, fof_group_t **fof_group);
void    free_halos(halo_t **halo);
galaxy_t* new_galaxy(run_globals_t *run_globals, int *unique_ID);
void    copy_halo_to_galaxy(halo_t *halo, galaxy_t *gal, int snapshot);
void    gas_infall(run_globals_t *run_globals, fof_group_t *FOFgroup, int snapshot);
double  calculate_merging_time(run_globals_t *run_globals, galaxy_t *gal, int snapshot);
void    merge_with_target(run_globals_t *run_globals, galaxy_t *gal, int *dead_gals);
void    form_stars_insitu(run_globals_t *run_globals, galaxy_t *gal, int snapshot);
void    prep_hdf5_file(run_globals_t *run_globals);
void    write_snapshot(run_globals_t *run_globals, int n_write, int i_out, int *last_n_write);
void    calc_hdf5_props(run_globals_t *run_globals);
void    prepare_galaxy_for_output(run_globals_t *run_globals, galaxy_t gal, galaxy_output_t *galout, int i_snap);
void    read_photometric_tables(run_globals_t *run_globals);
void    mpi_debug_here();
void    check_counts(run_globals_t *run_globals, fof_group_t *fof_group, int NGal, int NFof);
void    cn_quote();
int     get_corrected_snapshot(run_globals_t *run_globals, int snapshot);
double  Tvir_to_Mvir(run_globals_t *run_globals, double T, double z);
double  calculate_Mvir(run_globals_t *run_globals, halo_t *halo);
float   calculate_Rvir(run_globals_t *run_globals, halo_t *halo, double Mvir, int snapshot);
float   calculate_Vvir(run_globals_t *run_globals, double Mvir, float Rvir);

// Magnitude related
void    init_luminosities(run_globals_t *run_globals, galaxy_t *gal);
void    add_to_luminosities(run_globals_t *run_globals, galaxy_t *gal, double burst_mass, double metallicity, double burst_time);
double  lum_to_mag(double lum);
void    sum_luminosities(run_globals_t *run_globals, galaxy_t *parent, galaxy_t *gal, int outputbin);
void    prepare_magnitudes_for_output(run_globals_t *run_globals, galaxy_t gal, galaxy_output_t *galout, int i_snap);
void    apply_dust(int n_photo_bands, galaxy_t gal, double *LumDust, int outputbin);
void    cleanup_mags(run_globals_t *run_globals);

// Reionization related
// bool    check_reionization_cooling(run_globals_t *run_globals, halo_t *halo, int snapshot);
double  reionization_modifier(run_globals_t *run_globals, halo_t *halo, int snapshot);
double  global_ionizing_emmisivity(run_globals_t *run_globals);
#ifdef USE_TOCF
void    set_HII_eff_factor(run_globals_t *run_globals);
int     find_cell(double pos, double box_size);
void    malloc_reionization_grids(run_globals_t *run_globals);
void    free_reionization_grids(run_globals_t *run_globals);
void    construct_stellar_grids(run_globals_t *run_globals);
// void    assign_ionization_to_halos(run_globals_t *run_globals, halo_t *halo, int n_halos, float *xH_grid, int xH_dim);
int     read_dm_grid(run_globals_t *run_globals, int snapshot, int i_grid, float *grid);
void    calculate_Mvir_crit(run_globals_t *run_globals, double redshift);
void    call_find_HII_bubbles(run_globals_t *run_globals, int snapshot, int nout_gals);
void    save_tocf_grids(run_globals_t *run_globals, hid_t group_id, int snapshot);
void    check_if_reionization_complete(run_globals_t *run_globals);
#endif

#endif // _INIT_MERAXES
