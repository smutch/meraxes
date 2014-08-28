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

#ifndef NOUT
#define NOUT 1
#endif

#ifndef MAX_PHOTO_NBANDS
#define MAX_PHOTO_NBANDS 5
#endif
#ifndef N_PHOTO_JUMPS
#define N_PHOTO_JUMPS 1000
#endif

#ifndef N_HISTORY_SNAPS
#define N_HISTORY_SNAPS 5
#endif

#define ABORT(sigterm)                                                                 \
  do {                                                                                   \
    fprintf(stderr, "\nIn file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__);  \
    myexit(sigterm);                                                                     \
  } while (0)

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

#ifdef DEBUG
FILE *meraxes_debug_file;
#endif

/*
 * Structures
 */

//! Physics parameter values
struct physics_params_t {
  int Flag_ReionizationModifier;
  int Flag_BHFeedback;
  int Flag_SnDelay;

  double SfEfficiency;
  double SfRecycleFraction;
  double SnReheatEff;
  double SnEjectionEff;
  double ReincorporationEff;
  double Yield;
  double IMFSlope;
  double EnergyPerSN;
  double IMFNormConst;
  double RadioModeEff;
  double BlackHoleGrowthRate;

  double ThreshMajorMerger;
  double MinMergerRatioForBurst;
  double MergerBurstFactor;

  // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value
  double ReionNionPhotPerBary;
  double ReionEscapeFrac;
  double ReionMeanNRec; // Mean number of recombinations per baryon
  double ReionTcool;

  double ReionSobacchi_Zre;
  double ReionSobacchi_DeltaZre;
  double ReionSobacchi_DeltaZsc;
  double ReionSobacchi_T0;

  double ReionGnedin_z0;
  double ReionGnedin_zr;
};
typedef struct physics_params_t physics_params_t;

//! Run params
//! Everything in this structure is supplied by the user...
struct run_params_t {
  char             OutputDir[STRLEN];
  char             FileNameGalaxies[STRLEN];
  char             SimName[STRLEN];
  char             SimulationDir[STRLEN];
  char             CatalogFilePrefix[STRLEN];
  char             FileWithOutputSnaps[STRLEN];
  char             PhotometricTablesDir[STRLEN];
  char             CoolingFuncsDir[STRLEN];
  char             SSPModel[STRLEN];
  char             IMF[STRLEN];
  char             MagSystem[STRLEN];
  char             MagBands[STRLEN];
  int              FirstFile;
  int              LastFile;
  int              NSteps;
  double           BoxSize;
  double           VolumeFactor;
  double           Hubble_h;
  double           BaryonFrac;
  double           OmegaM;
  double           OmegaK;
  double           OmegaR;
  double           OmegaLambda;
  double           Sigma8;
  double           wLambda;
  double           SpectralIndex;
  double           PartMass;
  double           MergerTimeFactor;
  int              SnaplistLength;
  int              RandomSeed;
  char             ForestIDFile[STRLEN];
  physics_params_t physics;
  int              FlagInteractive;
  int              FlagGenDumpFile;
  int              FlagReadDumpFile;
  int              TOCF_Flag;
};
typedef struct run_params_t run_params_t;

struct run_units_t {
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

struct hdf5_output_t {
  size_t       dst_size;
  hid_t        array3f_tid;
  hid_t        array_nmag_f_tid;
  hid_t        array_nhist_f_tid;
  size_t      *dst_offsets;
  size_t      *dst_field_sizes;
  const char **field_names;
  hid_t       *field_types;
  int          n_props;
};
typedef struct hdf5_output_t hdf5_output_t;

struct phototabs_t {
  int    JumpTable[N_PHOTO_JUMPS];
  int    NAges;
  int    NBands;
  int    NMetals;
  float  JumpFactor;
  float *Table;
  float *Ages;
  float *Metals;
  char  (*MagBands)[5];
};
typedef struct phototabs_t phototabs_t;

#ifdef USE_TOCF
struct tocf_grids_t {
  float         *xH;
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
  float         *z_at_ionization;
  float         *J_21_at_ionization;
  float         *J_21;
  float         *Mvir_crit;
  float         *mfp;
  float          global_xH;
};
typedef struct tocf_grids_t tocf_grids_t;
#endif

//! The meraxis halo structure
typedef struct halo_t {
  long long           id_MBP; //!< ID of most bound particle
  int                 ID; //!< Halo ID
  int                 Type; //!< Type (0 for central, 1 for satellite)
  int                 SnapOffset; //!< Number of snapshots this halo skips before reappearing
  int                 DescIndex; //!< Index of descendant in next relevant snapshot
  int                 TreeFlags; //!< Bitwise flag indicating the type of match in the trees
  int                 ForestID;
  struct fof_group_t *FOFGroup;
  struct halo_t      *NextHaloInFOFGroup;
  struct galaxy_t    *Galaxy;
  double              Mvir; //!< virial mass [M_sol/h]
  int                 Len; //!< Number of particles in the structure
  float               Pos[3]; //!< Most bound particle position [Mpc/h]
  float               Vel[3]; //!< Centre-of-mass velocity [km/s]
  double              Rvir; //!< Virial radius [Mpc/h]
  float               Rhalo; //!< Distance of last halo particle from MBP [Mpc/h]
  float               Rmax; //!< Radius of maximum circular velocity [Mpc/h]
  double              Vvir; //!< Virial velocity [km/s]
  float               Vmax; //!< Maximum circular velocity [km/s]
  float               VelDisp; //!< Total 3D velocity dispersion [km/s]
  float               AngMom[3]; //!< Specific angular momentum vector [Mpc/h *km/s]
} halo_t;

struct fof_group_t {
  halo_t *FirstHalo;
  halo_t *FirstOccupiedHalo;
  double  Mvir;
  double  Rvir;
  double  Vvir;
  int     TotalSubhaloLen;
};
typedef struct fof_group_t fof_group_t;

struct galaxy_t {
  long long        id_MBP;
  int              ID;
  int              Type;
  int              OldType;
  int              SnapSkipCounter;
  int              HaloDescIndex;
  int              TreeFlags;
  int              LastIdentSnap;  //!< Last snapshot at which the halo in which this galaxy resides was identified
  bool             ghost_flag;
  struct halo_t   *Halo;
  struct galaxy_t *FirstGalInHalo;
  struct galaxy_t *NextGalInHalo;
  struct galaxy_t *Next;
  struct galaxy_t *MergerTarget;
  int              Len;
  double           dt; //!< Time between current snapshot and last identification

  // properties of subhalo at the last time this galaxy was a central galaxy
  float Pos[3];
  float Vel[3];
  double Mvir;
  double Rvir;
  double Vvir;
  double Vmax;
  double Spin;

  // baryonic reservoirs
  double HotGas;
  double MetalsHotGas;
  double ColdGas;
  double MetalsColdGas;
  double Mcool;
  double StellarMass;
  double MetalsStellarMass;
  double DiskScaleLength;
  double Sfr;
  double EjectedGas;
  double MetalsEjectedGas;
  double BlackHoleMass;

  // baryonic hostories
  double NewStars[N_HISTORY_SNAPS];

  // misc
  double Rcool;
  double Cos_Inc;
  double MergTime;
  double BaryonFracModifier;

  // write index
  int output_index;

#ifdef CALC_MAGS
  double Lum[MAX_PHOTO_NBANDS][NOUT];
#endif
};
typedef struct galaxy_t galaxy_t;

struct galaxy_output_t {
  long long id_MBP;
  int       ID;
  int       Type;
  int       CentralGal;
  int       GhostFlag;

  float Pos[3];
  float Vel[3];
  float Spin;
  int   Len;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float FOFMvir;

  // baryonic reservoirs
  float HotGas;
  float MetalsHotGas;
  float ColdGas;
  float MetalsColdGas;
  float Mcool;
  float DiskScaleLength;
  float StellarMass;
  float MetalsStellarMass;
  float Sfr;
  float EjectedGas;
  float MetalsEjectedGas;
  float BlackHoleMass;

  // baryonic histories
  float NewStars[N_HISTORY_SNAPS];

  // misc
  float Rcool;
  float Cos_Inc;
  float MergTime;
  float BaryonFracModifier;

#ifdef CALC_MAGS
  float Mag[MAX_PHOTO_NBANDS];
  float MagDust[MAX_PHOTO_NBANDS];
#endif
};
typedef struct galaxy_output_t galaxy_output_t;


//! Tree info struct
typedef struct trees_info_t {
  int n_step;
  int n_search;
  int n_halos;
  int n_halos_max;
  int max_tree_id;
  int n_fof_groups;
  int n_fof_groups_max;
  int unsampled_snapshot;
} trees_info_t;

//! Tree entry struct
typedef struct tree_entry_t {
  int    id;
  int    flags;
  int    desc_id;
  int    tree_id;
  int    file_offset;
  int    desc_index;
  int    central_index;
  int    forest_id;
  int    group_index;
} tree_entry_t;

//! This is the structure for a halo in the catalog files
typedef struct catalog_halo_t {
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
  float     ang_mom[3];                //!< Specific angular momentum vector        [Mpc/h*km/s]
  float     q_triaxial;                //!< Triaxial shape parameter q=b/a
  float     s_triaxial;                //!< Triaxial shape parameter s=c/a
  float     shape_eigen_vectors[3][3]; //!< Normalized triaxial shape eigenvectors
  char      padding[8];                //!< Alignment padding
} catalog_halo_t;

//! Global variables which will will be passed around
typedef struct run_globals_t {
  int                 LastOutputSnap;
  int                 ListOutputSnaps[NOUT];
  int                 NGhosts;
  double             *AA;
  double             *ZZ;
  double             *LTTime;
  double              Hubble;
  double              RhoCrit;
  double              G;
  char                FNameOut[STRLEN];
  int                 NHalosMax;
  int                 NFOFGroupsMax;
  int                 NRequestedForests;
  int                 TreesStep;
  int                 TreesScan;
  bool                SelectForestsSwitch;
  int                *RequestedForestId;
  int                 NStoreSnapshots;
  halo_t            **SnapshotHalo;
  fof_group_t       **SnapshotFOFGroup;
  int               **SnapshotIndexLookup;
  trees_info_t       *SnapshotTreesInfo;
  struct galaxy_t    *FirstGal;
  struct galaxy_t    *LastGal;
  gsl_rng            *random_generator;
  struct run_params_t params;
  struct run_units_t  units;
  hdf5_output_t       hdf5props;
  phototabs_t         photo;
#ifdef USE_TOCF
  tocf_grids_t tocf_grids;
#endif
} run_globals_t;


/*
 * Functions
 */

void         myexit(int signum);
void         read_parameter_file(run_globals_t *run_globals, char *fname, int mode);
void         init_meraxes(run_globals_t *run_globals);
void         continue_prompt(run_globals_t *run_globals, char *param_file);
void         free_halo_storage(run_globals_t *run_globals);
void         initialize_halo_storage(run_globals_t *run_globals);
void         dracarys(run_globals_t *run_globals);
int          evolve_galaxies(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot, int NGal, int NFof);
void         passively_evolve_ghost(run_globals_t *run_globals, galaxy_t *gal, int snapshot);
trees_info_t read_halos(run_globals_t *run_globals, int snapshot, halo_t **halo, fof_group_t **fof_group, int **index_lookup, trees_info_t *snapshot_trees_info);
galaxy_t   * new_galaxy(run_globals_t *run_globals, int snapshot, int halo_ID);
void         create_new_galaxy(run_globals_t *run_globals, int snapshot, halo_t *halo, int *NGal, int *new_gal_counter);
void         assign_galaxy_to_halo(galaxy_t *gal, halo_t *halo);
void         kill_galaxy(run_globals_t *run_globals, galaxy_t *gal, galaxy_t *prev_gal, int *NGal, int *kill_counter);
void         copy_halo_to_galaxy(halo_t *halo, galaxy_t *gal, int snapshot);
void         reset_galaxy_properties(galaxy_t *gal);
double       gas_infall(run_globals_t *run_globals, fof_group_t *FOFgroup, int snapshot);
void         add_infall_to_hot(galaxy_t *central, double infall_mass);
double       calculate_merging_time(run_globals_t *run_globals, galaxy_t *gal, int snapshot);
void         merge_with_target(run_globals_t *run_globals, galaxy_t *gal, int *dead_gals, int snapshot);
void         insitu_star_formation(run_globals_t *run_globals, galaxy_t *gal, int snapshot);
void         update_reservoirs_from_sf(run_globals_t *run_globals, galaxy_t *gal, double new_stars);
void         delayed_supernova_feedback(run_globals_t *run_globals, galaxy_t *gal, int snapshot);
void         contemporaneous_supernova_feedback(run_globals_t *run_globals, galaxy_t *gal, double *m_stars, int snapshot, double *m_reheat, double *m_eject, double *m_recycled, double *new_metals);
void         update_reservoirs_from_sn_feedback(galaxy_t *gal, double m_reheat, double m_eject, double m_recycled, double new_metals);
void         prep_hdf5_file(run_globals_t *run_globals);
void         create_master_file(run_globals_t *run_globals);
void         write_snapshot(run_globals_t *run_globals, int n_write, int i_out, int *last_n_write, trees_info_t *trees_info);
void         calc_hdf5_props(run_globals_t *run_globals);
void         prepare_galaxy_for_output(run_globals_t *run_globals, galaxy_t gal, galaxy_output_t *galout, int i_snap);
void         read_photometric_tables(run_globals_t *run_globals);
int          compare_ints(const void *a, const void *b);
void         mpi_debug_here(void);
void         check_counts(run_globals_t *run_globals, fof_group_t *fof_group, int NGal, int NFof);
void         cn_quote(void);
double       Tvir_to_Mvir(run_globals_t *run_globals, double T, double z);
double       calculate_Mvir(run_globals_t *run_globals, double Mvir, int len);
double       calculate_Rvir(run_globals_t *run_globals, double Mvir, int snapshot);
double       calculate_Vvir(run_globals_t *run_globals, double Mvir, double Rvir);
double       calculate_spin_param(halo_t *halo);
void         read_cooling_functions(run_globals_t *run_globals);
double       interpolate_cooling_rate(double logTemp, double logZ);
double       gas_cooling(run_globals_t *run_globals, galaxy_t *gal);
void         cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass);
double       calc_metallicity(double total_gas, double metals);
void         reincorporate_ejected_gas(run_globals_t *run_globals, galaxy_t *gal);
double       radio_mode_BH_heating(run_globals_t *run_globals, galaxy_t *gal, double cooling_mass);
void         merger_driven_BH_growth(run_globals_t *run_globals, galaxy_t *gal, double merger_ratio);


// Magnitude related
void   init_luminosities(run_globals_t *run_globals, galaxy_t *gal);
void   add_to_luminosities(run_globals_t *run_globals, galaxy_t *gal, double burst_mass, double metallicity, double burst_time);
double lum_to_mag(double lum);
void   sum_luminosities(run_globals_t *run_globals, galaxy_t *parent, galaxy_t *gal, int outputbin);
void   prepare_magnitudes_for_output(run_globals_t *run_globals, galaxy_t gal, galaxy_output_t *galout, int i_snap);
void   apply_dust(int n_photo_bands, galaxy_t gal, double *LumDust, int outputbin);
void   cleanup_mags(run_globals_t *run_globals);

// Reionization related
double reionization_modifier(run_globals_t *run_globals, double Mvir, float *Pos, int snapshot);
double sobacchi2013_modifier(run_globals_t *run_globals, double Mvir, double redshift);
double gnedin2000_modifer(run_globals_t *run_globals, double Mvir, double redshift);
double global_ionizing_emmisivity(run_globals_t *run_globals);
#ifdef USE_TOCF
double tocf_modifier(run_globals_t *run_globals, double Mvir, float *Pos, int snapshot);
void set_HII_eff_factor(run_globals_t *run_globals);
int  find_cell(float pos, double box_size);
void malloc_reionization_grids(run_globals_t *run_globals);
void free_reionization_grids(run_globals_t *run_globals);
void construct_stellar_grids(run_globals_t *run_globals);
// void    assign_ionization_to_halos(run_globals_t *run_globals, halo_t *halo, int n_halos, float *xH_grid, int xH_dim);
int  read_dm_grid(run_globals_t *run_globals, int snapshot, int i_grid, float *grid);
void calculate_Mvir_crit(run_globals_t *run_globals, double redshift);
void call_find_HII_bubbles(run_globals_t *run_globals, int snapshot, int nout_gals);
void save_tocf_grids(run_globals_t *run_globals, hid_t group_id, int snapshot);
void check_if_reionization_complete(run_globals_t *run_globals);
#endif

#ifdef DEBUG
int debug(const char * restrict format, ...);
void check_pointers(run_globals_t *run_globals, halo_t *halos, fof_group_t *fof_groups, trees_info_t *trees_info);
#endif
#endif // _INIT_MERAXES
