#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <gbpLib.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>
#include <hdf5.h>
#include <fftw3.h>

#ifndef _INIT_MERAXES
#define _INIT_MERAXES

/*
 * Definitions
 */


#ifdef CALC_MAGS
#ifndef NOUT
#define NOUT 1
#endif
#endif

#ifndef N_HISTORY_SNAPS
#define N_HISTORY_SNAPS 5
#endif

// ======================================================
// Don't change these unless you know what you are doing!
#define STRLEN  256  //!< Default string length
#ifndef MAX_PHOTO_NBANDS
#define MAX_PHOTO_NBANDS 5
#endif
#ifndef N_PHOTO_JUMPS
#define N_PHOTO_JUMPS 1000
#endif
// ======================================================

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
#define PLANCK           6.6262e-27  //! [erg/s]
#define PROTONMASS       1.6726e-24
#define HUBBLE           3.2407789e-18 //! [h/sec]
#define SEC_PER_MEGAYEAR 3.155e13
#define SEC_PER_YEAR     3.155e7
#define MPC              3.086e24
#define TCMB             2.728

// Constants
#define REL_TOL (float)1e-5
#define ABS_TOL (float)1e-8

#ifdef DEBUG
FILE *meraxes_debug_file;
#endif

#define PHYSICS_FLAG_MAXIMAL_COOLING 1
#define PHYSICS_FLAG_MAXIMAL_INSITU_SF 2
#define PHYSICS_FLAG_MAXIMAL_MERGER_SF 4

/*
 * Enums
 */
typedef enum index_type {
  INDEX_PADDED = 5674,
  INDEX_REAL,
  INDEX_COMPLEX_HERM,
} index_type;

/*
 * Structures
 */

//! Physics parameter values
typedef struct physics_params_t {
  double SfEfficiency;
  double SfCriticalSDNorm;
  double SfRecycleFraction;
  double SnReheatEff;
  double SnReheatLimit;
  double SnReheatScaling;
  double SnReheatNorm;
  double SnEjectionEff;
  double SnEjectionScaling;
  double SnEjectionNorm;
  double MaxCoolingMassFactor;
  double ReincorporationEff;
  double Yield;
  double IMFSlope;
  double EnergyPerSN;
  double IMFNormConst;
  double RadioModeEff;
  double BlackHoleGrowthRate;

  double ThreshMajorMerger;
  double MinMergerStellarMass;
  double MinMergerRatioForBurst;
  double MergerBurstScaling;
  double MergerBurstFactor;
  double MergerTimeFactor;

  // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value
  double ReionEfficiency;
  double ReionNionPhotPerBary;
  double ReionEscapeFrac;
  double ReionTcool;
  double Y_He;

  double  ReionGammaHaloBias;
  double  ReionAlphaUV;
  double  ReionRBubbleMin;
  double  ReionRBubbleMax;

  // global reionization prescription
  double ReionSobacchi_Zre;
  double ReionSobacchi_DeltaZre;
  double ReionSobacchi_DeltaZsc;
  double ReionSobacchi_T0;

  // global reionization prescription
  double ReionGnedin_z0;
  double ReionGnedin_zr;

  // filtering mass fit
  double  ReionSMParam_m0;
  double  ReionSMParam_a;
  double  ReionSMParam_b;
  double  ReionSMParam_c;
  double  ReionSMParam_d;

  // Flags
  int Flag_RedshiftDepEscFrac;
  int Flag_ReionizationModifier;
  int Flag_BHFeedback;
  int Flag_IRA;
  int Flag_FixDiskRadiusOnInfall;
  int Flag_FixVmaxOnInfall;
  int Flag_ReheatToFOFGroupTemp;
  int SfDiskVelOpt;

} physics_params_t;

//! Run params
//! Everything in this structure is supplied by the user...
typedef struct run_params_t {
  char             DefaultsFile[STRLEN];
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
  char             ForestIDFile[STRLEN];
  char             MvirCritFile[STRLEN];

  physics_params_t physics;

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
  long long        NPart;

  double          *MvirCrit;


  double           ReionDeltaRFactor;
  double           ReionPowerSpecDeltaK;
  int              ReionGridDim;
  int              ReionFilterType;
  int              ReionRtoMFilterType;
  int              ReionUVBFlag;

  int              FirstFile;
  int              LastFile;
  int              NSteps;
  int              SnaplistLength;
  int              RandomSeed;
  int              FlagSubhaloVirialProps;
  int              FlagInteractive;
  int              FlagGenDumpFile;
  int              FlagReadDumpFile;
  int              FlagMCMC;
  int              Flag_PatchyReion;
} run_params_t;


typedef struct run_units_t {
  double UnitTime_in_s;
  double UnitLength_in_cm;
  double UnitVelocity_in_cm_per_s;
  double UnitTime_in_Megayears;
  double UnitMass_in_g;
  double UnitDensity_in_cgs;
  double UnitPressure_in_cgs;
  double UnitCoolingRate_in_cgs;
  double UnitEnergy_in_cgs;
  // TOTAL : 72  (must be multiple of 8)
} run_units_t;

typedef struct hdf5_output_t {
  char       **params_tag;
  void       **params_addr;
  int         *params_type;
  size_t      *dst_offsets;
  size_t      *dst_field_sizes;
  const char **field_names;
  const char **field_units;
  const char **field_h_conv;
  hid_t       *field_types;
  size_t       dst_size;
  hid_t        array3f_tid;  // sizeof(hid_t) = 4
  hid_t        array_nmag_f_tid;
  hid_t        array_nhist_f_tid;
  int          n_props;
  int          params_count;

  // TOTAL : 52 + 4 padding (must be multiple of 8)
} hdf5_output_t;

typedef struct phototabs_t {
  int    JumpTable[N_PHOTO_JUMPS];
  char  (*MagBands)[5];
  float *Table;
  float *Ages;
  float *Metals;
  float  JumpFactor;
  int    NAges;
  int    NBands;
  int    NMetals;
} phototabs_t;


typedef struct gal_to_slab_t {
  int index;
  struct galaxy_t *galaxy;
  int slab_ind;
} gal_to_slab_t;


typedef struct reion_grids_t {
  ptrdiff_t     *slab_nix;
  ptrdiff_t     *slab_ix_start;
  ptrdiff_t     *slab_n_complex;

  float         *buffer;
  float         *stars;
  fftwf_complex *stars_unfiltered;
  fftwf_complex *stars_filtered;
  float         *deltax;
  fftwf_complex *deltax_unfiltered;
  fftwf_complex *deltax_filtered;
  float         *sfr;
  fftwf_complex *sfr_unfiltered;
  fftwf_complex *sfr_filtered;
  float         *xH;
  float         *z_at_ionization;
  float         *J_21_at_ionization;
  float         *J_21;
  float         *Mvir_crit;
  float         *r_bubble;
  gal_to_slab_t *galaxy_to_slab_map;

  double         global_xH;
  bool           reion_complete;
  int            buffer_size;
} reion_grids_t;

//! The meraxis halo structure
typedef struct halo_t {
  long long           id_MBP; //!< ID of most bound particle
  struct fof_group_t *FOFGroup;
  struct halo_t      *NextHaloInFOFGroup;
  struct galaxy_t    *Galaxy;

  float               Pos[3]; //!< Most bound particle position [Mpc/h]
  float               Vel[3]; //!< Centre-of-mass velocity [km/s]
  float               AngMom[3]; //!< Specific angular momentum vector [Mpc/h *km/s]
  float               Rhalo; //!< Distance of last halo particle from MBP [Mpc/h]

  double              Mvir; //!< virial mass [M_sol/h]
  double              Rvir; //!< Virial radius [Mpc/h]
  double              Vvir; //!< Virial velocity [km/s]

  float               Rmax; //!< Radius of maximum circular velocity [Mpc/h]
  float               Vmax; //!< Maximum circular velocity [km/s]
  float               VelDisp; //!< Total 3D velocity dispersion [km/s]
  int                 ID; //!< Halo ID
  int                 Type; //!< Type (0 for central, 1 for satellite)
  int                 SnapOffset; //!< Number of snapshots this halo skips before reappearing
  int                 DescIndex; //!< Index of descendant in next relevant snapshot
  int                 TreeFlags; //!< Bitwise flag indicating the type of match in the trees
  int                 ForestID;
  int                 Len; //!< Number of particles in the structure
} halo_t;

typedef struct fof_group_t {
  halo_t *FirstHalo;
  halo_t *FirstOccupiedHalo;
  double  Mvir;
  double  Rvir;
  double  Vvir;
  int     TotalSubhaloLen;
} fof_group_t;

typedef struct galaxy_t {
#ifdef CALC_MAGS
  double Lum[MAX_PHOTO_NBANDS][NOUT];
#endif

  double NewStars[N_HISTORY_SNAPS];

  // Unique ID for the galaxy
  long long ID;

  // properties of subhalo at the last time this galaxy was a central galaxy
  float  Pos[3];
  float  Vel[3];
  double Mvir;
  double Rvir;
  double Vvir;
  double Vmax;
  double Spin;

  double dt; //!< Time between current snapshot and last identification

  long long id_MBP;
  struct halo_t   *Halo;
  struct galaxy_t *FirstGalInHalo;
  struct galaxy_t *NextGalInHalo;
  struct galaxy_t *Next;
  struct galaxy_t *MergerTarget;

  // baryonic reservoirs
  double HotGas;
  double MetalsHotGas;
  double ColdGas;
  double MetalsColdGas;
  double Mcool;
  double StellarMass;
  double GrossStellarMass;
  double MetalsStellarMass;
  double DiskScaleLength;
  double Sfr;
  double EjectedGas;
  double MetalsEjectedGas;
  double BlackHoleMass;

  // baryonic hostories
  double mwmsa_num;
  double mwmsa_denom;

  // misc
  double Rcool;
  double Cos_Inc;
  double MergTime;
  double MergerStartRadius;
  double BaryonFracModifier;
  double MvirCrit;

  int  Type;
  int  OldType;
  int  Len;
  int  MaxLen;
  int  SnapSkipCounter;
  int  HaloDescIndex;
  int  TreeFlags;
  int  LastIdentSnap;     //!< Last snapshot at which the halo in which this galaxy resides was identified
  int  output_index;  //!< write index
  int  PhysicsFlags;

  bool ghost_flag;

  // N.B. There will be padding present at the end of this struct, but amount
  // is dependent on CALC_MAGS, MAX_PHOTO_NBANDS, NOUT and N_HISTORY_SNAPS.
} galaxy_t;

typedef struct galaxy_output_t {
  long long id_MBP;

  // Unique ID for the galaxy
  long long ID;

#ifdef CALC_MAGS
  float Mag[MAX_PHOTO_NBANDS];
  float MagDust[MAX_PHOTO_NBANDS];
#endif

  int Type;
  int CentralGal;
  int GhostFlag;
  int Len;
  int MaxLen;

  float Pos[3];
  float Vel[3];
  float Spin;
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
  float GrossStellarMass;
  float MetalsStellarMass;
  float Sfr;
  float EjectedGas;
  float MetalsEjectedGas;
  float BlackHoleMass;

  // misc
  float Rcool;
  float Cos_Inc;
  float MergTime;
  float MergerStartRadius;
  float BaryonFracModifier;
  float MvirCrit;
  int   PhysicsFlags;

  // baryonic histories
  float MWMSA;  // Mass weighted mean stellar age
  float NewStars[N_HISTORY_SNAPS];
} galaxy_output_t;


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
  int id;
  int flags;
  int desc_id;
  int tree_id;
  int file_offset;
  int desc_index;
  int central_index;
  int forest_id;
  int group_index;
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
  struct run_params_t params;
  char                FNameOut[STRLEN];
  reion_grids_t       reion_grids;
  struct run_units_t  units;
  hdf5_output_t       hdf5props;

  double             *AA;
  double             *ZZ;
  double             *LTTime;
  int                *RequestedForestId;
  int                *ListOutputSnaps;
  halo_t            **SnapshotHalo;
  fof_group_t       **SnapshotFOFGroup;
  int               **SnapshotIndexLookup;
  float             **SnapshotDeltax;
  trees_info_t       *SnapshotTreesInfo;
  phototabs_t        *photo;
  struct galaxy_t    *FirstGal;
  struct galaxy_t    *LastGal;
  gsl_rng            *random_generator;
  void               *mhysa_self;
  double              Hubble;
  double              RhoCrit;
  double              G;

  int                 NOutputSnaps;
  int                 LastOutputSnap;
  int                 NGhosts;
  int                 NHalosMax;
  int                 NFOFGroupsMax;
  int                 NRequestedForests;
  int                 TreesStep;
  int                 TreesScan;
  int                 NStoreSnapshots;

  bool                SelectForestsSwitch;
} run_globals_t;
#ifdef _MAIN
run_globals_t run_globals;
#else
extern run_globals_t run_globals;
#endif


/*
 * Functions
 */

void         myexit(int signum);
void         cleanup();
void         read_parameter_file(char *fname, int mode);
void         init_meraxes();
void         set_units();
void         continue_prompt(char *param_file);
void         free_halo_storage();
void         initialize_halo_storage();
void         dracarys();
int          evolve_galaxies(fof_group_t *fof_group, int snapshot, int NGal, int NFof);
void         passively_evolve_ghost(galaxy_t *gal, int snapshot);
trees_info_t read_halos(int snapshot, halo_t **halo, fof_group_t **fof_group, int **index_lookup, trees_info_t *snapshot_trees_info);
galaxy_t   * new_galaxy(int snapshot, int halo_ID);
void         create_new_galaxy(int snapshot, halo_t *halo, int *NGal, int *new_gal_counter);
void         assign_galaxy_to_halo(galaxy_t *gal, halo_t *halo);
void         kill_galaxy(galaxy_t *gal, galaxy_t *prev_gal, int *NGal, int *kill_counter);
void         copy_halo_to_galaxy(halo_t *halo, galaxy_t *gal, int snapshot);
void         reset_galaxy_properties(galaxy_t *gal, int snapshot);
double       gas_infall(fof_group_t *FOFgroup, int snapshot);
void         add_infall_to_hot(galaxy_t *central, double infall_mass);
double       calculate_merging_time(galaxy_t *gal, int snapshot);
void         merge_with_target(galaxy_t *gal, int *dead_gals, int snapshot);
void         insitu_star_formation(galaxy_t *gal, int snapshot);
void         update_reservoirs_from_sf(galaxy_t *gal, double new_stars);
double       sn_m_low(double log_dt);
double       calc_recycled_frac(double m_high, double m_low, double *burst_mass_frac);
void         delayed_supernova_feedback(galaxy_t *gal, int snapshot);
void         evolve_stellar_pops(galaxy_t *gal, int snapshot);
void         contemporaneous_supernova_feedback(galaxy_t *gal, double *m_stars, int snapshot, double *m_reheat, double *m_eject, double *m_recycled, double *new_metals);
void         update_reservoirs_from_sn_feedback(galaxy_t *gal, double m_reheat, double m_eject, double m_recycled, double new_metals);
void         prep_hdf5_file();
void         create_master_file();
void         write_snapshot(int n_write, int i_out, int *last_n_write, trees_info_t *trees_info);
void         calc_hdf5_props();
void         prepare_galaxy_for_output(galaxy_t gal, galaxy_output_t *galout, int i_snap);
void         read_photometric_tables();
int          compare_ints(const void *a, const void *b);
int          compare_floats(const void *a, const void *b);
int          compare_ptrdiff(const void *a, const void *b);
int          compare_slab_assign(const void *a, const void *b);
int          searchsorted(void *val,
                          void *arr,
                          int count,
                          size_t size,
                          int(*compare)(const void *a, const void *b),
                          int imin,
                          int imax);
float        comoving_distance(float a[3], float b[3]);
int          pos_to_ngp(double x, double side, int nx);
float        apply_pbc_pos(float x);
double       accurate_sumf(float *arr, int n);
int          grid_index(int i, int j, int k, int dim, int type);
void         mpi_debug_here(void);
int 				 isclosef(float a, float b, float rel_tol, float abs_tol);
void         printProgress (double percentage);
void         check_counts(fof_group_t *fof_group, int NGal, int NFof);
void         cn_quote(void);
double       Tvir_to_Mvir(double T, double z);
double       hubble_at_snapshot(int snapshot);
double       hubble_time(int snapshot);
double       calculate_Mvir(double Mvir, int len);
double       calculate_Rvir(double Mvir, int snapshot);
double       calculate_Vvir(double Mvir, double Rvir);
double       calculate_spin_param(halo_t *halo);
void         read_cooling_functions();
double       interpolate_cooling_rate(double logTemp, double logZ);
double       gas_cooling(galaxy_t *gal);
void         cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass);
double       calc_metallicity(double total_gas, double metals);
void         reincorporate_ejected_gas(galaxy_t *gal);
double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass);
void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio);

// Magnitude related
void   init_luminosities(galaxy_t *gal);
void   add_to_luminosities(galaxy_t *gal, double burst_mass, double metallicity, double burst_time);
double lum_to_mag(double lum);
void   sum_luminosities(galaxy_t *parent, galaxy_t *gal, int outputbin);
void   prepare_magnitudes_for_output(galaxy_t gal, galaxy_output_t *galout, int i_snap);
void   apply_dust(int n_photo_bands, galaxy_t gal, double *LumDust, int outputbin);
void   cleanup_mags();

// Reionization related
void   read_Mcrit_table();
double reionization_modifier(galaxy_t *gal, double Mvir, int snapshot);
double sobacchi2013_modifier(double Mvir, double redshift);
double gnedin2000_modifer(double Mvir, double redshift);
void   assign_slabs();

void   filter(fftwf_complex *box, int local_ix_start, int slab_nx, int grid_dim, float R);
double find_HII_bubbles(float redshift);
double tocf_modifier(galaxy_t *gal, double Mvir);
void   set_ReionEfficiency();
int    find_cell(float pos, double box_size);
void   malloc_reionization_grids();
void   free_reionization_grids();
int    map_galaxies_to_slabs(int ngals);
void   assign_Mvir_crit_to_galaxies(int ngals_in_slabs);
void   construct_baryon_grids(int snapshot, int ngals);
void   gen_grids_fname(int snapshot, char *name, bool relative);
void   create_grids_file();
int    read_dm_grid(int snapshot, int i_grid, float *grid);
void   free_grids_cache();
void   calculate_Mvir_crit(double redshift);
void   call_find_HII_bubbles(int snapshot, int unsampled_snapshot, int nout_gals);
void   save_reion_input_grids(int snapshot);
void   save_reion_output_grids(int snapshot);
bool   check_if_reionization_complete();
void   write_single_grid(const char *fname, float *grid, const char *grid_name, bool padded_flag, bool create_file_flag);


// MCMC related
// meraxes_mhysa_hook must be implemented by the calling code (Mhysa)!
int (*meraxes_mhysa_hook)(void *self, int snapshot, int ngals);

#ifdef DEBUG
int  debug(const char * restrict format, ...);
void check_pointers(halo_t *halos, fof_group_t *fof_groups, trees_info_t *trees_info);
#endif
#endif // _INIT_MERAXES
