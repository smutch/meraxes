#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gsl/gsl_rng.h>

/*
 * Definitions
 */

#define STRLEN  256  //!< Default string length
#define MAXTAGS 50   //!< Maximum number of allowed tags in input file

// TODO: This should not be hard coded!
#define MAXSNAPS 10  //!< Maximum number of snapshots

#define ABORT(sigterm)                                                                 \
do {                                                                                   \
  SID_log_error("in file: %s\tfunc: %s\tline: %i", __FILE__, __FUNCTION__, __LINE__); \
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
struct physics_params_struct{
  int    funcprop;
  double peak;
  double sigma;
  double stellarfrac;
  double peak_evo;
  double sigma_evo;
  double stellarfrac_evo;
  double bhgrowthfactor;
};
typedef struct physics_params_struct physics_params_struct;

//! Run params
//! Everything in this structure is supplied by the user...
struct run_params_struct{
  char                  filename[STRLEN];
  char                  OutputDir[STRLEN];
  char                  FileNameGalaxies[STRLEN];
  char                  SimulationDir[STRLEN];
  char                  PhotometricTabsDir[STRLEN];
  char                  CoolFunctionsDir[STRLEN];
  char                  SimulationFilePrefix[STRLEN];
  char                  FileWithOutputSnaps[STRLEN];
  char                  FileWithSnapList[STRLEN];
  int                   FilesPerSnapshot;
  int                   LastSnapShotNr;
  int                   FirstFile;
  int                   LastFile;
  double                BoxSize;
  double                VolumeFactor;
  double                ThreshMajorMerger;
  double                RecycleFraction;
  double                UnitVelocity_in_cm_per_s;
  double                UnitLength_in_cm;
  double                UnitMass_in_g;
  double                SimHubble_h;
  double                ObsHubble_h;
  int                   DiskInstabilityOn;
  double                BaryonFrac;
  double                Omega;
  double                OmegaLambda;
  double                PartMass;
  double                MergerTimeFactor;
  int                   SnaplistLength;
  physics_params_struct physics;
};
typedef struct run_params_struct run_params_struct;

//! Global variables which will will be passed around
struct run_globals_struct{
  gsl_rng               *random_generator;
  double                AA[MAXSNAPS];
  double                ZZ[MAXSNAPS];
  double                Age[MAXSNAPS];
  run_params_struct params;
};
typedef struct run_globals_struct run_globals_struct;

//! The header from the input tree files.
struct trees_header_struct{
  int n_groups;
  int n_subgroups;
  int n_halos_max;
  int n_trees_subgroup;
  int n_trees_group;
};
typedef struct trees_header_struct trees_header_struct;

//! This is the structure for a halo in the catalog files
struct catalog_halo_struct{
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
typedef struct catalog_halo_struct catalog_halo_struct;


//! The meraxis halo structure
struct halo_struct{
  int    id;             //!< Halo ID
  int    type;           //!< Type (0 for central, 1 for satellite)
  int    desc_id;        //!< Descendant ID
  int    file_offset;    //!< Number of snapshots until the descendant of this halo reappears
  int    file_index;     //!< Index of descendant in next relevant snapshot
  int    tree_flags;     //!< Bitwise flag indicating the type of match in the trees
  int    n_subgroups;    //!< Number of subgroups belonging to this type 0 (=-1 if type=1)
  int    len;            //!< Number of satellites belonging to this halo (-1 if halo is satellite itself)
  double Mvir;           //!< Bryan &Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int    n_particles;    //!< Number of particles in the structure
  float  position[3];    //!< Most bound particle position [Mpc/h]
  float  velocity[3];    //!< Centre-of-mass velocity [km/s]
  float  Rvir;           //!< Virial radius [Mpc/h]
  float  Rhalo;          //!< Distance of last halo particle from MBP [Mpc/h]
  float  Rmax;           //!< Radius of maximum circular velocity [Mpc/h]
  float  Vmax;           //!< Maximum circular velocity [km/s]
  float  VelDisp;        //!< Total 3D velocity dispersion [km/s]
  float  spin[3];        //!< Specific angular momentum vector [Mpc/h *km/s]
};
typedef struct halo_struct halo_struct;


struct galaxy_struct
{
  int   Type;
  int   CentralGal;
  float CentralMvir;

  // properties of subhalo at the last time this galaxy was a central galaxy
  float Pos[3];
  float Vel[3];
  int   Len;
  float Mvir;
  float dM;
  float dMdt;
  float Rvir;
  float Vvir;
  float Vmax;

  // baryonic reservoirs
  float StellarMass;
  float BulgeMass;
  float BlackHoleMass;

  // misc
  float Sfr[NOUT];
  float SfrBulge[NOUT];
  float DiskRadius;
  float Cos_Inc;
  float MergTime;
};
typedef struct galaxy_struct galaxy_struct;


struct galaxy_output_struct
{
  int   Type;
  int   HaloIndex;
  int   SnapNum;
  int   CentralGal;
  float CentralMvir;

  // properties of subhalo at the last time this galaxy was a central galaxy
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
  float VelDisp;

  // baryonic reservoirs
  float StellarMass;
  float BulgeMass;
  float BlackHoleMass;

  // misc
  float Sfr;
  float SfrBulge;
  float DiskRadius;
  float Cos_Inc;
  float MergTime;
};
typedef struct galaxy_output_struct galaxy_output_struct;


/*
 * Functions
 */

void myexit(int signum);
void read_parameter_file(run_globals_struct *run_globals, char *fname);
void init_meraxis(run_globals_struct *run_globals);
void dracarys(run_globals_struct *run_globals);
