#include <stdio.h>
#include <stdlib.h>
#include <gbpLib.h>

/*
 * Definitions
 */

#define STRLEN  256  //!< Default string length
#define MAXTAGS 50   //!< Maximum number of allowed tags in input file

#define ABORT(sigterm)                                                                 \
do {                                                                                   \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
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
#define CM_PER_MPC       3.085678e24
#define PROTONMASS       1.6726e-24
#define HUBBLE           3.2407789e-18 //! [h/sec]
#define SEC_PER_MEGAYEAR 3.155e13
#define SEC_PER_YEAR     3.155e7


/*
 * Structures
 */

//! Physics parameter values
struct physics_params{
  double peak;
  double sigma;
  double stellarfrac;
  double peak_evo;
  double sigma_evo;
  double stellarfrac_evo;
};

//! Run params
struct run_params_struct{
  char OutputDir[STRLEN];
  char FileNameGalaxies[STRLEN];
  char SimulationDir[STRLEN];
  char PhotometricTabsDir[STRLEN];
  char CoolFunctionsDir[STRLEN];
  char SimulationFilePrefix[STRLEN];
  char FileWithOutputSnaps[STRLEN];
  char FileWithSnapList[STRLEN];
  int FilesPerSnapshot;
  int LastSnapShotNr;
  int FirstFile;
  int LastFile;
  double BoxSize;
  double VolumeFactor;
  double ThreshMajorMerger;
  double RecycleFraction;
  double UnitVelocity_in_cm_per_s;
  double UnitLength_in_cm;
  double UnitMass_in_g;
  double SimHubble_h;
  double ObsHubble_h;
  int DiskInstabilityOn;
  double BaryonFrac;
  double Omega;
  double OmegaLambda;
  double PartMass;
  double MergerTimeFactor;
  int funcprop;
  double peak;
  double peak_evo;
  double sigma;
  double sigma_evo;
  double stellarfrac;
  double stellarfrac_evo;
  double bhgrowthfactor;
  struct physics_params physics; 
};

#ifdef _MAIN
struct run_params_struct run_params;
#else
extern struct run_params_struct run_params;
#endif

//! The header from the input tree files.
struct trees_header_struct{
  int n_groups;
  int n_subgroups;
  int n_halos_max;
  int n_trees_subgroup;
  int n_trees_group;
};
typedef struct trees_header_struct trees_header_struct;


//! The halo structure
struct halo_struct{
  int    id;             //!< Halo ID
  int    type;           //!< Type (0 for central, 1 for satellite)
  int    desc_id;        //!< Descendant ID
  int    file_offset;    //!< Number of snapshots until the descendant of this halo reappears
  int    file_index;     //!< Index of descendant in next relevant snapshot
  int    tree_flags;     //!< Bitwise flag indicating the type of match in the trees
  int    n_satellites;   //!< Number of satellites belonging to this halo (-1 if halo is satellite itself)
  double M_vir;          //!< Bryan &Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int    n_particles;    //!< Number of particles in the structure
  float  position[3];    //!< Most bound particle position [Mpc/h]
  float  velocity[3];    //!< Centre-of-mass velocity [km/s]
  float  R_vir;          //!< Virial radius [Mpc/h]
  float  R_halo;         //!< Distance of last halo particle from MBP [Mpc/h]
  float  R_max;          //!< Radius of maximum circular velocity [Mpc/h]
  float  V_max;          //!< Maximum circular velocity [km/s]
  float  sigma_v;        //!< Total 3D velocity dispersion [km/s]
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

