#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stddef.h>
#include <gbpSID.h>
#include <cn_exceptions.h>

#define STR_LEN 256

// Galaxy types
#define CENTRALT 0
#define SATELLITET 1
#define ORPHANT 2
#define DEADT 4

// Tree identification flags
#define TREE_CASE_SIMPLE                        1
#define TREE_CASE_STRAYED                       2
#define TREE_CASE_SPUTTERED                     4
#define TREE_CASE_DROPPED                       8
#define TREE_CASE_MERGER                        16
#define TREE_CASE_BRIDGED                       32
#define TREE_CASE_EMERGED                       64
#define TREE_CASE_BRIDGE_PROGENITOR             128
#define TREE_CASE_BRIDGE_PROGENITOR_UNPROCESSED 256
#define TREE_CASE_BRIDGE_FINALIZE               512
#define TREE_CASE_BRIDGE_DEFAULT                1024
#define TREE_CASE_FOUND                         2048
#define TREE_CASE_MAIN_PROGENITOR               4096
#define TREE_CASE_UNPROCESSED                   8192
#define TREE_CASE_INVALID                       16384

/*
 * STRUCTURES
 */
typedef struct TreesHeader TreesHeader;
struct TreesHeader{
  int n_groups;
  int n_subgroups;
  int n_halos_max;
  int n_trees_subgroup;
  int n_trees_group;  
};

// TODO: Check if this struct would benefit from padding...
typedef struct Halo Halo;
struct Halo{
  int id;
  int type;
  int desc_id;
  int tree_id;
  int file_offset;
  int file_index;
  int tree_flags;
  int n_subgroups;
  long long id_MBP;                    // ID of most bound particle in structure
  double    M_vir;                     // Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int       n_particles;               // Number of particles in the structure
  float     position_COM[3];           // Centre-of-mass position      [Mpc/h]
  float     position_MBP[3];           // Most bound particle position [Mpc/h]
  float     velocity_COM[3];           // Centre-of-mass velocity      [km/s]
  float     velocity_MBP[3];           // Most bound particle velocity [km/s]
  float     R_vir;                     // Virial radius [Mpc/h]
  float     R_halo;                    // Distance of last halo particle from MBP [Mpc/h]
  float     R_max;                     // Radius of maximum circular velocity     [Mpc/h]
  float     V_max;                     // Maximum circular velocity               [km/s]
  float     sigma_v;                   // Total 3D velocity dispersion            [km/s]
  float     spin[3];                   // Specific angular momentum vector        [Mpc/h*km/s]
  float     q_triaxial;                // Triaxial shape parameter q=b/a
  float     s_triaxial;                // Triaxial shape parameter s=c/a
  float     shape_eigen_vectors[3][3]; // Normalized triaxial shape eigenvectors
};

typedef struct Group Group;
struct Group{
  int id;
  int type;
  int desc_id;
  int tree_id;
  int file_offset;
  int file_index;
  int n_subgroups;
  Halo *halo;
};

typedef struct Subgroup Subgroup;
struct Subgroup{
  int id;
  int type;
  int desc_id;
  int tree_id;
  int file_offset;
  int file_index;
  Halo *halo;
};

typedef struct Galaxy Galaxy;
struct Galaxy {
  // Tree pointers
  Group *group;
  Subgroup *subgroup;
  Halo *halo;

  // Galaxy properties
  int type;
  double cold_gas;
 };


/*
 * FUNCTIONS
 */
TreesHeader read_trees(char *sim, int total_sim_snaps, int n_every_snaps, int n_scan_snaps, int snapshot,
                       Halo **halos);
void free_trees(Halo **halos);

