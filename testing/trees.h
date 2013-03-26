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
#define TREE_CASE_SIMPLE                        1      // Set when a halo has a file_offset=1 
#define TREE_CASE_MAIN_PROGENITOR               2      // Set for the progenitor with the highest match score
#define TREE_CASE_MERGER                        4      // Set when new IDs are created (ie. last point the halo was seen)
#define TREE_CASE_DROPPED                       8      // Set if file_offset>1 and TREE_CASE_MATCHED_TO_BRIDGE is not set
#define TREE_CASE_STRAYED                       16     // Set for halos for which a descendant was not found
#define TREE_CASE_SPUTTERED                     32     // Set for halos whose descendant was not given a valid ID
#define TREE_CASE_BRIDGED                       64     // Set for halos with multiple back-matches from halos with unique IDs
#define TREE_CASE_EMERGED_CANDIDATE             128    // Set when a halo is identified as a unique back-match to a halo marked TREE_CASE_BRIDGED 
                                                       //    and is not identified as the BRIDGE's main descendant
#define TREE_CASE_FOUND                         256    // Set when a halo has a progenitor with a file_offset>1
#define TREE_CASE_NO_PROGENITORS                512    // Set for halos that have no progenitors.
#define TREE_CASE_FRAGMENTED_LOST               1024   // Set for halos that are marked TREE_CASE_EMERGED_CANDIDATE, and whose 
                                                       //    decendant_id!=a valid id (ie they are not a progenitor of anything)
#define TREE_CASE_FRAGMENTED_RETURNED           2048   // Set for halos that are marked TREE_CASE_EMERGED_CANDIDATE, and whose 
                                                       //    decendant_id==the id of the halo they are emerged from
#define TREE_CASE_FRAGMENTED_EXCHANGED          4096   // Set for halos that are marked TREE_CASE_EMERGED_CANDIDATE, and whose 
                                                       //    decendant_id!=the id of the halo they are emerged but is nevertheless valid 
                                                       //    (ie. they are still a progenitor of something)
#define TREE_CASE_MATCHED_TO_BRIDGE             8192   // Set when a halo is matched to one with TREE_CASE_BRIDGED set
#define TREE_CASE_BRIDGE_DEFAULT                16384  // Set when a halo matched to a bridge is not matched to any emergent halos
#define TREE_CASE_MATCHED_TO_BRIDGE_UNPROCESSED 32768  // For internal use.  This should never be seen in the output.
#define TREE_CASE_BRIDGE_FINALIZE               65536  // For internal use.  This should never be seen in the output.
#define TREE_CASE_UNPROCESSED                   131072 // For internal use.  This should never be seen in the output.
#define TREE_CASE_INVALID                       262144 // For internal use.  This should never be seen in the output.
#define TREE_CASE_EMERGED                       (TREE_CASE_EMERGED_CANDIDATE+TREE_CASE_FOUND)
#define TREE_CASE_FRAGMENTED                    (TREE_CASE_EMERGED_CANDIDATE+TREE_CASE_NO_PROGENITORS)

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

