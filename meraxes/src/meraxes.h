#include <stdio.h>
#include <stdlib.h>

/*
 * Structures
 */


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
  int id;                 //! Halo ID
  int type;               //! Type (0 for central, 1 for satellite)
  int desc_id;            //! Descendant ID
  int file_offset;        //! Number of snapshots until the descendant of this halo reappears
  int file_index;         //! Index of descendant in next relevant snapshot
  int tree_flags;         //! Bitwise flag indicating the type of match in the trees
  int n_satellites;       //! Number of satellites belonging to this halo (-1 if halo is satellite itself)
  double    M_vir;        //! Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int       n_particles;  //! Number of particles in the structure
  float     position[3];  //! Most bound particle position [Mpc/h]
  float     velocity[3];  //! Centre-of-mass velocity      [km/s]
  float     R_vir;        //! Virial radius [Mpc/h]
  float     R_halo;       //! Distance of last halo particle from MBP [Mpc/h]
  float     R_max;        //! Radius of maximum circular velocity [Mpc/h]
  float     V_max;        //! Maximum circular velocity [km/s]
  float     sigma_v;      //! Total 3D velocity dispersion [km/s]
  float     spin[3];      //! Specific angular momentum vector [Mpc/h*km/s]
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

