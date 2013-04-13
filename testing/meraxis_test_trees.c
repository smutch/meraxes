#include <stdio.h>
#include <string.h>

#define STRLEN 256

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

static void write_trees_header(
  FILE *fout,            
  int   n_groups,        
  int   n_subgroups,     
  int   n_halos_max,     
  int   n_trees_subgroup,
  int   n_trees_group)   
{
  fwrite(&n_groups         , 1, sizeof(int), fout);
  fwrite(&n_subgroups      , 1, sizeof(int), fout);
  fwrite(&n_halos_max      , 1, sizeof(int), fout);
  fwrite(&n_trees_subgroup , 1, sizeof(int), fout);
  fwrite(&n_trees_group    , 1, sizeof(int), fout);
}

static void write_catalogs_header(
    FILE *fout,
    int i_file,
    int N_files,
    int N_halos_file,
    int N_halos_total)
{
  fwrite(&i_file       , 1, sizeof(int), fout);
  fwrite(&N_files      , 1, sizeof(int), fout);
  fwrite(&N_halos_file , 1, sizeof(int), fout);
  fwrite(&N_halos_total, 1, sizeof(int), fout);
}

static void write_group(
    FILE *fout,
    int id,
    int tree_flags,
    int desc_id,
    int file_offset,
    int file_index,
    int n_subgroups)
{
  int dummy = 0;

  fwrite(&id         , 1, sizeof(int), fout);
  fwrite(&tree_flags , 1, sizeof(int), fout);
  fwrite(&desc_id    , 1, sizeof(int), fout);
  fwrite(&dummy      , 1, sizeof(int), fout);
  fwrite(&file_offset, 1, sizeof(int), fout);
  fwrite(&file_index , 1, sizeof(int), fout);
  fwrite(&n_subgroups, 1, sizeof(int), fout);
}

static void write_subgroup(
    FILE *fout,
    int id,
    int tree_flags,
    int desc_id,
    int file_offset,
    int file_index)
{
  int dummy = 0;

  fwrite(&id         , 1, sizeof(int), fout);
  fwrite(&tree_flags , 1, sizeof(int), fout);
  fwrite(&desc_id    , 1, sizeof(int), fout);
  fwrite(&dummy      , 1, sizeof(int), fout);
  fwrite(&file_offset, 1, sizeof(int), fout);
  fwrite(&file_index , 1, sizeof(int), fout);
}

static void write_catalog_entry(
  FILE      *fout,       
  long long  id_MBP,     
  double     Mvir,       
  int        n_particles,
  float      R_vir,      
  float      V_max)      
{
  catalog_halo_struct entry;

  entry.id_MBP                    = id_MBP;
  entry.M_vir                     = Mvir;
  entry.n_particles               = n_particles;
  entry.position_COM[0]           = 0;
  entry.position_COM[1]           = 1;
  entry.position_COM[2]           = 2;
  entry.position_MBP[0]           = 0;
  entry.position_MBP[1]           = 1;
  entry.position_MBP[2]           = 2;
  entry.velocity_COM[0]           = 0;
  entry.velocity_COM[1]           = 1;
  entry.velocity_COM[2]           = 2;
  entry.velocity_MBP[0]           = 0;
  entry.velocity_MBP[1]           = 1;
  entry.velocity_MBP[2]           = 2;
  entry.R_vir                     = R_vir;
  entry.R_halo                    = R_vir;
  entry.R_max                     = R_vir;
  entry.V_max                     = V_max;
  entry.sigma_v                   = 100.0;
  entry.spin[0]                   = 0;
  entry.spin[1]                   = 1;
  entry.spin[2]                   = 2;
  entry.q_triaxial                = 1;
  entry.s_triaxial                = 1;
  entry.shape_eigen_vectors[0][0] = 0;
  entry.shape_eigen_vectors[0][1] = 0;
  entry.shape_eigen_vectors[0][2] = 0;
  entry.shape_eigen_vectors[1][0] = 0;
  entry.shape_eigen_vectors[1][1] = 0;
  entry.shape_eigen_vectors[1][2] = 0;
  entry.shape_eigen_vectors[2][0] = 0;
  entry.shape_eigen_vectors[2][1] = 0;
  entry.shape_eigen_vectors[2][2] = 0;

  fwrite(&entry, 1, sizeof(catalog_halo_struct), fout);
}

int main(int argc, char const* argv[])
{

  int snapshot;
  char trees_fname_base[STRLEN] = "halos/test/trees/horizontal/trees/test_step001_scan001.trees_horizontal";
  char fname[STRLEN];
  FILE *fout;
  FILE *fout_groups;
  FILE *fout_subgroups;

  /*
   * Snapshot 0
   */
  snapshot = 0;

  // Start with the trees
  sprintf(fname, "%s_%d", trees_fname_base, snapshot);
  fout = fopen(fname, "wb");

  write_trees_header(fout, 3, 3, 6, 3, 3);
  write_group(fout   , 0, 1, 0, 1, 0  , 1);
  write_subgroup(fout, 0, 1, 0, 1, 0);
  write_group(fout   , 1, 1, 1, 1, 1  , 1);
  write_subgroup(fout, 1, 1, 1, 1, 1);
  write_group(fout   , 2, 1, 2, 1, 2  , 1);
  write_subgroup(fout, 2, 1, 2, 1, 2);
  fclose(fout);

  // now deal with the catalogs
  sprintf(fname, "halos/test/catalogs/subfind_%03d.catalog_groups_properties", snapshot);
  fout_groups = fopen(fname, "wb");
  sprintf(fname, "halos/test/catalogs/subfind_%03d.catalog_subgroups_properties", snapshot);
  fout_subgroups = fopen(fname, "wb");

  write_catalogs_header(fout_groups   , 0, 1, 3, 3);
  write_catalogs_header(fout_subgroups, 0, 1, 3, 3);
  write_catalog_entry(fout_groups   , 100 , 6e9, 600, 1, 80);
  write_catalog_entry(fout_subgroups, 100 , 5e9, 500, 1, 80);
  write_catalog_entry(fout_groups   , 1300, 8e9, 800, 1, 90);
  write_catalog_entry(fout_subgroups, 1300, 7e9, 700, 1, 90);
  write_catalog_entry(fout_groups   , 2200, 9e9, 900, 1, 95);
  write_catalog_entry(fout_subgroups, 2200, 8e9, 800, 1, 95);
  fclose(fout_groups);
  fclose(fout_subgroups);

  /*
   * Snapshot 1
   */
  snapshot = 1;

  // Start with the trees
  sprintf(fname, "%s_%d", trees_fname_base, snapshot);
  fout = fopen(fname, "wb");

  write_trees_header(fout, 2, 3  , 6, 3, 3);
  write_group(fout       , 0, 1  , 0, 1, 0  , 1);
  write_subgroup(fout    , 0, 1  , 0, 1, 0);
  write_group(fout       , 1, 1  , 1, 1, 1  , 2);
  write_subgroup(fout    , 1, 1|2, 1, 1, 1);
  write_subgroup(fout    , 2, 1|2, 1, 1, 1);
  fclose(fout);

  // now deal with the catalogs
  sprintf(fname, "halos/test/catalogs/subfind_%03d.catalog_groups_properties", snapshot);
  fout_groups = fopen(fname, "wb");
  sprintf(fname, "halos/test/catalogs/subfind_%03d.catalog_subgroups_properties", snapshot);
  fout_subgroups = fopen(fname, "wb");

  write_catalogs_header(fout_groups   , 0, 1, 2, 2);
  write_catalogs_header(fout_subgroups, 0, 1, 3, 3);
  write_catalog_entry(fout_groups   , 100 , 6e9   , 600 , 1, 80);
  write_catalog_entry(fout_subgroups, 100 , 5e9   , 500 , 1, 80);
  write_catalog_entry(fout_groups   , 1300, 1.7e10, 1800, 1, 120);
  write_catalog_entry(fout_subgroups, 1300, 8e9   , 800 , 1, 95);
  write_catalog_entry(fout_subgroups, 2200, 7e9   , 700 , 1, 90);
  fclose(fout_groups);
  fclose(fout_subgroups);

  /*
   * Snapshot 2
   */
  snapshot = 2;

  // Start with the trees
  sprintf(fname, "%s_%d", trees_fname_base, snapshot);
  fout = fopen(fname, "wb");

  write_trees_header(fout, 2, 2, 6, 3, 3);
  write_group(fout   , 0, 1, 0, 1, 0  , 1);
  write_subgroup(fout, 0, 1, 0, 1, 0);
  write_group(fout   , 1, 1, 1, 1, 1  , 2);
  write_subgroup(fout, 1, 1, 1, 1, 1);
  fclose(fout);

  // now deal with the catalogs
  sprintf(fname, "halos/test/catalogs/subfind_%03d.catalog_groups_properties", snapshot);
  fout_groups = fopen(fname, "wb");
  sprintf(fname, "halos/test/catalogs/subfind_%03d.catalog_subgroups_properties", snapshot);
  fout_subgroups = fopen(fname, "wb");

  write_catalogs_header(fout_groups   , 0, 1, 2, 2);
  write_catalogs_header(fout_subgroups, 0, 1, 2, 2);
  write_catalog_entry(fout_groups   , 100 , 6e9   , 600 , 1, 80);
  write_catalog_entry(fout_subgroups, 100 , 5e9   , 500 , 1, 80);
  write_catalog_entry(fout_groups   , 1300, 1.7e10, 1800, 1, 120);
  write_catalog_entry(fout_subgroups, 100 , 1.5e10, 1500, 1, 120);
  fclose(fout_groups);
  fclose(fout_subgroups);

  return 0;
}
