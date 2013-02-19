#include "trees.h"

typedef struct RawHalo RawHalo;
struct RawHalo{
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


static void halo_catalog_filename(char *root, char *sim, int snapshot, char *group_type, int sub, int *i_layout, char *fname)
{
  bool flag_success = false;
  FILE *fin;

  // If we need to determine the filename structure...
  if (*i_layout==-1){
    for (*i_layout=0; (*i_layout<4) && (flag_success==false); (*i_layout)++){
      if (*i_layout==0)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, snapshot, group_type, sub);
      else if (*i_layout==1)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type, snapshot, group_type);
      else if (*i_layout==2)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, sub);
      else if (*i_layout==3)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type);
      
      // SID_log("i_layout = %d, Trying : %s", SID_LOG_COMMENT, *i_layout, fname);

      if ((fin = fopen(fname, "rb"))!=NULL){
        flag_success = true;
        fclose(fin);
        break;
      }
    }
  } 

  // Ensure we have a valid i_layout value.
  if (*i_layout<0 && *i_layout>3){
    fprintf(stderr, "Cannot resolve catalogue filename.\n");
    ABORT(EXIT_FAILURE);
  }
   
  // Provide the correct filename
  if (*i_layout==0)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, snapshot, group_type, sub);
  else if (*i_layout==1)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type, snapshot, group_type);
  else if (*i_layout==2)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, sub);
  else if (*i_layout==3)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type);

  // SID_log("FNAME : %s", SID_LOG_COMMENT, fname);

}

static void inline read_catalogs_header(FILE *fin, int *i_file, int *N_files, int *N_halos_file, int *N_halos_total)
{
  fread(i_file, sizeof(int), 1, fin);
  fread(N_files, sizeof(int), 1, fin);
  fread(N_halos_file, sizeof(int), 1, fin);
  fread(N_halos_total, sizeof(int), 1, fin);
}


static void inline read_halo(FILE** fin, char *root, char *sim, int snapshot, char *group_type, int *flayout_switch, 
                             int *i_file, int *N_halos_file,
                             int *i_halo, Halo *halos, int N_files, int *halo_count)
{
  char fname[STR_LEN];
  int dummy;
  RawHalo halo_in;

  // Is this the first read?
  if((*fin)==NULL){
    halo_catalog_filename(root, sim, snapshot, group_type, *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin==NULL)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, N_halos_file, &dummy);
  }

  // Have we already read all the halos in this file?
  if((*i_halo)>=(*N_halos_file)){
    fprintf(stderr, "i_halo = %d, N_halos_file = %d\n", (*i_halo), (*N_halos_file));
    fclose(*fin);
    (*i_file)++;
    (*i_halo) = 0;
    halo_catalog_filename(root, sim, snapshot, group_type, *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin==NULL)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, N_halos_file, &dummy);
  }

  // Read in a halo and then paste it into our storage array appropriately
  fread(&halo_in, sizeof(RawHalo), 1, *fin);
  memcpy(&(halos[*halo_count].id_MBP), &halo_in, offsetof(RawHalo, padding)); 
  
  // Update the halo counters
  (*halo_count)++;
}

static void inline read_trees_header(FILE *fin, TreesHeader *header)
{
  fread(&(header->n_groups)        , sizeof(int), 1, fin);
  fread(&(header->n_subgroups)     , sizeof(int), 1, fin);
  fread(&(header->n_halos_max)     , sizeof(int), 1, fin);
  fread(&(header->n_trees_subgroup), sizeof(int), 1, fin);
  fread(&(header->n_trees_group)   , sizeof(int), 1, fin);
}

static void inline read_group(FILE *fin, Halo *halos, int i_halo)
{
  fread(&(halos[i_halo].id)         , sizeof(int), 1, fin);
  fread(&(halos[i_halo].type)       , sizeof(int), 1, fin);
  fread(&(halos[i_halo].desc_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].tree_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].file_offset), sizeof(int), 1, fin);
  fread(&(halos[i_halo].n_subgroups), sizeof(int), 1, fin);

// #ifdef DEBUG
//   printf("INDEX %d - Group:\n", i_halo);
//   printf("\tid: %d\n", halos[i_halo].id);
//   printf("\ttype: %d\n", halos[i_halo].type);
//   printf("\tdesc_id: %d\n", halos[i_halo].desc_id);
//   printf("\ttree_id: %d\n", halos[i_halo].tree_id);
//   printf("\tfile_offset: %d\n", halos[i_halo].file_offset);
//   printf("\tn_subgroups: %d\n", halos[i_halo].n_subgroups);
// #endif
}

static void inline read_subgroup(FILE *fin, Halo *halos, int i_halo)
{
  fread(&(halos[i_halo].id)         , sizeof(int), 1, fin);
  fread(&(halos[i_halo].type)       , sizeof(int), 1, fin);
  fread(&(halos[i_halo].desc_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].tree_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].file_offset), sizeof(int), 1, fin);
  halos[i_halo].n_subgroups = -1;

// #ifdef DEBUG
//   printf("INDEX %d - Subgroup:\n", i_halo);
//   printf("\tid: %d\n", halos[i_halo].id);
//   printf("\ttype: %d\n", halos[i_halo].type);
//   printf("\tdesc_id: %d\n", halos[i_halo].desc_id);
//   printf("\ttree_id: %d\n", halos[i_halo].tree_id);
//   printf("\tfile_offset: %d\n", halos[i_halo].file_offset);
// #endif
}


TreesHeader read_trees(char *sim, int total_sim_snaps, int n_every_snaps, int n_scan_snaps, int snapshot,   // IN
                       Halo **halos)            // OUT
{

  int N_halos, N_halos_groups, N_halos_subgroups, N_groups_files, N_subgroups_files, dummy;
  TreesHeader header;
 
  // TODO: Sanity checks should go here...

  char sim_variant[18];
  sprintf(sim_variant, "step_%03d_scan_%03d", n_every_snaps, n_scan_snaps);
  SID_log("Reading snapshot %d (%s:%s) trees and halos...", SID_LOG_OPEN|SID_LOG_TIMER, snapshot, sim, sim_variant);


  // Calculate the actual (unsampled) simulation snapshot.
  int corrected_snapshot;
  if (n_every_snaps>1)
    corrected_snapshot = total_sim_snaps-1 - ((int)(total_sim_snaps/n_every_snaps))*n_every_snaps + snapshot*n_every_snaps;
  else
    corrected_snapshot = snapshot;

  // HALOS
  // Read the header info
  char fname[STR_LEN];
  FILE *fin, *fin_trees;
  int catalog_groups_flayout = -1;
  int catalog_subgroups_flayout = -1;
  halo_catalog_filename("data", sim, corrected_snapshot, "groups", 0, &catalog_groups_flayout, fname);
  // SID_log("DEBUG: fname = %s", SID_LOG_COMMENT, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(34494);
  }
  read_catalogs_header(fin, &dummy, &N_groups_files, &dummy, &N_halos_groups);
  fclose(fin);

  halo_catalog_filename("data", sim, corrected_snapshot, "subgroups", 0, &catalog_subgroups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(34494);
  }
  read_catalogs_header(fin, &dummy, &N_subgroups_files, &dummy, &N_halos_subgroups);
  fclose(fin);

  // Allocate the halo array
  // N_halos = N_halos_groups + N_halos_subgroups;
  N_halos = N_halos_subgroups;
  SID_log("N_halos_groups = %d", SID_LOG_COMMENT, N_halos_groups);
  SID_log("N_halos_subgroups = %d", SID_LOG_COMMENT, N_halos_subgroups);
  SID_log("N_halos = %d", SID_LOG_COMMENT, N_halos);
  if (N_halos>0)
    *halos = malloc(sizeof(Halo) * N_halos);

  // TREES
  SID_log("Reading in trees...", SID_LOG_COMMENT);
  sprintf(fname, "data/%s/trees/%s_%s/horizontal/trees/%s_%s.trees_horizontal_%d", sim, sim, sim_variant, sim, sim_variant, corrected_snapshot);
  // SID_log("DEBUG: fname = %s", SID_LOG_COMMENT, fname);
  fin_trees = fopen(fname, "rb");
  if (fin_trees==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(34494);
  }

  // Read the header info
  read_trees_header(fin_trees, &header);
  if (N_halos<1)
  {
    SID_log("No halos in this file... skipping...", SID_LOG_COMMENT);
    return header;
  }

  // Initialise everything we need to read the halos in
  FILE *fin_group_halos        = NULL;
  FILE *fin_subgroup_halos     = NULL;
  int   i_group_file           = 0;
  int   i_subgroup_file        = 0;
  int   halo_count             = 0;
  int   N_halos_groups_file    = 0;
  int   N_halos_subgroups_file = 0;
  int   group_count_infile     = 0;
  int   subgroup_count_infile  = 0;
  int   n_subgroups            = 0;
  int   group_count            = 0;

  // Loop through the groups and subgroups and read them in
  Halo group_halos[1];
  for (int i_group=0; i_group<header.n_groups; i_group++){
    read_group(fin_trees, group_halos, group_count);
    n_subgroups = group_halos[group_count].n_subgroups;
    read_halo(&fin_group_halos, "data", sim, corrected_snapshot, "groups", &catalog_groups_flayout, 
              &i_group_file, &N_halos_groups_file, &group_count_infile, group_halos, N_groups_files, &group_count);
    group_count=0; // Reset this after every group read as we are using a dummy 1 element array for group_halos

    if(n_subgroups <= 0)
      SID_log_warning("Just read in a group with n_subgroups=%d... Group ID=%d, Descendent ID=%d", SID_LOG_COMMENT, 
          n_subgroups, group_halos[0].id, group_halos[0].desc_id);

    // Groups with n_subgroups=0 are spurious halos identified by the halo finder which should be ignored...
    if(n_subgroups > 0)
    {
      // The first subhalo is actually the FOF halo but with the substructure
      // removed.  We want to restore this to be just the FOF halo with no
      // alterations.
      read_subgroup(fin_trees, *halos, halo_count);
      read_halo(&fin_subgroup_halos, "data", sim, corrected_snapshot, "subgroups", &catalog_subgroups_flayout, 
          &i_subgroup_file, &N_halos_subgroups_file, &subgroup_count_infile, *halos, N_subgroups_files, &halo_count);
      // Copy the relevant FOF group data over the top...
      memcpy(&((*halos)[halo_count-1].id_MBP), &(group_halos[0].id_MBP), sizeof(Halo)-offsetof(Halo, id_MBP)); 
      (*halos)[halo_count-1].n_subgroups = group_halos[0].n_subgroups-1;
      (*halos)[halo_count-1].type = 0;

      // Deal with any remaining subhalos
      for (int i_subgroup=1; i_subgroup<n_subgroups; i_subgroup++){
        read_subgroup(fin_trees, *halos, halo_count);
        read_halo(&fin_subgroup_halos, "data", sim, corrected_snapshot, "subgroups", &catalog_subgroups_flayout, 
            &i_subgroup_file, &N_halos_subgroups_file, &subgroup_count_infile, *halos, N_subgroups_files, &halo_count);
        (*halos)[halo_count-1].type = 1;
      }
    }
  }
  
  // Close the files
  fclose(fin_group_halos);
  fclose(fin_subgroup_halos);
  fclose(fin_trees);

  if (halo_count!=N_halos){
    fprintf(stderr, "ERROR==> halo_count != N_halos <==ERROR\n");
    ABORT(EXIT_FAILURE);
  }

  SID_log("...done", SID_LOG_CLOSE);

  return header;
}

void free_trees(Halo **halos){
  // Free allocated arrays
  SID_free(SID_FARG *halos);
}

