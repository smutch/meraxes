#include "meraxes.h"

static void halo_catalog_filename(
  char *root,      
  char *sim,       
  int   snapshot,  
  char *group_type,
  int   sub,       
  int  *i_layout,  
  char *fname)     
{

  bool flag_success = false;
  FILE *fin;

  // if we need to determine the filename structure...
  if (*i_layout==-1)
  {
    for (*i_layout=0; (*i_layout<4) && (flag_success==false); (*i_layout)++)
    {
      if (*i_layout==0)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, snapshot, group_type, sub);
      else if (*i_layout==1)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type, snapshot, group_type);
      else if (*i_layout==2)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, sub);
      else if (*i_layout==3)
        sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type);
      
      if ((fin = fopen(fname, "rb"))!=NULL)
      {
        flag_success = true;
        fclose(fin);
        break;
      }
    }
  } 

  // ensure we have a valid i_layout value.
  if (*i_layout<0 && *i_layout>3)
  {
    fprintf(stderr, "cannot resolve catalogue filename.\n");
    ABORT(EXIT_FAILURE);
  }
   
  // provide the correct filename
  if (*i_layout==0)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, snapshot, group_type, sub);
  else if (*i_layout==1)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type, snapshot, group_type);
  else if (*i_layout==2)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties.%d", root, sim, snapshot, group_type, sub);
  else if (*i_layout==3)
    sprintf(fname, "%s/%s/catalogs/subfind_%03d.catalog_%s_properties", root, sim, snapshot, group_type);

}


static void inline read_catalogs_header(
  FILE *fin,          
  int  *i_file,       
  int  *N_files,      
  int  *N_halos_file, 
  int  *N_halos_total )
{
  fread(i_file       , sizeof(int), 1, fin);
  fread(N_files      , sizeof(int), 1, fin);
  fread(N_halos_file , sizeof(int), 1, fin);
  fread(N_halos_total, sizeof(int), 1, fin);
}


static void inline read_catalog_halo(
  FILE        **fin,           
  char         *root,          
  char         *sim,           
  int           snapshot,      
  char         *group_type,    
  int          *flayout_switch,
  int          *i_file,        
  int          *N_halos_file,  
  int          *i_halo,        
  halo_struct  *halo,         
  int           N_files,       
  int          *halo_count)    
{

  char                 fname[STRLEN];
  int                  dummy;
  catalog_halo_struct  halo_in;
  halo_struct         *cur_model_halo;

  // Is this the first read?
  if((*fin)==NULL)
  {
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
  if((*i_halo)>=(*N_halos_file))
  {
    // SID_log("***i_halo = %d, N_halos_file = %d", SID_LOG_COMMENT, (*i_halo), (*N_halos_file));
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
  fread(&halo_in, sizeof(catalog_halo_struct), 1, *fin);
  cur_model_halo = &(halo[*halo_count]);

  // Copy over the properties we want to keep
  cur_model_halo->id          = *halo_count;
  cur_model_halo->Mvir        = halo_in.M_vir;
  cur_model_halo->len         = halo_in.n_particles;
  cur_model_halo->position[0] = halo_in.position_MBP[0];
  cur_model_halo->position[1] = halo_in.position_MBP[1];
  cur_model_halo->position[2] = halo_in.position_MBP[2];
  cur_model_halo->velocity[0] = halo_in.velocity_COM[0];
  cur_model_halo->velocity[1] = halo_in.velocity_COM[1];
  cur_model_halo->velocity[2] = halo_in.velocity_COM[2];
  cur_model_halo->Rvir        = halo_in.R_vir;
  cur_model_halo->Rhalo       = halo_in.R_halo;
  cur_model_halo->Rmax        = halo_in.R_max;
  cur_model_halo->Vmax        = halo_in.V_max;
  cur_model_halo->VelDisp     = halo_in.sigma_v;
  cur_model_halo->spin[0]     = halo_in.spin[0];
  cur_model_halo->spin[1]     = halo_in.spin[1];
  cur_model_halo->spin[2]     = halo_in.spin[2];

  
  // Update the halo counters
  (*halo_count)++;
  (*i_halo)++;
}


static void inline read_trees_header(FILE *fin, trees_header_struct *header)
{
  fread(&(header->n_groups)        , sizeof(int), 1, fin);
  fread(&(header->n_subgroups)     , sizeof(int), 1, fin);
  fread(&(header->n_halos_max)     , sizeof(int), 1, fin);
  fread(&(header->n_trees_subgroup), sizeof(int), 1, fin);
  fread(&(header->n_trees_group)   , sizeof(int), 1, fin);
}

static void inline read_group(FILE *fin, halo_struct *halo, int i_halo)
{
  int dummy;
  fread(&(halo[i_halo].id)         , sizeof(int), 1, fin);
  fread(&(halo[i_halo].tree_flags) , sizeof(int), 1, fin);
  fread(&(halo[i_halo].desc_id)    , sizeof(int), 1, fin);
  fread(&dummy                     , sizeof(int), 1, fin);
  fread(&(halo[i_halo].file_offset), sizeof(int), 1, fin);
  fread(&(halo[i_halo].file_index) , sizeof(int), 1, fin);
  fread(&(halo[i_halo].n_subgroups), sizeof(int), 1, fin);
}

static void inline read_subgroup(FILE *fin, halo_struct *halo, int i_halo)
{
  int dummy;
  fread(&(halo[i_halo].id)         , sizeof(int), 1, fin);
  fread(&(halo[i_halo].tree_flags) , sizeof(int), 1, fin);
  fread(&(halo[i_halo].desc_id)    , sizeof(int), 1, fin);
  fread(&dummy                     , sizeof(int), 1, fin);
  fread(&(halo[i_halo].file_offset), sizeof(int), 1, fin);
  fread(&(halo[i_halo].file_index) , sizeof(int), 1, fin);
  halo[i_halo].n_subgroups = -1;
}

trees_header_struct read_halos(
  run_globals_struct  *run_globals,
  int                  snapshot,   
  halo_struct        **halo)      
{

  char                sim[STRLEN];
  int                 N_halos;
  int                 N_halos_groups;
  int                 N_halos_subgroups;
  int                 N_groups_files;
  int                 N_subgroups_files;
  int                 dummy;
  trees_header_struct header;

  int  total_sim_snaps = run_globals->params.LastSnapShotNr +1;
  int  n_every_snaps   = run_globals->params.NEverySnap;
  int  n_scan_snaps    = run_globals->params.NScanSnap;

  strcpy(sim, run_globals->params.SimName);

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
  char  fname[STRLEN];
  FILE *fin;
  FILE *fin_trees;
  int   catalog_groups_flayout    = -1;
  int   catalog_subgroups_flayout = -1;
  halo_catalog_filename("data", sim, corrected_snapshot, "groups", 0, &catalog_groups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(34494);
  }
  read_catalogs_header(fin, &dummy, &N_groups_files, &dummy, &N_halos_groups);
  SID_log("N_groups_files = %d", SID_LOG_COMMENT, N_groups_files);
  fclose(fin);

  halo_catalog_filename("data", sim, corrected_snapshot, "subgroups", 0, &catalog_subgroups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(34494);
  }
  read_catalogs_header(fin, &dummy, &N_subgroups_files, &dummy, &N_halos_subgroups);
  SID_log("N_subgroups_files = %d", SID_LOG_COMMENT, N_subgroups_files);
  fclose(fin);

  // Allocate the halo array
  N_halos = N_halos_subgroups;
  SID_log("N_halos_groups = %d", SID_LOG_COMMENT, N_halos_groups);
  SID_log("N_halos_subgroups = %d", SID_LOG_COMMENT, N_halos_subgroups);
  SID_log("N_halos = %d", SID_LOG_COMMENT, N_halos);
  if (N_halos>0)
    *halo = malloc(sizeof(halo_struct) * N_halos);

  // TREES
  SID_log("Reading in trees...", SID_LOG_COMMENT);
  sprintf(fname, "data/%s/trees/%s_%s/horizontal/trees/%s_%s.trees_horizontal_%d", sim, sim, sim_variant, sim, sim_variant, corrected_snapshot);
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
    SID_log("No halos in this file... skipping...", SID_LOG_CLOSE);
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
  int   phantom_group_count    = 0;

  // Loop through the groups and subgroups and read them in
  halo_struct group_halos[1];
  for (int i_group=0; i_group<header.n_groups; i_group++){
    read_group(fin_trees, group_halos, group_count);
    n_subgroups = group_halos[group_count].n_subgroups;
    read_catalog_halo(&fin_group_halos, "data", sim, corrected_snapshot, "groups", &catalog_groups_flayout, 
        &i_group_file, &N_halos_groups_file, &group_count_infile, group_halos, N_groups_files, &group_count);
    group_count=0; // Reset this after every group read as we are using a dummy 1 element array for group_halos

    if(n_subgroups <= 0)
      phantom_group_count++;

    // Groups with n_subgroups=0 are spurious halos identified by the halo finder which should be ignored...
    if(n_subgroups > 0)
    {
      // The first subhalo is actually the FOF halo but with the substructure
      // removed.  We want to restore this to be just the FOF halo with no
      // alterations.
      read_subgroup(fin_trees, *halo, halo_count);
      read_catalog_halo(&fin_subgroup_halos, "data", sim, corrected_snapshot, "subgroups", &catalog_subgroups_flayout, 
          &i_subgroup_file, &N_halos_subgroups_file, &subgroup_count_infile, *halo, N_subgroups_files, &halo_count);
      // Copy the relevant FOF group data over the top...
      memcpy(&((*halo)[halo_count-1].Mvir), &(group_halos[0].Mvir), sizeof(halo_struct)-offsetof(halo_struct, Mvir)); 
      (*halo)[halo_count-1].n_subgroups = group_halos[0].n_subgroups-1;
      (*halo)[halo_count-1].type = 0;

      // Deal with any remaining subhalos
      for (int i_subgroup=1; i_subgroup<n_subgroups; i_subgroup++){
        read_subgroup(fin_trees, *halo, halo_count);
        read_catalog_halo(&fin_subgroup_halos, "data", sim, corrected_snapshot, "subgroups", &catalog_subgroups_flayout, 
            &i_subgroup_file, &N_halos_subgroups_file, &subgroup_count_infile, *halo, N_subgroups_files, &halo_count);
        (*halo)[halo_count-1].type = 1;
      }
    }
  }


  // SID_log_warning("Skipped %d phantom groups...", SID_LOG_COMMENT, phantom_group_count);
  
  // Close the files
  fclose(fin_group_halos);
  fclose(fin_subgroup_halos);
  fclose(fin_trees);

  if (halo_count!=N_halos){
    SID_log_error("halo_count != N_halos\n");
    ABORT(EXIT_FAILURE);
  }

  SID_log("...done", SID_LOG_CLOSE);

  return header;
}

void free_halos(halo_struct **halo){
  // Free allocated arrays
  SID_free(SID_FARG *halo);
}

