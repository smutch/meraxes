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
  file *fin;

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
      
      if ((fin = fopen(fname, "rb"))!=null)
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
    abort(exit_failure);
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
  FILE **fin,           
  char  *root,          
  char  *sim,           
  int    snapshot,      
  char  *group_type,    
  int   *flayout_switch,
  int   *i_file,        
  int   *N_halos_file,  
  int   *i_halo,        
  Halo  *halos,         
  int    N_files,       
  int   *halo_count )    
{

  char            fname[STR_LEN];
  int             dummy;
  raw_halo_struct halo_in;
  halo_struct    *cur_model_halo;

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
  fread(&halo_in, sizeof(RawHalo), 1, *fin);
  cur_model_halo = &(halos[*halo_count]);

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
  cur_model_halo->sigma_v     = halo_in.VelDisp;
  cur_model_halo->spin[0]     = halo_in.spin[0];
  cur_model_halo->spin[1]     = halo_in.spin[1];
  cur_model_halo->spin[2]     = halo_in.spin[2];

  
  // Update the halo counters
  (*halo_count)++;
  (*i_halo)++;
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
  fread(&(halos[i_halo].tree_flags) , sizeof(int), 1, fin);
  fread(&(halos[i_halo].desc_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].tree_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].file_offset), sizeof(int), 1, fin);
  fread(&(halos[i_halo].file_index) , sizeof(int), 1, fin);
  fread(&(halos[i_halo].n_subgroups), sizeof(int), 1, fin);
}

static void inline read_subgroup(FILE *fin, Halo *halos, int i_halo)
{
  fread(&(halos[i_halo].id)         , sizeof(int), 1, fin);
  fread(&(halos[i_halo].tree_flags) , sizeof(int), 1, fin);
  fread(&(halos[i_halo].desc_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].tree_id)    , sizeof(int), 1, fin);
  fread(&(halos[i_halo].file_offset), sizeof(int), 1, fin);
  fread(&(halos[i_halo].file_index) , sizeof(int), 1, fin);
  halos[i_halo].n_subgroups = -1;
}


