#include "meraxes.h"
#include "tree_flags.h"

static inline bool is_ghost(int flags)
{
 if ((flags & TREE_CASE_GHOST)==TREE_CASE_GHOST)
   return true;
 else
   return false;
}

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

static void ghost_catalog_filename(
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
        sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties/ghosts_%03d.catalog_%s_properties.%d", root, sim, sim, snapshot, group_type, snapshot, group_type, sub);
      else if (*i_layout==1)
        sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties/ghosts_%03d.catalog_%s_properties", root, sim, sim, snapshot, group_type, snapshot, group_type);
      else if (*i_layout==2)
        sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties.%d", root, sim, sim, snapshot, group_type, sub);
      else if (*i_layout==3)
        sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties", root, sim, sim, snapshot, group_type);
      
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
    sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties/ghosts_%03d.catalog_%s_properties.%d", root, sim, sim, snapshot, group_type, snapshot, group_type, sub);
  else if (*i_layout==1)
    sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties/ghosts_%03d.catalog_%s_properties", root, sim, sim, snapshot, group_type, snapshot, group_type);
  else if (*i_layout==2)
    sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties.%d", root, sim, sim, snapshot, group_type, sub);
  else if (*i_layout==3)
    sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/ghost_catalogs/ghosts_%03d.catalog_%s_properties", root, sim, sim, snapshot, group_type);


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
  // cur_model_halo->ID        = *halo_count;
  cur_model_halo->Mvir               = halo_in.M_vir;
  cur_model_halo->Len                = halo_in.n_particles;
  cur_model_halo->Pos[0]             = halo_in.position_MBP[0];
  cur_model_halo->Pos[1]             = halo_in.position_MBP[1];
  cur_model_halo->Pos[2]             = halo_in.position_MBP[2];
  cur_model_halo->Vel[0]             = halo_in.velocity_COM[0];
  cur_model_halo->Vel[1]             = halo_in.velocity_COM[1];
  cur_model_halo->Vel[2]             = halo_in.velocity_COM[2];
  cur_model_halo->Rvir               = halo_in.R_vir;
  cur_model_halo->Rhalo              = halo_in.R_halo;
  cur_model_halo->Rmax               = halo_in.R_max;
  cur_model_halo->Vmax               = halo_in.V_max;
  cur_model_halo->VelDisp            = halo_in.sigma_v;
  cur_model_halo->Spin[0]            = halo_in.spin[0];
  cur_model_halo->Spin[1]            = halo_in.spin[1];
  cur_model_halo->Spin[2]            = halo_in.spin[2];
  cur_model_halo->NextHaloInFOFGroup = NULL;
  cur_model_halo->Galaxy             = NULL;

  // Update the counters
  (*halo_count)++;
  (*i_halo)++;
}


static void inline read_trees_header(FILE *fin, trees_header_struct *header)
{
  fread(&(header->n_step)              , sizeof(int), 1, fin);
  fread(&(header->n_search)            , sizeof(int), 1, fin);
  fread(&(header->n_groups)            , sizeof(int), 1, fin);
  fread(&(header->n_subgroups)         , sizeof(int), 1, fin);
  fread(&(header->n_groups_max)        , sizeof(int), 1, fin);
  fread(&(header->n_subgroups_max)     , sizeof(int), 1, fin);
  fread(&(header->max_tree_id_subgroup), sizeof(int), 1, fin);
  fread(&(header->max_tree_id_group)   , sizeof(int), 1, fin);
}

static void inline read_group(FILE *fin, halo_struct *halo, int i_halo)
{
  int dummy;
  fread(&(halo[i_halo].ID)         , sizeof(int), 1, fin);
  fread(&(halo[i_halo].TreeFlags)  , sizeof(int), 1, fin);
  fread(&dummy                     , sizeof(int), 1, fin);
  fread(&dummy                     , sizeof(int), 1, fin);
  fread(&(halo[i_halo].DescIndex)  , sizeof(int), 1, fin);
  fread(&(halo[i_halo].NSubgroups) , sizeof(int), 1, fin);
}

static void inline read_subgroup(FILE *fin, halo_struct *halo, int i_halo)
{
  int dummy;
  fread(&(halo[i_halo].ID)         , sizeof(int), 1, fin);
  fread(&(halo[i_halo].TreeFlags)  , sizeof(int), 1, fin);
  fread(&dummy                     , sizeof(int), 1, fin);
  fread(&dummy                     , sizeof(int), 1, fin);
  fread(&(halo[i_halo].DescIndex)  , sizeof(int), 1, fin);
  halo[i_halo].NSubgroups = -1;
}


static void inline convert_input_halo_units(halo_struct *halo)
{
  halo->Mvir /= 1.0e10;
}

trees_header_struct read_halos(
  run_globals_struct  *run_globals,
  int                  snapshot,   
  halo_struct        **halo,
  fof_group_struct   **fof_group)      
{

  int                  N_halos;                               //!< Number of halos including ghosts
  int                  N_groups;                              //!< Number of FOF groups including ghosts
  int                  N_halos_groups;                        //!< Number of group halos excluding ghosts
  int                  N_ghosts_groups;                       //!< Number of group ghost halos
  int                  N_halos_subgroups;                     //!< Number of subgroup halos excluding ghosts
  int                  N_ghosts_subgroups;                    //!< Number of ghost subgroup halos
  int                  N_groups_files;                        //!< Number of (real) group files
  int                  N_subgroups_files;                     //!< Number of (real) subgroup files
  int                  N_ghosts_groups_files;                 //!< Number of ghost group files
  int                  N_ghosts_subgroups_files;              //!< Number of ghost subgroup files
  int                  dummy;
  trees_header_struct  header;
  char                 fname[STRLEN];
  FILE                *fin;
  FILE                *fin_trees;
  int                  catalog_groups_flayout           = -1;
  int                  catalog_subgroups_flayout        = -1;
  int                  catalog_ghosts_groups_flayout    = -1;
  int                  catalog_ghosts_subgroups_flayout = -1;
  int                  central_index;                         //!< Index of the central halo associated with this one
  int                  corrected_snapshot;
  char                 sim_variant[18];

  int  total_sim_snaps = run_globals->params.LastSnapShotNr +1;
  int  n_every_snaps   = run_globals->params.NEverySnap;
  int  n_scan_snaps    = run_globals->params.NScanSnap;

  sprintf(sim_variant, "step_%03d_scan_%03d", n_every_snaps, n_scan_snaps);
  SID_log("Reading snapshot %d (%s:%s) trees and halos...", SID_LOG_OPEN|SID_LOG_TIMER, snapshot, run_globals->params.SimName, sim_variant);

  // Calculate the actual (unsampled) simulation snapshot.
  if (n_every_snaps>1)
    corrected_snapshot = total_sim_snaps-1 - ((int)(total_sim_snaps/n_every_snaps))*n_every_snaps + snapshot*n_every_snaps;
  else
    corrected_snapshot = snapshot;

  // HALOS
  // Read the header info
  halo_catalog_filename(run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "groups", 0, &catalog_groups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }
  read_catalogs_header(fin, &dummy, &N_groups_files, &dummy, &N_halos_groups);
  SID_log("N_groups_files = %d", SID_LOG_COMMENT, N_groups_files);
  fclose(fin);

  ghost_catalog_filename(run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "groups", 0, &catalog_ghosts_groups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }
  read_catalogs_header(fin, &dummy, &N_ghosts_groups_files, &dummy, &N_ghosts_groups);
  SID_log("N_ghosts_groups_files = %d", SID_LOG_COMMENT, N_ghosts_groups_files);
  fclose(fin);

  halo_catalog_filename(run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "subgroups", 0, &catalog_subgroups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }
  read_catalogs_header(fin, &dummy, &N_subgroups_files, &dummy, &N_halos_subgroups);
  SID_log("N_subgroups_files = %d", SID_LOG_COMMENT, N_subgroups_files);
  fclose(fin);

  ghost_catalog_filename(run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "subgroups", 0, &catalog_ghosts_subgroups_flayout, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }
  read_catalogs_header(fin, &dummy, &N_ghosts_subgroups_files, &dummy, &N_ghosts_subgroups);
  SID_log("N_ghosts_subgroups_files = %d", SID_LOG_COMMENT, N_ghosts_subgroups_files);
  fclose(fin);

  // TREES
  SID_log("Reading in trees...", SID_LOG_COMMENT);
  // sprintf(fname, "%s/%s/trees/horizontal/trees/%s_%s.trees_horizontal_%d",
      // run_globals->params.SimulationDir, run_globals->params.SimName, run_globals->params.SimName, sim_variant, corrected_snapshot);
  sprintf(fname, "%s/%s/trees/%s_test_ghosts/horizontal/trees/horizontal_trees_ghosts_%03d.dat",
      run_globals->params.SimulationDir, run_globals->params.SimName, run_globals->params.SimName, corrected_snapshot);
  fin_trees = fopen(fname, "rb");
  if (fin_trees==NULL)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }

  // Read the header info
  read_trees_header(fin_trees, &header);

  // // Temp fix for faulty input catalogs
  // if (header.n_subgroups>(N_halos_subgroups+N_ghosts_subgroups))
  //   N_halos_subgroups=header.n_subgroups;
  // if (header.n_groups>(N_halos_groups+N_ghosts_groups))
  //   N_halos_groups = header.n_groups;
  
  N_halos = N_halos_subgroups+N_ghosts_subgroups;
  N_groups = N_halos_groups+N_ghosts_groups;

  if(header.n_subgroups != N_halos)
  {
    SID_log_error("Trees header n_subgroups (%d) != catalog value (%d)", header.n_subgroups, N_halos);
    ABORT(EXIT_FAILURE);
  }
  if(header.n_groups != N_groups)
  {
    SID_log_error("Trees header n_groups (%d) != catalog value (%d)", header.n_groups, N_groups);
    ABORT(EXIT_FAILURE);
  }

  SID_log("N_halos_groups = %d"    , SID_LOG_COMMENT, N_halos_groups);
  SID_log("N_ghosts_groups = %d"   , SID_LOG_COMMENT, N_ghosts_groups);
  SID_log("N_halos_subgroups = %d" , SID_LOG_COMMENT, N_halos_subgroups);
  SID_log("N_ghosts_subgroups = %d", SID_LOG_COMMENT, N_ghosts_subgroups);
  SID_log("N_halos = %d"           , SID_LOG_COMMENT, N_halos);

  // Allocate the halo array
  if (N_halos>0)
    *halo = SID_malloc(sizeof(halo_struct) * N_halos);

  // Allocate the fof_group array
  if (N_groups>0)
  {
    *fof_group = SID_malloc(sizeof(fof_group_struct) * N_groups);
    for(int ii=0; ii<N_groups; ii++)
    {
      (*fof_group)[ii].FirstHalo  = NULL;
    }
  }

  if (N_halos<1)
  {
    SID_log("No halos in this file... skipping...", SID_LOG_CLOSE);
    fclose(fin_trees);
    return header;
  }

  // Initialise everything we need to read the halos in
  FILE *fin_group_halos              = NULL;
  FILE *fin_subgroup_halos           = NULL;
  FILE *fin_group_ghosts              = NULL;
  FILE *fin_subgroup_ghosts           = NULL;
  int   i_group_file                 = 0;
  int   i_subgroup_file              = 0;
  int   i_ghosts_group_file          = 0;
  int   i_ghosts_subgroup_file       = 0;
  int   halo_count                   = 0;
  int   i_halo                       = 0;
  int   N_halos_groups_file          = 0;
  int   N_halos_subgroups_file       = 0;
  int   N_ghosts_groups_file         = 0;
  int   N_ghosts_subgroups_file      = 0;
  int   group_count_infile           = 0;
  int   subgroup_count_infile        = 0;
  int   ghosts_group_count_infile    = 0;
  int   ghosts_subgroup_count_infile = 0;
  int   n_subgroups                  = 0;
  int   group_count                  = 0;
  int   phantom_group_count          = 0;

  // Loop through the groups and subgroups and read them in
  halo_struct group_halos[1];
  for (int i_group=0; i_group<header.n_groups; i_group++){
    read_group(fin_trees, group_halos, group_count);
    n_subgroups = group_halos[group_count].NSubgroups;
    if(is_ghost(group_halos[group_count].TreeFlags))
      read_catalog_halo(&fin_group_ghosts, run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "groups", &catalog_ghosts_groups_flayout, 
          &i_ghosts_group_file, &N_ghosts_groups_file, &ghosts_group_count_infile, group_halos, N_ghosts_groups_files, &group_count);
    else
      read_catalog_halo(&fin_group_halos, run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "groups", &catalog_groups_flayout, 
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
      if(is_ghost((*halo)[halo_count].TreeFlags))
        read_catalog_halo(&fin_subgroup_ghosts, run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "subgroups", &catalog_ghosts_subgroups_flayout, 
            &i_ghosts_subgroup_file, &N_ghosts_subgroups_file, &ghosts_subgroup_count_infile, *halo, N_ghosts_subgroups_files, &halo_count);
      else
        read_catalog_halo(&fin_subgroup_halos, run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "subgroups", &catalog_subgroups_flayout, 
            &i_subgroup_file, &N_halos_subgroups_file, &subgroup_count_infile, *halo, N_subgroups_files, &halo_count);
      i_halo = halo_count-1;
      // Copy the relevant FOF group data over the top...
      memcpy(&((*halo)[i_halo].Mvir), &(group_halos[0].Mvir), sizeof(halo_struct)-offsetof(halo_struct, Mvir)); 
      (*halo)[i_halo].NSubgroups = group_halos[0].NSubgroups-1;
      (*halo)[i_halo].Type = 0;
      convert_input_halo_units(&((*halo)[i_halo]));
      central_index = i_halo;
      (*fof_group)[i_group-phantom_group_count].FirstHalo = &((*halo)[central_index]);
      (*halo)[i_halo].FOFGroup = &((*fof_group)[i_group-phantom_group_count]);

      // DEBUG
      if((*halo)[i_halo].Mvir <= 1.e-9)
      {
        SID_log("Just read a type 0 subhalo with Mvir=%.2e", SID_LOG_OPEN, (*halo)[i_halo].Mvir);
        SID_log("i_group = %d", SID_LOG_COMMENT, i_group);
        SID_log("i_subgroup = %d", SID_LOG_COMMENT, i_halo);
        SID_log("DescIndex = %d", SID_LOG_COMMENT, (*halo)[i_halo].DescIndex);
        SID_log("TreeFlags = %d", SID_LOG_COMMENT, (*halo)[i_halo].TreeFlags);
        SID_log("... need to check this out!", SID_LOG_CLOSE);
      }

      // Deal with any remaining subhalos
      for (int i_subgroup=1; i_subgroup<n_subgroups; i_subgroup++){
        read_subgroup(fin_trees, *halo, halo_count);
        if(is_ghost((*halo)[halo_count].TreeFlags))
          read_catalog_halo(&fin_subgroup_ghosts, run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "subgroups", &catalog_ghosts_subgroups_flayout, 
              &i_ghosts_subgroup_file, &N_ghosts_subgroups_file, &ghosts_subgroup_count_infile, *halo, N_ghosts_subgroups_files, &halo_count);
        else
          read_catalog_halo(&fin_subgroup_halos, run_globals->params.SimulationDir, run_globals->params.SimName, corrected_snapshot, "subgroups", &catalog_subgroups_flayout, 
              &i_subgroup_file, &N_halos_subgroups_file, &subgroup_count_infile, *halo, N_subgroups_files, &halo_count);
        i_halo = halo_count-1;
        (*halo)[i_halo].Type = 1;
        convert_input_halo_units(&((*halo)[i_halo]));
        (*halo)[i_halo].FOFGroup = &((*fof_group)[i_group-phantom_group_count]);
        (*halo)[halo_count-2].NextHaloInFOFGroup = &((*halo)[i_halo]);
      }
    }
  }

  // Close the files
  fclose(fin_group_halos);
  fclose(fin_group_ghosts);
  fclose(fin_subgroup_halos);
  fclose(fin_subgroup_ghosts);
  fclose(fin_trees);

  if (halo_count!=N_halos){
    SID_log_error("halo_count != N_halos\n");
    ABORT(EXIT_FAILURE);
  }

  // Update the header n_groups to take into account phantoms
  header.n_groups -= phantom_group_count;

  SID_log("...done", SID_LOG_CLOSE);

  return header;
}

void free_halos(halo_struct **halo){
  // Free allocated arrays
  SID_free(SID_FARG *halo);
}

