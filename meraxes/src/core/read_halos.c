#include "meraxes.h"
#include "tree_flags.h"
#include <hdf5.h>
#include <hdf5_hl.h>

int get_corrected_snapshot(run_globals_t *run_globals, int snapshot)
{
  int  total_sim_snaps = run_globals->params.TotalSimSnaps;
  int  n_every_snaps   = run_globals->params.NEverySnap;

  // Calculate the actual (unsampled) simulation snapshot.
  if (n_every_snaps>1)
    return total_sim_snaps - ((int)((total_sim_snaps+1)/n_every_snaps))*n_every_snaps + snapshot*n_every_snaps;
  else
    return snapshot;
}

static inline bool is_ghost(int flags)
{
 if ((flags & TREE_CASE_GHOST)==TREE_CASE_GHOST)
   return true;
 else
   return false;
}

static void halo_catalog_filename(
  char *simulation_dir,
  char *catalog_file_prefix,
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
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
      else if (*i_layout==1)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
      else if (*i_layout==2)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
      else if (*i_layout==3)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);

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
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
  else if (*i_layout==1)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
  else if (*i_layout==2)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
  else if (*i_layout==3)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);

}


static void inline read_catalogs_header(
  FILE *fin,
  int  *i_file,
  int  *n_files,
  int  *n_halos_file,
  int  *n_halos_total )
{
  fread(i_file       , sizeof(int), 1, fin);
  fread(n_files      , sizeof(int), 1, fin);
  fread(n_halos_file , sizeof(int), 1, fin);
  fread(n_halos_total, sizeof(int), 1, fin);
}


static void read_catalog_halos(
    FILE           **fin,
    char           *simulation_dir,
    char           *catalog_file_prefix,
    int             snapshot,
    int            *flayout_switch,
    int            *i_file,
    int            *n_halos_file,
    int            *i_halo,
    catalog_halo_t *halo,
    int             n_to_read)
{

  char                 fname[STRLEN];
  int                  dummy;
  int                  n_from_this_file;

  // Is this the first read?
  if((*fin)==NULL)
  {
    halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, "subgroups", *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin==NULL)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, n_halos_file, &dummy);
  }

  // Have we already read all the halos in this file?
  if((*i_halo)>=(*n_halos_file))
  {
    // SID_log("***i_halo = %d, n_halos_file = %d", SID_LOG_COMMENT, (*i_halo), (*n_halos_file));
    fclose(*fin);
    (*i_file)++;
    (*i_halo) = 0;
    halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, "subgroups", *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin==NULL)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, n_halos_file, &dummy);
  }

  // Read in as many halos as we can from this file
  if((*i_halo + n_to_read) <= *n_halos_file)
  {
    fread(halo, sizeof(catalog_halo_t), n_to_read, *fin);
    *i_halo += n_to_read;
  }
  else
  {
    // read in as many as we can from this file and then get the rest from the next file
    n_from_this_file = (*n_halos_file)- *i_halo;

    // DEBUG
    // SID_log("Spilling over to next file: i_halo = %d, n_to_read = %d,"\
    //     "n_halos_file = %d, n_from_this_file = %d", SID_LOG_COMMENT, *i_halo,
    //     n_to_read, *n_halos_file, n_from_this_file);

    fread(halo, sizeof(catalog_halo_t), n_from_this_file, *fin);
    *i_halo += n_from_this_file;
    n_to_read -= n_from_this_file;
    read_catalog_halos(fin, simulation_dir, catalog_file_prefix, snapshot, flayout_switch, i_file, n_halos_file, i_halo, halo, n_to_read);
  }

}


static void inline convert_input_halo_units(run_globals_t *run_globals, halo_t *halo, int snapshot)
{
  halo->Mvir /= 1.0e10;

  // Update the virial properties
  halo->Mvir = calculate_Mvir(run_globals, halo);
  halo->Rvir = calculate_Rvir(run_globals, halo, halo->Mvir, snapshot);
  halo->Vvir = calculate_Vvir(run_globals, halo->Mvir, halo->Rvir);
}


//! Buffered read of hdf5 trees into halo structures
static void read_trees_and_catalogs(
  run_globals_t *run_globals,
  int            snapshot,
  hid_t          fd,
  halo_t        *halo,
  int            n_halos,
  fof_group_t   *fof_group,
  int            n_fof_groups,
  int            n_requested_forests,
  int           *n_halos_kept,
  int           *n_fof_groups_kept,
  int           *index_lookup)
{
  // I guess this should ideally be equal to the chunk size of the input hdf5 file...
  int buffer_size = 1000;
  int n_read = 0;
  int n_to_read = 0;
  bool keep_flag;

  FILE *fin_catalogs = NULL;
  int flayout_switch = -1;
  int i_catalog_file = 0;
  int n_halos_in_catalog_file = 0;
  int i_halo_in_catalog_file = 0;
  char catalog_file_prefix[50];
  char simulation_dir[STRLEN];
  catalog_halo_t *catalog_buffer;

  sprintf(catalog_file_prefix, "%s", run_globals->params.CatalogFilePrefix);
  sprintf(simulation_dir, "%s", run_globals->params.SimulationDir);

  tree_entry_t *tree_buffer;
  tree_buffer = SID_malloc(sizeof(tree_entry_t) * buffer_size);
  catalog_buffer = SID_malloc(sizeof(catalog_halo_t) * buffer_size);

  *n_halos_kept = 0;
  *n_fof_groups_kept = 0;

  // DEBUG
  // SID_log("Calling read_trees_and_catalogs() with:", SID_LOG_OPEN);
  // SID_log("snapshot = %d", SID_LOG_COMMENT, snapshot);
  // SID_log("n_halos = %d", SID_LOG_COMMENT, n_halos);
  // SID_log("n_fof_groups = %d", SID_LOG_COMMENT, n_fof_groups);
  // SID_log("n_requested_forests = %d", SID_LOG_COMMENT, n_requested_forests);
  // SID_log("---", SID_LOG_CLOSE);

  size_t dst_size = sizeof(tree_entry_t);
  size_t dst_offsets[9] = {
    HOFFSET(tree_entry_t, id),
    HOFFSET(tree_entry_t, flags),
    HOFFSET(tree_entry_t, desc_id),
    HOFFSET(tree_entry_t, tree_id),
    HOFFSET(tree_entry_t, file_offset),
    HOFFSET(tree_entry_t, desc_index),
    HOFFSET(tree_entry_t, central_index),
    HOFFSET(tree_entry_t, forest_id),
    HOFFSET(tree_entry_t, fof_mvir) };
  size_t dst_sizes[9] = {
    sizeof(tree_buffer[0].id),
    sizeof(tree_buffer[0].flags),
    sizeof(tree_buffer[0].desc_id),
    sizeof(tree_buffer[0].tree_id),
    sizeof(tree_buffer[0].file_offset),
    sizeof(tree_buffer[0].desc_index),
    sizeof(tree_buffer[0].central_index),
    sizeof(tree_buffer[0].forest_id),
    sizeof(tree_buffer[0].fof_mvir) };


  keep_flag = true;
  while(n_read<n_halos)
  {

    if((n_halos-n_read) >= buffer_size)
      n_to_read = buffer_size;
    else
      n_to_read = n_halos-n_read;

    // read in a tree_buffer of the trees
    if(SID.My_rank == 0)
      H5TBread_records(fd, "trees", n_read, (hsize_t)n_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);
    SID_Bcast(tree_buffer, n_to_read * sizeof(tree_entry_t), 0, SID.COMM_WORLD);

    // read in the corresponding catalog entrys
    if(SID.My_rank == 0)
      read_catalog_halos(&fin_catalogs, simulation_dir, catalog_file_prefix,
          snapshot, &flayout_switch, &i_catalog_file, &n_halos_in_catalog_file,
          &i_halo_in_catalog_file, catalog_buffer, n_to_read);
    SID_Bcast(catalog_buffer, n_to_read * sizeof(catalog_halo_t), 0, SID.COMM_WORLD);

    // paste the data into the halo structures
    for(int jj=0; jj<n_to_read; jj++)
    {

      if(run_globals->RequestedForestId!=NULL)
      {
        if(bsearch(&(tree_buffer[jj].forest_id), run_globals->RequestedForestId,
              (size_t)n_requested_forests, sizeof(int), compare_ints) != NULL)
          keep_flag = true;
        else
          keep_flag = false;
      }

      if(keep_flag)
      {
        halo[*n_halos_kept].ID = tree_buffer[jj].id;
        halo[*n_halos_kept].TreeFlags = tree_buffer[jj].flags;
        halo[*n_halos_kept].SnapOffset = tree_buffer[jj].file_offset;
        halo[*n_halos_kept].DescIndex = tree_buffer[jj].desc_index;
        halo[*n_halos_kept].Mvir = tree_buffer[jj].fof_mvir;  // this will be overwritten for type>0 halos later
        halo[*n_halos_kept].NextHaloInFOFGroup = NULL;
        halo[*n_halos_kept].ForestID = tree_buffer[jj].forest_id;

        if(index_lookup)
          index_lookup[*n_halos_kept] = n_read+jj;

        if(n_read+jj == tree_buffer[jj].central_index)
        {
          halo[*n_halos_kept].Type = 0;
          fof_group[(*n_fof_groups_kept)++].FirstHalo = &(halo[*n_halos_kept]);
        }
        else
        {
          halo[*n_halos_kept].Type = 1;
          halo[(*n_halos_kept)-1].NextHaloInFOFGroup = &(halo[*n_halos_kept]);
        }

        halo[*n_halos_kept].FOFGroup = &(fof_group[(*n_fof_groups_kept)-1]);

        // paste in the halo properties
        halo[*n_halos_kept].id_MBP             = catalog_buffer[jj].id_MBP;
        halo[*n_halos_kept].Len                = catalog_buffer[jj].n_particles;
        halo[*n_halos_kept].Pos[0]             = catalog_buffer[jj].position_MBP[0];
        halo[*n_halos_kept].Pos[1]             = catalog_buffer[jj].position_MBP[1];
        halo[*n_halos_kept].Pos[2]             = catalog_buffer[jj].position_MBP[2];
        halo[*n_halos_kept].Vel[0]             = catalog_buffer[jj].velocity_COM[0];
        halo[*n_halos_kept].Vel[1]             = catalog_buffer[jj].velocity_COM[1];
        halo[*n_halos_kept].Vel[2]             = catalog_buffer[jj].velocity_COM[2];
        halo[*n_halos_kept].Rvir               = catalog_buffer[jj].R_vir;
        halo[*n_halos_kept].Rhalo              = catalog_buffer[jj].R_halo;
        halo[*n_halos_kept].Rmax               = catalog_buffer[jj].R_max;
        halo[*n_halos_kept].Vmax               = catalog_buffer[jj].V_max;
        halo[*n_halos_kept].VelDisp            = catalog_buffer[jj].sigma_v;
        halo[*n_halos_kept].Spin[0]            = catalog_buffer[jj].spin[0];
        halo[*n_halos_kept].Spin[1]            = catalog_buffer[jj].spin[1];
        halo[*n_halos_kept].Spin[2]            = catalog_buffer[jj].spin[2];
        halo[*n_halos_kept].Galaxy             = NULL;
        if(halo[*n_halos_kept].Type > 0)
        {
          halo[*n_halos_kept].Mvir             = catalog_buffer[jj].M_vir;
        }

        convert_input_halo_units(run_globals, &(halo[*n_halos_kept]), snapshot);

        (*n_halos_kept)++;
      }

    }

    n_read += n_to_read;
  }

  // close the catalogs file
  if(SID.My_rank == 0)
    if(fin_catalogs)
      fclose(fin_catalogs);

  // free the buffers
  SID_free(SID_FARG catalog_buffer);
  SID_free(SID_FARG tree_buffer);

}


static trees_info_t read_trees_info(hid_t fd)
{

  trees_info_t trees_info;

  H5LTget_attribute_int(fd, "trees", "n_step", &(trees_info.n_step));
  H5LTget_attribute_int(fd, "trees", "n_search", &(trees_info.n_search));
  H5LTget_attribute_int(fd, "trees", "n_halos", &(trees_info.n_halos));
  H5LTget_attribute_int(fd, "trees", "n_halos_max", &(trees_info.n_halos_max));
  H5LTget_attribute_int(fd, "trees", "max_tree_id", &(trees_info.max_tree_id));
  H5LTget_attribute_int(fd, "trees", "n_fof_groups", &(trees_info.n_fof_groups));
  H5LTget_attribute_int(fd, "trees", "n_fof_groups_max", &(trees_info.n_fof_groups_max));

  return trees_info;

}


static void select_forests(run_globals_t *run_globals)
{

  char fname[STRLEN];
  hid_t fin;
  int i_forest;
  int i_req;
  int n_forests;
  int n_halos_tot;
  int n_halos_unused;
  int *forest_id;
  int *max_contemp_halo;
  int *max_contemp_fof;
  int *n_halos;
  int max_halos;
  int max_fof_groups;
  int *n_requested_forests = &(run_globals->NRequestedForests);
  int *requested_ind;
  bool sample_forests;

  SID_log("Calling select_forests()...", SID_LOG_COMMENT);

  // are we sampling the forests or just dividing them amoungst cores?
  if(run_globals->RequestedForestId == NULL)
    sample_forests = false;
  else
    sample_forests = true;

  // if this is the master rank then read in the forest info
  if(SID.My_rank == 0)
  {
    sprintf(fname, "%s/trees/forests_info.hdf5", run_globals->params.SimulationDir);
    if ((fin = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(EXIT_FAILURE);
    }

    // find out how many forests there are
    H5LTget_attribute_int(fin, "info", "n_forests", &n_forests);

    // find out how many halos in the whole simulation there are
    H5LTget_attribute_int(fin, "info", "n_subgroups", &n_halos_tot);
    H5LTget_attribute_int(fin, "info", "n_unused", &n_halos_unused);
    n_halos_tot += n_halos_unused;

    // allocate the arrays
    forest_id        = (int *)SID_malloc(sizeof(int) * n_forests);
    max_contemp_halo = (int *)SID_malloc(sizeof(int) * n_forests);
    max_contemp_fof  = (int *)SID_malloc(sizeof(int) * n_forests);
    n_halos          = (int *)SID_malloc(sizeof(int) * n_forests);

    // read in the max number of contemporaneous halos and groups and the forest ids
    H5LTread_dataset_int(fin, "info/max_contemporaneous_halos"     , max_contemp_halo);
    H5LTread_dataset_int(fin, "info/max_contemporaneous_fof_groups", max_contemp_fof);
    H5LTread_dataset_int(fin, "info/forest_id"                     , forest_id);
    H5LTread_dataset_int(fin, "info/n_halos"                       , n_halos);

    // close the file
    H5Fclose(fin);
  }

  // broadcast the forest info
  SID_Bcast(&n_forests, sizeof(int), 0, SID.COMM_WORLD);
  SID_Bcast(&n_halos_tot, sizeof(int), 0, SID.COMM_WORLD);
  SID_Bcast(&n_halos_unused, sizeof(int), 0, SID.COMM_WORLD);
  if(SID.My_rank > 0)
  {
    forest_id        = (int *)SID_malloc(sizeof(int) * n_forests);
    max_contemp_halo = (int *)SID_malloc(sizeof(int) * n_forests);
    max_contemp_fof  = (int *)SID_malloc(sizeof(int) * n_forests);
    n_halos          = (int *)SID_malloc(sizeof(int) * n_forests);
  }
  SID_Bcast(forest_id       , sizeof(int) * n_forests, 0, SID.COMM_WORLD);
  SID_Bcast(max_contemp_halo, sizeof(int) * n_forests, 0, SID.COMM_WORLD);
  SID_Bcast(max_contemp_fof , sizeof(int) * n_forests, 0, SID.COMM_WORLD);
  SID_Bcast(n_halos         , sizeof(int) * n_forests, 0, SID.COMM_WORLD);

  // if we have read in a list of requested forest IDs then use these to create
  // an array of indices pointing to the elements we want
  if(sample_forests)
  {
    requested_ind = (int *)SID_malloc(sizeof(int) * (*n_requested_forests));
    n_halos_tot = 0;
    for(i_forest=0, i_req=0; (i_forest<n_forests) && (i_req<(*n_requested_forests)); i_forest++)
    {
      if(forest_id[i_forest] == run_globals->RequestedForestId[i_req])
      {
        requested_ind[i_req] = i_forest;
        n_halos_tot+= n_halos[i_forest];
        i_req++;
      }
    }
    n_forests = (*n_requested_forests);
  }
  else
  {
    // if we haven't asked for any specific forest IDs then just fill the
    // requested ind array sequentially
    requested_ind = (int *)SID_malloc(sizeof(int) * n_forests);
    for(i_req=0; i_req<n_forests; i_req++)
      requested_ind[i_req] = i_req;
  }

  // if we are running the code with more than one core, let's subselect the forests on each one
  if(SID.n_proc > 1)
  {
    if(n_forests < SID.n_proc)
    {
      SID_log_error("There are fewer processors than there are forests to be processed (%d).  Try again with fewer cores!", n_forests);
      ABORT(EXIT_FAILURE);
    }

    // malloc the arrays we need for keeping track of the load balancing
    int *rank_first_forest = (int *)SID_malloc(sizeof(int) * SID.n_proc);
    int *rank_last_forest  = (int *)SID_malloc(sizeof(int) * SID.n_proc);
    int *rank_n_forests    = (int *)SID_malloc(sizeof(int) * SID.n_proc);
    int *rank_n_halos      = (int *)SID_malloc(sizeof(int) * SID.n_proc);
    int i_rank;
    int n_halos_used;
    int n_halos_target;

    // initialise
    for(int ii=0; ii<SID.n_proc; ii++)
    {
      rank_first_forest[ii] = 0;
      rank_last_forest[ii]  = 0;
      rank_n_forests[ii]    = 0;
      rank_n_halos[ii]      = 0;
    }

    // loop through each rank
    for(i_rank=0, n_halos_used=0, i_forest=0, n_halos_target=0; i_rank < SID.n_proc; i_rank++, i_forest++)
    {
      // start with the next forest (or first forest if this is the first rank)
      rank_first_forest[i_rank] = i_forest;
      rank_last_forest[i_rank]  = i_forest;

      // if we still have forests left to check from the total list of forests
      if(i_forest < n_forests)
      {
        // add this forest's halos to this rank's total
        rank_n_halos[i_rank] += n_halos[requested_ind[i_forest]];

        // adjust our target to smooth out variations as much as possible
        n_halos_target = (n_halos_tot - n_halos_used)/(SID.n_proc - i_rank);

        // increment the counter of the number of forests on this rank
        rank_n_forests[i_rank]++;

        // keep adding forests until we have reached (or exceeded) the current target
        while((rank_n_halos[i_rank] < n_halos_target) && (i_forest < (n_forests-1)))
        {
          i_forest++;
          rank_n_forests[i_rank]++;
          rank_last_forest[i_rank] = i_forest;
          rank_n_halos[i_rank]    += n_halos[requested_ind[i_forest]];
        }

        // updated the total number of used halos
        n_halos_used += rank_n_halos[i_rank];
      }

      // now double back and make sure we haven't got way more halos on this rank than any previous rank
      for(int j_rank=i_rank-1; j_rank > -1; j_rank--)
      {
        while((rank_n_halos[j_rank+1] > 1.05*rank_n_halos[j_rank]) && (rank_n_halos[j_rank+1] > 1))
        {
          rank_n_halos[j_rank+1] -= n_halos[rank_first_forest[j_rank+1]];
          rank_n_halos[j_rank]   += n_halos[rank_first_forest[j_rank+1]];

          rank_first_forest[j_rank+1]++;
          rank_last_forest[j_rank]++;

          rank_n_forests[j_rank+1]--;
          rank_n_forests[j_rank]++;
        }
      }
    }

    // add any uncounted forests to the last process
    for(i_forest++ ;i_forest < n_forests; i_forest++)
    {
      rank_last_forest[SID.n_proc-1] = i_forest;
      rank_n_halos[SID.n_proc-1] += n_halos[requested_ind[i_forest]];
      rank_n_forests[SID.n_proc-1]++;
    }

    // create our list of forest_ids for this rank
    if((SID.My_rank < SID.n_proc-1) || sample_forests)
      *n_requested_forests = rank_n_forests[SID.My_rank];
    else
      *n_requested_forests = rank_n_forests[SID.My_rank]+1;

    if(run_globals->RequestedForestId != NULL)
      run_globals->RequestedForestId = SID_realloc(run_globals->RequestedForestId, *n_requested_forests * sizeof(int));
    else
      run_globals->RequestedForestId = (int *)SID_malloc(sizeof(int) * (*n_requested_forests));

    for(int ii=rank_first_forest[SID.My_rank], jj=0; ii <= rank_last_forest[SID.My_rank]; ii++, jj++)
      run_globals->RequestedForestId[jj] = forest_id[requested_ind[ii]];

    // note that when we actually read in the halos, the last rank will always
    // take any halos that have a forest_id of -1 unless we provided a list of
    // requested forest IDs.  This fact should have been taken in to account
    // above when we load balanced the forests.
    if((SID.My_rank == SID.n_proc-1) && !sample_forests)
      run_globals->RequestedForestId[*n_requested_forests -1] = -1;

    // DEBUG
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("DEBUGGING INFO FOR select_forests()", SID_LOG_COMMENT);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("n_requested_forests = %d", SID_LOG_COMMENT, *n_requested_forests);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("run_globals->RequestedForestId[0] = %d", SID_LOG_COMMENT, run_globals->RequestedForestId[0]);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("run_globals->RequestedForestId[-1] = %d", SID_LOG_COMMENT, run_globals->RequestedForestId[*n_requested_forests-1]);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("rank_n_halos[%d] = %d", SID_LOG_COMMENT, SID.My_rank, rank_n_halos[SID.My_rank]);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("rank_n_forests[%d] = %d", SID_LOG_COMMENT, SID.My_rank, rank_n_forests[SID.My_rank]);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("rank_first_forest = %d", SID_LOG_COMMENT, rank_first_forest[SID.My_rank]);
    // SID_Barrier(SID.COMM_WORLD);
    // SID_log("rank_last_forest = %d", SID_LOG_COMMENT, rank_last_forest[SID.My_rank]);
    // SID_Barrier(SID.COMM_WORLD);


    // free the arrays
    SID_free(SID_FARG rank_n_halos);
    SID_free(SID_FARG rank_n_forests);
    SID_free(SID_FARG rank_last_forest);
    SID_free(SID_FARG rank_first_forest);

  }


  // loop through and tot up the max number of halos and fof_groups we will need to allocate
  if(sample_forests)
  {
    max_halos = 0;
    max_fof_groups = 0;
    for(int i_forest=0, i_req=0; (i_forest<n_forests) && (i_req<(*n_requested_forests)); i_forest++)
    {
      if(forest_id[requested_ind[i_forest]] == run_globals->RequestedForestId[i_req])
      {
        max_halos += max_contemp_halo[requested_ind[i_forest]];
        max_fof_groups += max_contemp_fof[requested_ind[i_forest]];
        i_req++;
      }
    }

    // store the maximum number of halos and fof groups needed at any one snapshot
    run_globals->NHalosMax = max_halos;
    run_globals->NFOFGroupsMax = max_fof_groups;
  }

  // {
  //   // DEBUG
  //   hid_t fd;
  //   char fname[50];
  //   hsize_t dims;
  //   sprintf(fname, "debug_%d.hdf5", SID.My_rank);
  //   fd = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //   dims = (hsize_t)(*n_requested_forests);
  //   H5LTmake_dataset_int(fd, "forest_id", 1, &dims, run_globals->RequestedForestId);
  //   dims = (hsize_t)n_forests;
  //   H5LTmake_dataset_int(fd, "requested_ind", 1, &dims, requested_ind);
  //   H5LTmake_dataset_int(fd, "n_halos", 1, &dims, n_halos);
  //   H5Fclose(fd);
  // }

  // free the arrays
  SID_free(SID_FARG requested_ind);
  SID_free(SID_FARG n_halos);
  SID_free(SID_FARG forest_id);
  SID_free(SID_FARG max_contemp_fof);
  SID_free(SID_FARG max_contemp_halo);

}


trees_info_t read_halos(
  run_globals_t  *run_globals,
  int                  snapshot,
  halo_t        **halo,
  fof_group_t   **fof_group,
  int           **index_lookup)
{

  int             n_halos;                 //!< Number of halos
  char            fname[STRLEN];
  int             corrected_snapshot;
  int             n_fof_groups;
  trees_info_t    trees_info;
  hid_t           fin_trees;

  int             n_halos_kept;
  int             n_fof_groups_kept;
  int             n_requested_forests = run_globals->NRequestedForests;

  SID_log("Reading snapshot %d (z=%.2f) trees and halos...", SID_LOG_OPEN|SID_LOG_TIMER, snapshot, run_globals->ZZ[snapshot]);

  corrected_snapshot = get_corrected_snapshot(run_globals, snapshot);

  // open the tree file
  if(SID.My_rank == 0)
  {
    sprintf(fname, "%s/trees/horizontal_trees_%03d.hdf5", run_globals->params.SimulationDir, corrected_snapshot);
    if ((fin_trees = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(EXIT_FAILURE);
    }

    // read the info attributes
    trees_info = read_trees_info(fin_trees);
  }

  // if necessary, broadcast the tree file info
  SID_Bcast(&trees_info, sizeof(trees_info_t), 0, SID.COMM_WORLD);

  n_halos = trees_info.n_halos;
  n_fof_groups = trees_info.n_fof_groups;

  // If necessary, allocate the halo array
  if(*halo == NULL)
  {
    run_globals->NHalosMax = trees_info.n_halos_max;
    run_globals->NFOFGroupsMax = trees_info.n_fof_groups_max;

    // if required, select forests and calculate the maximum number of halos and fof groups
    if((n_requested_forests > -1) || (SID.n_proc > 1))
    {
      select_forests(run_globals);
      *index_lookup = SID_malloc(sizeof(int) * run_globals->NHalosMax);
    }

    SID_log("Allocating halo array with %d elements...", SID_LOG_COMMENT, run_globals->NHalosMax);
    *halo = SID_malloc(sizeof(halo_t) * run_globals->NHalosMax);
  }

  // Allocate the fof_group array if necessary
  if(*fof_group == NULL)
  {
    SID_log("Allocating fof_group array with %d elements...", SID_LOG_COMMENT, run_globals->NFOFGroupsMax);
    *fof_group = SID_malloc(sizeof(fof_group_t) * run_globals->NFOFGroupsMax);
  }

  // reset the fof group pointers and index lookup (if necessary)
  for(int ii=0; ii<run_globals->NFOFGroupsMax; ii++)
    (*fof_group)[ii].FirstHalo  = NULL;
  if(n_requested_forests > -1)
    for(int ii=0; ii<run_globals->NHalosMax; ii++)
      (*index_lookup)[ii] = -1;

  if (n_halos<1)
  {
    SID_log("No halos in this file... skipping...", SID_LOG_CLOSE);
    if(SID.My_rank == 0)
      H5Fclose(fin_trees);
    return trees_info;
  }

  // read in the trees
  read_trees_and_catalogs(run_globals, snapshot, fin_trees, *halo, n_halos,
      *fof_group, n_fof_groups, n_requested_forests,
      &n_halos_kept, &n_fof_groups_kept, *index_lookup);

  // close the tree file
  if(SID.My_rank == 0)
    H5Fclose(fin_trees);

  // if subsampling the trees, then update the trees_info to reflect what we now have
  if(n_requested_forests > -1)
  {
    trees_info.n_halos = n_halos_kept;
    trees_info.n_fof_groups = n_fof_groups_kept;
  }

  SID_log("Read %d halos in %d fof_groups.", SID_LOG_COMMENT, trees_info.n_halos, trees_info.n_fof_groups);

  if(SID.n_proc > 0)
  {
    SID_Allreduce(SID_IN_PLACE, &n_halos_kept, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
    SID_Allreduce(SID_IN_PLACE, &n_fof_groups_kept, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
    SID_log("Read %d halos in %d fof_groups in total.", SID_LOG_COMMENT, n_halos_kept, n_fof_groups_kept);
  }

  SID_log("...done", SID_LOG_CLOSE);

  return trees_info;
}


