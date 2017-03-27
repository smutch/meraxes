#include "meraxes.h"
#include "tree_flags.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_sort_int.h>
#include <assert.h>
#include <math.h>

static void halo_catalog_filename(
  char *simulation_dir,
  char *catalog_file_prefix,
  int   snapshot,
  char *group_type,
  int   sub,
  int  *i_layout,
  char *fname)
{
  bool  flag_success = false;
  FILE *fin;

  // if we need to determine the filename structure...
  if (*i_layout == -1)
    for (*i_layout = 0; (*i_layout < 4) && (flag_success == false); (*i_layout)++)
    {
      if (*i_layout == 0)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
      else if (*i_layout == 1)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
      else if (*i_layout == 2)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
      else if (*i_layout == 3)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);

      if ((fin = fopen(fname, "rb")) != NULL)
      {
        flag_success = true;
        fclose(fin);
        break;
      }
    }

  // ensure we have a valid i_layout value.
  if (*i_layout < 0 || *i_layout > 3)
  {
    fprintf(stderr, "cannot resolve catalogue filename.\n");
    ABORT(EXIT_FAILURE);
  }

  // provide the correct filename
  if (*i_layout == 0)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
  else if (*i_layout == 1)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
  else if (*i_layout == 2)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
  else if (*i_layout == 3)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);
}


static void inline read_catalogs_header(
  FILE *fin,
  int  *i_file,
  int  *n_files,
  int  *n_halos_file,
  int  *n_halos_total)
{
  fread(i_file, sizeof(int), 1, fin);
  fread(n_files, sizeof(int), 1, fin);
  fread(n_halos_file, sizeof(int), 1, fin);
  fread(n_halos_total, sizeof(int), 1, fin);
}


static void read_catalog_halos(
  FILE          **fin,
  char           *simulation_dir,
  char           *catalog_file_prefix,
  int             snapshot,
  int            *flayout_switch,
  int            *i_file,
  int            *n_halos_file,
  int            *i_halo_in_file,
  int            *i_halo,
  catalog_halo_t *halo,
  int             n_to_read,
  int             type_flag)
{
  char fname[STRLEN];
  char halo_type[10];
  int  dummy;
  int  n_from_this_file;

  switch (type_flag)
  {
    case 0:
      sprintf(halo_type, "groups");
      break;

    case 1:
      sprintf(halo_type, "subgroups");
      break;

    default:
      mlog_error("Unrecognised type_flag in read_catalog_halos()");
      break;
  }

  // Is this the first read?
  if ((*fin) == NULL)
  {
    halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, halo_type, *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin == NULL)
    {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, n_halos_file, &dummy);
  }

  // Have we already read all the halos in this file?
  if ((*i_halo_in_file) >= (*n_halos_file))
  {
    fclose(*fin);
    (*i_file)++;
    (*i_halo_in_file) = 0;
    halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, halo_type, *i_file, flayout_switch, fname);
    *fin              = fopen(fname, "rb");
    if (*fin == NULL)
    {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, n_halos_file, &dummy);
  }

  // Read in as many halos as we can from this file
  if ((*i_halo_in_file + n_to_read) <= *n_halos_file)
  {
    fread(&(halo[*i_halo]), sizeof(catalog_halo_t), n_to_read, *fin);
    *i_halo         += n_to_read;
    *i_halo_in_file += n_to_read;
  }
  else
  {
    // read in as many as we can from this file and then get the rest from the next file
    n_from_this_file = (*n_halos_file) - *i_halo_in_file;

    fread(&(halo[*i_halo]), sizeof(catalog_halo_t), n_from_this_file, *fin);
    *i_halo         += n_from_this_file;
    *i_halo_in_file += n_from_this_file;
    n_to_read       -= n_from_this_file;
    read_catalog_halos(fin, simulation_dir, catalog_file_prefix, snapshot, flayout_switch, i_file, n_halos_file, i_halo_in_file, i_halo, halo, n_to_read, type_flag);
  }
}


static void inline convert_input_virial_props(double *Mvir, double *Rvir, double *Vvir, double *FOFMvirModifier, int len, int snapshot)
{
  if (len >= 0){
    // Update the virial properties for subhalos
    *Mvir = calculate_Mvir(*Mvir, len);         
    *Rvir = calculate_Rvir(*Mvir, snapshot);    
  }
  else{
    // Convert the mass unit for FoFs
    *Mvir /= 1.0e10;
    if (run_globals.RequestedMassRatioModifier == 1){
      // Modifier the FoF mass and update the virial radius
      *FOFMvirModifier = interpolate_modifier(run_globals.mass_ratio_modifier, log10(*Mvir / run_globals.params.Hubble_h)+10.0);
      *Mvir *= *FOFMvirModifier;
      *Rvir  = calculate_Rvir(*Mvir, snapshot);
    }
  }
  *Vvir  = calculate_Vvir(*Mvir, *Rvir);
}


//! Buffered read of hdf5 trees into halo structures
static void read_trees_and_catalogs(
  int          snapshot,
  int          unsampled_snapshot,
  hid_t        fd,
  halo_t      *halo,
  int          n_halos,
  fof_group_t *fof_group,
  int          n_fof_groups,
  int          n_requested_forests,
  int         *n_halos_kept,
  int         *n_fof_groups_kept,
  int         *index_lookup)
{
  // I guess this should ideally be equal to the chunk size of the input hdf5 file...
  int             buffer_size              = 5000;
  int             group_buffer_size        = 0; // This will be calculated below
  int             n_read                   = 0;
  int             n_to_read                = 0;
  bool            keep_flag;

  FILE           *fin_catalogs             = NULL;
  int             flayout_switch           = -1;
  int             i_catalog_file           = 0;
  int             n_halos_in_catalog_file  = 0;
  int             i_halo_in_catalog_file   = 0;
  int             i_halo                   = 0;
  int             Len;

  FILE           *fin_groups               = NULL;
  int             group_flayout_switch     = -1;
  int             i_group_file             = 0;
  int             n_groups_in_catalog_file = 0;
  int             i_group_in_catalog_file  = 0;
  int             i_group                  = 0;
  int             n_groups                 = 0;
  int             first_group_index        = 0;
  int             last_group_index         = 1; // N.B. Must be init to different value from first_group_index


  char            catalog_file_prefix[50];
  char            simulation_dir[STRLEN];
  catalog_halo_t *catalog_buffer;
  catalog_halo_t *group_buffer;

  mlog("Doing read...", MLOG_OPEN);

  sprintf(catalog_file_prefix, "%s", run_globals.params.CatalogFilePrefix);
  sprintf(simulation_dir, "%s", run_globals.params.SimulationDir);

  tree_entry_t *tree_buffer;
  tree_buffer        = malloc(sizeof(tree_entry_t) * buffer_size);
  catalog_buffer     = malloc(sizeof(catalog_halo_t) * buffer_size);

  *n_halos_kept      = 0;
  *n_fof_groups_kept = 0;

  size_t dst_size       = sizeof(tree_entry_t);
  size_t dst_offsets[10] = {
    HOFFSET(tree_entry_t, id),
    HOFFSET(tree_entry_t, flags),
    HOFFSET(tree_entry_t, desc_id),
    HOFFSET(tree_entry_t, tree_id),
    HOFFSET(tree_entry_t, file_offset),
    HOFFSET(tree_entry_t, desc_index),
    HOFFSET(tree_entry_t, n_particle_peak),
    HOFFSET(tree_entry_t, central_index),
    HOFFSET(tree_entry_t, forest_id),
    HOFFSET(tree_entry_t, group_index)
  };
  size_t dst_sizes[10] = {
    sizeof(tree_buffer[0].id),
    sizeof(tree_buffer[0].flags),
    sizeof(tree_buffer[0].desc_id),
    sizeof(tree_buffer[0].tree_id),
    sizeof(tree_buffer[0].file_offset),
    sizeof(tree_buffer[0].desc_index),
    sizeof(tree_buffer[0].n_particle_peak),
    sizeof(tree_buffer[0].central_index),
    sizeof(tree_buffer[0].forest_id),
    sizeof(tree_buffer[0].group_index)
  };

  // Calculate the maximum group buffer size we require to hold a single
  // buffer worth of subgroups.  This is necessary as subfind can produce
  // many groups with no subgroups in them and so the group buffer may need
  // to be larger than the subgroup buffer.
  if (run_globals.mpi_rank == 0)
    while (n_read < n_halos)
    {
      if ((n_halos - n_read) >= buffer_size)
        n_to_read = buffer_size;
      else
        n_to_read = n_halos - n_read;

      H5TBread_records(fd, "trees", n_read, (hsize_t)n_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);
      n_read += n_to_read;

      int tmp_size = tree_buffer[n_to_read - 1].group_index - tree_buffer[0].group_index + 1;
      if (tmp_size > group_buffer_size)
        group_buffer_size = tmp_size;
    }
  MPI_Bcast(&group_buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  mlog("Using group buffer size = %d", MLOG_MESG, group_buffer_size);
  group_buffer = malloc(sizeof(catalog_halo_t) * group_buffer_size);

  // Now actually do the buffered read...
  n_read       = 0;
  n_to_read    = 0;
  keep_flag    = true;
  while (n_read < n_halos)
  {
    if ((n_halos - n_read) >= buffer_size)
      n_to_read = buffer_size;
    else
      n_to_read = n_halos - n_read;

    // read in a tree_buffer of the trees
    if (run_globals.mpi_rank == 0)
      H5TBread_records(fd, "trees", n_read, (hsize_t)n_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);
    MPI_Bcast(tree_buffer, n_to_read * sizeof(tree_entry_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    first_group_index = tree_buffer[0].group_index;
    if (first_group_index == last_group_index)
    {
      i_group = 1;
      memcpy(&(group_buffer[0]), &(group_buffer[n_groups - 1]), sizeof(catalog_halo_t));
    }
    else
      i_group = 0;
    last_group_index = tree_buffer[n_to_read - 1].group_index;
    n_groups         = last_group_index - first_group_index + 1;

    // read in the corresponding catalog entrys
    if (run_globals.mpi_rank == 0)
    {
      i_halo = 0;
      read_catalog_halos(&fin_catalogs, simulation_dir, catalog_file_prefix,
                         unsampled_snapshot, &flayout_switch, &i_catalog_file, &n_halos_in_catalog_file,
                         &i_halo_in_catalog_file, &i_halo, catalog_buffer, n_to_read, 1);
      read_catalog_halos(&fin_groups, simulation_dir, catalog_file_prefix,
                         unsampled_snapshot, &group_flayout_switch, &i_group_file, &n_groups_in_catalog_file,
                         &i_group_in_catalog_file, &i_group, group_buffer, n_groups - i_group, 0);
    }
    MPI_Bcast(catalog_buffer, n_to_read * sizeof(catalog_halo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(group_buffer, n_groups * sizeof(catalog_halo_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    // paste the data into the halo structures
    for (int jj = 0; jj < n_to_read; jj++)
    {
      if (run_globals.RequestedForestId != NULL)
      {
        if (bsearch(&(tree_buffer[jj].forest_id), run_globals.RequestedForestId,
                    (size_t)n_requested_forests, sizeof(int), compare_ints) != NULL)
          keep_flag = true;
        else
          keep_flag = false;
      }

      if (keep_flag)
      {
        halo_t         *cur_halo       = &(halo[*n_halos_kept]);
        catalog_halo_t *cur_cat_halo   = &(catalog_buffer[jj]);
        tree_entry_t   *cur_tree_entry = &(tree_buffer[jj]);

        cur_halo->ID                 = cur_tree_entry->id;
        cur_halo->TreeFlags          = cur_tree_entry->flags;
        cur_halo->SnapOffset         = cur_tree_entry->file_offset;
        cur_halo->DescIndex          = cur_tree_entry->desc_index;
        cur_halo->NextHaloInFOFGroup = NULL;
        cur_halo->ForestID           = cur_tree_entry->forest_id;

        if (index_lookup)
          index_lookup[*n_halos_kept] = n_read + jj;

        if (n_read + jj == tree_buffer[jj].central_index)
        {
          cur_halo->Type = 0;

          assert((*n_fof_groups_kept) < run_globals.NFOFGroupsMax);
          assert((tree_buffer[jj].group_index - first_group_index) < n_groups);
          catalog_halo_t *cur_cat_group = &(group_buffer[tree_buffer[jj].group_index - first_group_index]);
          fof_group_t    *cur_group     = &(fof_group[*n_fof_groups_kept]);

          cur_group->Mvir = cur_cat_group->M_vir;
          cur_group->Rvir = cur_cat_group->R_vir;
          cur_group->FOFMvirModifier = 1.0;

          convert_input_virial_props(&(cur_group->Mvir),
                                     &(cur_group->Rvir),
                                     &(cur_group->Vvir),
                                     &(cur_group->FOFMvirModifier),
                                     -1,
                                     snapshot);

          fof_group[(*n_fof_groups_kept)++].FirstHalo = &(halo[*n_halos_kept]);
        }
        else
        {
          cur_halo->Type                               = 1;
          halo[(*n_halos_kept) - 1].NextHaloInFOFGroup = &(halo[*n_halos_kept]);
        }

        cur_halo->FOFGroup  = &(fof_group[(*n_fof_groups_kept) - 1]);

        // paste in the halo properties
        cur_halo->id_MBP    = cur_cat_halo->id_MBP;
        cur_halo->Len       = cur_cat_halo->n_particles;
        cur_halo->Pos[0]    = cur_cat_halo->position_MBP[0];
        cur_halo->Pos[1]    = cur_cat_halo->position_MBP[1];
        cur_halo->Pos[2]    = cur_cat_halo->position_MBP[2];
        cur_halo->Vel[0]    = cur_cat_halo->velocity_COM[0];
        cur_halo->Vel[1]    = cur_cat_halo->velocity_COM[1];
        cur_halo->Vel[2]    = cur_cat_halo->velocity_COM[2];
        cur_halo->Rvir      = cur_cat_halo->R_vir;
        cur_halo->Rhalo     = cur_cat_halo->R_halo;
        cur_halo->Rmax      = cur_cat_halo->R_max;
        cur_halo->Vmax      = cur_cat_halo->V_max;
        cur_halo->VelDisp   = cur_cat_halo->sigma_v;
        cur_halo->AngMom[0] = cur_cat_halo->ang_mom[0];
        cur_halo->AngMom[1] = cur_cat_halo->ang_mom[1];
        cur_halo->AngMom[2] = cur_cat_halo->ang_mom[2];
        cur_halo->Galaxy    = NULL;
        cur_halo->Mvir      = cur_cat_halo->M_vir;

        // double check that PBC conditions are met!
        cur_halo->Pos[0]    = apply_pbc_pos(cur_halo->Pos[0]);
        cur_halo->Pos[1]    = apply_pbc_pos(cur_halo->Pos[1]);
        cur_halo->Pos[2]    = apply_pbc_pos(cur_halo->Pos[2]);

        // TODO: sort this out once and for all!
        if ((cur_halo->Type == 0) && run_globals.params.FlagSubhaloVirialProps)
          Len = -1;
        else
          Len = cur_halo->Len;

        convert_input_virial_props(&(cur_halo->Mvir),
                                   &(cur_halo->Rvir),
                                   &(cur_halo->Vvir),
                                   NULL,
                                   Len,
                                   snapshot);

        // // Replace the virial properties of the FOF group by those of the first
        // // subgroup
        // if (cur_halo->Type == 0)
        // {
        //   cur_halo->FOFGroup->Mvir = cur_halo->Mvir;
        //   cur_halo->FOFGroup->Rvir = cur_halo->Rvir;
        //   cur_halo->FOFGroup->Vvir = cur_halo->Vvir;
        // }

        (*n_halos_kept)++;
      }
    }

    n_read += n_to_read;
  }

  // close the catalogs file
  if (run_globals.mpi_rank == 0)
    if (fin_catalogs)
      fclose(fin_catalogs);

  // free the buffers
  free(group_buffer);
  free(catalog_buffer);
  free(tree_buffer);

  mlog(" ...done", MLOG_CLOSE);
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
  H5LTget_attribute_int(fd, "trees", "unsampled_snapshot", &(trees_info.unsampled_snapshot));

  return trees_info;
}


static void reorder_forest_array(int *arr, size_t *sort_ind, int n_forests, int *temp)
{
  memcpy(temp, arr, sizeof(int) * n_forests);
  for(int ii=0, jj=n_forests-1; ii < n_forests; ii++, jj--)
    arr[ii] = temp[sort_ind[jj]];
}


static void select_forests()
{
  char  fname[STRLEN];
  hid_t fin;
  int   i_forest;
  int   i_req;
  int   n_forests;
  int   n_halos_tot;
  int   n_halos_unused;
  int  *forest_id;
  int  *max_contemp_halo;
  int  *max_contemp_fof;
  int  *n_halos;
  int   max_halos;
  int   max_fof_groups;
  int  *n_requested_forests = &(run_globals.NRequestedForests);
  int  *requested_ind;
  bool  sample_forests;

  mlog("Calling select_forests()...", MLOG_MESG);

  // are we sampling the forests or just dividing them amoungst cores?
  if (run_globals.RequestedForestId == NULL)
    sample_forests = false;
  else
    sample_forests = true;

  // if this is the master rank then read in the forest info
  if (run_globals.mpi_rank == 0)
  {
    sprintf(fname, "%s/trees/forests_info.hdf5", run_globals.params.SimulationDir);
    if ((fin = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(EXIT_FAILURE);
    }

    // find out how many forests there are
    H5LTget_attribute_int(fin, "info", "n_forests", &n_forests);

    // find out how many halos in the whole simulation there are
    H5LTget_attribute_int(fin, "info", "n_subgroups", &n_halos_tot);
    H5LTget_attribute_int(fin, "info", "n_unused", &n_halos_unused);
    n_halos_tot     += n_halos_unused;

    // allocate the arrays
    forest_id        = (int*)malloc(sizeof(int) * n_forests);
    max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
    max_contemp_fof  = (int*)malloc(sizeof(int) * n_forests);
    n_halos          = (int*)malloc(sizeof(int) * n_forests);

    // read in the max number of contemporaneous halos and groups and the forest ids
    H5LTread_dataset_int(fin, "info/max_contemporaneous_halos", max_contemp_halo);
    H5LTread_dataset_int(fin, "info/max_contemporaneous_fof_groups", max_contemp_fof);
    H5LTread_dataset_int(fin, "info/forest_id", forest_id);
    H5LTread_dataset_int(fin, "info/n_halos", n_halos);

    // close the file
    H5Fclose(fin);
  }

  // broadcast the forest info
  MPI_Bcast(&n_forests, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_halos_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_halos_unused, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (run_globals.mpi_rank > 0)
  {
    forest_id        = (int*)malloc(sizeof(int) * n_forests);
    max_contemp_halo = (int*)malloc(sizeof(int) * n_forests);
    max_contemp_fof  = (int*)malloc(sizeof(int) * n_forests);
    n_halos          = (int*)malloc(sizeof(int) * n_forests);
  }
  MPI_Bcast(forest_id, n_forests, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(max_contemp_halo, n_forests, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(max_contemp_fof, n_forests, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_halos, n_forests, MPI_INT, 0, MPI_COMM_WORLD);

  // sort the forests by the number of halos in each one
  size_t *sort_ind = malloc(sizeof(size_t) * n_forests);
  gsl_sort_int_index(sort_ind, n_halos, 1, n_forests);
  int *temp = malloc(sizeof(int) * n_forests);

  reorder_forest_array(forest_id, sort_ind, n_forests, temp);
  reorder_forest_array(max_contemp_halo, sort_ind, n_forests, temp);
  reorder_forest_array(max_contemp_fof, sort_ind, n_forests, temp);
  reorder_forest_array(n_halos, sort_ind, n_forests, temp);

  free(temp);
  free(sort_ind);
  
  // if we have read in a list of requested forest IDs then use these to create
  // an array of indices pointing to the elements we want
  if (sample_forests)
  {
    requested_ind = (int*)malloc(sizeof(int) * (*n_requested_forests));
    n_halos_tot   = 0;
    for (i_forest = 0, i_req = 0; (i_forest < n_forests) && (i_req < (*n_requested_forests)); i_forest++)
      if (forest_id[i_forest] == run_globals.RequestedForestId[i_req])
      {
        requested_ind[i_req] = i_forest;
        n_halos_tot         += n_halos[i_forest];
        i_req++;
      }
    n_forests = (*n_requested_forests);
  }
  else
  {
    // if we haven't asked for any specific forest IDs then just fill the
    // requested ind array sequentially
    requested_ind = (int*)malloc(sizeof(int) * n_forests);
    for (i_req = 0; i_req < n_forests; i_req++)
      requested_ind[i_req] = i_req;
  }

  // if we are running the code with more than one core, let's subselect the forests on each one
  if (run_globals.mpi_size > 1)
  {
    if (n_forests < run_globals.mpi_size)
    {
      mlog_error("There are fewer processors than there are forests to be processed (%d).  Try again with fewer cores!", n_forests);
      ABORT(EXIT_FAILURE);
    }

    // TODO: Load balancing still seems pretty bad.
    // This should be looked into to see if there are any improvements that can
    // be made.  This would hopefully improve runtimes...

    // malloc the arrays we need for keeping track of the load balancing
    int *rank_first_forest = (int*)malloc(sizeof(int) * run_globals.mpi_size);
    int *rank_last_forest  = (int*)malloc(sizeof(int) * run_globals.mpi_size);
    int *rank_n_forests    = (int*)malloc(sizeof(int) * run_globals.mpi_size);
    int *rank_n_halos      = (int*)malloc(sizeof(int) * run_globals.mpi_size);
    int  i_rank;
    int  n_halos_used;
    int  n_halos_target;

    // initialise
    for (int ii = 0; ii < run_globals.mpi_size; ii++)
    {
      rank_first_forest[ii] = 0;
      rank_last_forest[ii]  = 0;
      rank_n_forests[ii]    = 0;
      rank_n_halos[ii]      = 0;
    }

    // loop through each rank
    for (i_rank = 0, n_halos_used = 0, i_forest = 0, n_halos_target = 0; i_rank < run_globals.mpi_size; i_rank++, i_forest++)
    {
      // start with the next forest (or first forest if this is the first rank)
      rank_first_forest[i_rank] = i_forest;
      rank_last_forest[i_rank]  = i_forest;

      // if we still have forests left to check from the total list of forests
      if (i_forest < n_forests)
      {
        // add this forest's halos to this rank's total
        rank_n_halos[i_rank] += n_halos[requested_ind[i_forest]];

        // adjust our target to smooth out variations as much as possible
        n_halos_target        = (n_halos_tot - n_halos_used) / (run_globals.mpi_size - i_rank);

        // increment the counter of the number of forests on this rank
        rank_n_forests[i_rank]++;

        // keep adding forests until we have reached (or exceeded) the current target
        while ((rank_n_halos[i_rank] < n_halos_target) && (i_forest < (n_forests - 1)))
        {
          i_forest++;
          rank_n_forests[i_rank]++;
          rank_last_forest[i_rank] = i_forest;
          rank_n_halos[i_rank]    += n_halos[requested_ind[i_forest]];
        }

        // updated the total number of used halos
        n_halos_used += rank_n_halos[i_rank];
      }
    }

    // add any uncounted forests to the last process
    for (i_forest++; i_forest < n_forests; i_forest++)
    {
      rank_last_forest[run_globals.mpi_size - 1] = i_forest;
      rank_n_halos[run_globals.mpi_size - 1]    += n_halos[requested_ind[i_forest]];
      rank_n_forests[run_globals.mpi_size - 1]++;
    }

    // create our list of forest_ids for this rank
    *n_requested_forests = rank_n_forests[run_globals.mpi_rank];
    if (run_globals.RequestedForestId != NULL)
      run_globals.RequestedForestId = reallocf(run_globals.RequestedForestId, *n_requested_forests * sizeof(int));
    else
      run_globals.RequestedForestId = (int*)malloc(sizeof(int) * (*n_requested_forests));

    for (int ii = rank_first_forest[run_globals.mpi_rank], jj = 0; ii <= rank_last_forest[run_globals.mpi_rank]; ii++, jj++)
      run_globals.RequestedForestId[jj] = forest_id[requested_ind[ii]];

    // free the arrays
    free(rank_n_halos);
    free(rank_n_forests);
    free(rank_last_forest);
    free(rank_first_forest);
  }


  // loop through and tot up the max number of halos and fof_groups we will need to allocate
  max_halos      = 0;
  max_fof_groups = 0;
  for (int i_forest = 0, i_req = 0; (i_forest < n_forests) && (i_req < (*n_requested_forests)); i_forest++)
    if (forest_id[requested_ind[i_forest]] == run_globals.RequestedForestId[i_req])
    {
      max_halos      += max_contemp_halo[requested_ind[i_forest]];
      max_fof_groups += max_contemp_fof[requested_ind[i_forest]];
      i_req++;
    }

  // store the maximum number of halos and fof groups needed at any one snapshot
  run_globals.NHalosMax     = max_halos;
  run_globals.NFOFGroupsMax = max_fof_groups;

  // sort the requested forest ids so that they can be bsearch'd later
  qsort(run_globals.RequestedForestId, *n_requested_forests, sizeof(int), compare_ints);

  // {
  //   // DEBUG
  //   hid_t fd;
  //   char fname[50];
  //   hsize_t dims;
  //   sprintf(fname, "debug_%d.hdf5", run_globals.mpi_rank);
  //   fd = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //   dims = (hsize_t)(*n_requested_forests);
  //   H5LTmake_dataset_int(fd, "forest_id", 1, &dims, run_globals.RequestedForestId);
  //   dims = (hsize_t)n_forests;
  //   H5LTmake_dataset_int(fd, "requested_ind", 1, &dims, requested_ind);
  //   H5LTmake_dataset_int(fd, "n_halos", 1, &dims, n_halos);
  //   H5Fclose(fd);
  // }

  // free the arrays
  free(requested_ind);
  free(n_halos);
  free(forest_id);
  free(max_contemp_fof);
  free(max_contemp_halo);
}


static void inline gen_dump_fname(char *fname, int snapshot)
{
  sprintf(fname, "%s/dump_files/%s-rank%d.%d",
          run_globals.params.OutputDir, run_globals.params.FileNameGalaxies,
          run_globals.mpi_rank, snapshot);
}


static void inline update_pointers_from_offsets(
  int          n_fof_groups_kept,
  fof_group_t *fof_group,
  size_t      *fof_FirstHalo_os,
  int          n_halos_kept,
  halo_t      *halo,
  size_t      *halo_FOFGroup_os,
  size_t      *halo_NextHaloInFOFGroup_os)
{
  for (int ii = 0; ii < n_halos_kept; ii++)
  {
    halo[ii].FOFGroup = &(fof_group[halo_FOFGroup_os[ii]]);
    if (halo_NextHaloInFOFGroup_os[ii] != -1)
      halo[ii].NextHaloInFOFGroup = &(halo[halo_NextHaloInFOFGroup_os[ii]]);
    else
      halo[ii].NextHaloInFOFGroup = NULL;
  }
  for (int ii = 0; ii < n_fof_groups_kept; ii++)
    fof_group[ii].FirstHalo = &(halo[fof_FirstHalo_os[ii]]);
}


trees_info_t read_halos(
  int           snapshot,
  halo_t      **halo,
  fof_group_t **fof_group,
  int         **index_lookup,
  trees_info_t *snapshot_trees_info)
{
  int          n_halos;                    //!< Number of halos
  char         fname[STRLEN];
  int          n_fof_groups;
  trees_info_t trees_info;
  hid_t        fin_trees;

  int          n_halos_kept;
  int          n_fof_groups_kept;
  int          n_requested_forests = run_globals.NRequestedForests;

  // if we are doing multiple runs and have already read in this snapshot then we don't need to do read anything else
  if ((run_globals.params.FlagInteractive || run_globals.params.FlagMCMC) && (snapshot_trees_info[snapshot].n_halos != -1))
  {
    mlog("Snapshot %d has alread been read in... (n_halos = %d)", MLOG_MESG, snapshot, snapshot_trees_info[snapshot].n_halos);
    return snapshot_trees_info[snapshot];
  }

  mlog("Reading snapshot %d (z = %.2f) trees and halos...", MLOG_OPEN | MLOG_TIMERSTART, snapshot, run_globals.ZZ[snapshot]);

  // Read mass ratio modifiers and baryon fraction modifiers if required
  if (run_globals.RequestedMassRatioModifier == 1)
    read_mass_ratio_modifiers(snapshot);

  if (run_globals.RequestedBaryonFracModifier == 1)
    read_baryon_frac_modifiers(snapshot);

  // open the tree file
  if (run_globals.mpi_rank == 0)
  {
    sprintf(fname, "%s/trees/horizontal_trees_%03d.hdf5", run_globals.params.SimulationDir, snapshot);
    if ((fin_trees = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(EXIT_FAILURE);
    }

    // read the info attributes
    trees_info = read_trees_info(fin_trees);
  }

  // if necessary, broadcast the tree file info
  MPI_Bcast(&trees_info, sizeof(trees_info_t), MPI_BYTE, 0, MPI_COMM_WORLD);

  n_halos      = trees_info.n_halos;
  n_fof_groups = trees_info.n_fof_groups;

  // if we haven't already, store the TreesStep & TreesScan parameters
  if (run_globals.TreesStep == -1)
    run_globals.TreesStep = trees_info.n_step;
  if (run_globals.TreesScan == -1)
    run_globals.TreesScan = trees_info.n_search;

  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC)
    if (n_halos < 1)
    {
      mlog("No halos in this file... skipping...", MLOG_CLOSE | MLOG_TIMERSTOP);
      if (run_globals.mpi_rank == 0)
        H5Fclose(fin_trees);
      snapshot_trees_info[snapshot] = trees_info;
      return trees_info;
    }

  // If necessary, allocate the halo array
  if (*halo == NULL)
  {
    // if required, select forests and calculate the maximum number of halos and fof groups
    if ((n_requested_forests > -1) || (run_globals.mpi_size > 1))
    {
      if (run_globals.SelectForestsSwitch == true)
      {
        select_forests();
        n_requested_forests             = run_globals.NRequestedForests;

        // Toggle the SelectForestsSwitch so that we know we don't need to call this again
        run_globals.SelectForestsSwitch = false;
      }
      *index_lookup = malloc(sizeof(int) * run_globals.NHalosMax);
    }
    else
    {
      run_globals.NHalosMax     = trees_info.n_halos_max;
      run_globals.NFOFGroupsMax = trees_info.n_fof_groups_max;
    }

    mlog("Allocating halo array with %d elements...", MLOG_MESG, run_globals.NHalosMax);
    *halo = malloc(sizeof(halo_t) * run_globals.NHalosMax);
  }

  // Allocate the fof_group array if necessary
  if (*fof_group == NULL)
  {
    mlog("Allocating fof_group array with %d elements...", MLOG_MESG, run_globals.NFOFGroupsMax);
    *fof_group = malloc(sizeof(fof_group_t) * run_globals.NFOFGroupsMax);
  }

  // reset the fof group pointers and index lookup (if necessary)
  for (int ii = 0; ii < run_globals.NFOFGroupsMax; ii++)
  {
    (*fof_group)[ii].FirstHalo         = NULL;
    (*fof_group)[ii].FirstOccupiedHalo = NULL;
    (*fof_group)[ii].Mvir              = 0.0;
  }
  if (n_requested_forests > -1)
    for (int ii = 0; ii < run_globals.NHalosMax; ii++)
      (*index_lookup)[ii] = -1;

  if (n_halos < 1)
  {
    mlog("No halos in this file... skipping...", MLOG_CLOSE | MLOG_TIMERSTOP);
    if (run_globals.mpi_rank == 0)
      H5Fclose(fin_trees);
    return trees_info;
  }


  if (run_globals.params.FlagReadDumpFile)
  {
    FILE *fdump = NULL;
    char  fname[STRLEN];

    gen_dump_fname(fname, snapshot);

    if ((fdump = fopen(fname, "rb")) != NULL)
    {
      fread(&n_fof_groups_kept, sizeof(int), 1, fdump);
      size_t *fof_FirstHalo_os = malloc(sizeof(size_t) * n_fof_groups_kept);
      fread(*fof_group, sizeof(fof_group_t), n_fof_groups_kept, fdump);
      fread(fof_FirstHalo_os, sizeof(size_t), n_fof_groups_kept, fdump);

      fread(&n_halos_kept, sizeof(int), 1, fdump);
      size_t *halo_FOFGroup_os           = malloc(sizeof(size_t) * n_halos_kept);
      size_t *halo_NextHaloInFOFGroup_os = malloc(sizeof(size_t) * n_halos_kept);
      fread(*halo, sizeof(halo_t), n_halos_kept, fdump);
      fread(halo_FOFGroup_os, sizeof(size_t), n_halos_kept, fdump);
      fread(halo_NextHaloInFOFGroup_os, sizeof(size_t), n_halos_kept, fdump);
      fread(*index_lookup, sizeof(int), n_halos_kept, fdump);

      update_pointers_from_offsets(n_fof_groups_kept, *fof_group, fof_FirstHalo_os,
                                   n_halos_kept, *halo, halo_FOFGroup_os, halo_NextHaloInFOFGroup_os);

      free(fof_FirstHalo_os);
      free(halo_NextHaloInFOFGroup_os);
      free(halo_FOFGroup_os);

      fclose(fdump);
    }
    else
    {
      mlog_error("FlagReadDumpFile = 1, but there is no dumpfile at:\n%s", fname);
      ABORT(EXIT_FAILURE);
    }
  }
  else
    // read in the trees
    read_trees_and_catalogs(snapshot, trees_info.unsampled_snapshot,
                            fin_trees, *halo, n_halos, *fof_group, n_fof_groups, n_requested_forests,
                            &n_halos_kept, &n_fof_groups_kept, *index_lookup);

  // close the tree file
  if (run_globals.mpi_rank == 0)
    H5Fclose(fin_trees);

  // if subsampling the trees, then update the trees_info to reflect what we now have
  if (n_requested_forests > -1)
  {
    trees_info.n_halos      = n_halos_kept;
    trees_info.n_fof_groups = n_fof_groups_kept;
  }

  // if we are doing multiple runs then resize the arrays to save space and store the trees_info
  if (run_globals.params.FlagInteractive || run_globals.params.FlagGenDumpFile || run_globals.params.FlagMCMC)
  {
    // Ok - what follows here is hacky as hell.  By calling realloc on these
    // arrays, there is a good chance that the actual array will be moved and
    // there is no way to prevent this.  A side effect will be that all of the
    // pointers that refer to any of these arrays will be broken.  The simplest
    // way to deal with this is to calculate the array offsets which each
    // pointer refers to, and then reset all of the pointers in the realloc'd
    // arrays.
    //
    // In future, this should definitely be changed. All pointers to these
    // array members should be made integer offsets instead.  That will require
    // trawling through the code and making the relevant updates in a number of
    // places (so we'll leave that till later!).

    mlog("Calculating pointer offsets...", MLOG_MESG);

    size_t *halo_FOFGroup_os           = malloc(sizeof(size_t) * n_halos_kept);
    size_t *halo_NextHaloInFOFGroup_os = malloc(sizeof(size_t) * n_halos_kept);
    size_t *fof_FirstHalo_os           = malloc(sizeof(size_t) * n_fof_groups_kept);

    for (int ii = 0; ii < n_halos_kept; ii++)
    {
      halo_FOFGroup_os[ii] = (size_t)((*halo)[ii].FOFGroup - (*fof_group));
      if ((*halo)[ii].NextHaloInFOFGroup != NULL)
        halo_NextHaloInFOFGroup_os[ii] = (size_t)((*halo)[ii].NextHaloInFOFGroup - (*halo));
      else
        halo_NextHaloInFOFGroup_os[ii] = -1;
    }
    for (int ii = 0; ii < n_fof_groups_kept; ii++)
      fof_FirstHalo_os[ii] = (size_t)((*fof_group)[ii].FirstHalo - (*halo));

    if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC)
    {
      mlog("Reallocing halo storage arrays...", MLOG_OPEN);

      *halo      = (halo_t*)realloc(*halo, sizeof(halo_t) * n_halos_kept);
      *fof_group = (fof_group_t*)realloc(*fof_group, sizeof(fof_group_t) * n_fof_groups_kept);

      // Only realloc `index_lookup` if we are actually using it (`n_proc` > 1
      // or we are subsampling the trees).
      if (*index_lookup)
        *index_lookup = (int*)realloc(*index_lookup, sizeof(int) * n_halos_kept);
      update_pointers_from_offsets(n_fof_groups_kept, *fof_group, fof_FirstHalo_os,
                                   n_halos_kept, *halo, halo_FOFGroup_os, halo_NextHaloInFOFGroup_os);

      // save the trees_info for this snapshot as well...
      snapshot_trees_info[snapshot] = trees_info;

      mlog(" ...done (resized to %d halos)", MLOG_CLOSE, n_halos_kept);
    }
    else
    {
      mlog("Generating dumpfile...", MLOG_OPEN);

      FILE *fdump = NULL;
      char  fname[STRLEN];

      gen_dump_fname(fname, snapshot);

      if ((fdump = fopen(fname, "wb")) != NULL)
      {
        fwrite(&n_fof_groups_kept, sizeof(int), 1, fdump);
        fwrite(*fof_group, sizeof(fof_group_t), n_fof_groups_kept, fdump);
        fwrite(fof_FirstHalo_os, sizeof(size_t), n_fof_groups_kept, fdump);

        fwrite(&n_halos_kept, sizeof(int), 1, fdump);
        fwrite(*halo, sizeof(halo_t), n_halos_kept, fdump);
        fwrite(halo_FOFGroup_os, sizeof(size_t), n_halos_kept, fdump);
        fwrite(halo_NextHaloInFOFGroup_os, sizeof(size_t), n_halos_kept, fdump);
        fwrite(*index_lookup, sizeof(int), n_halos_kept, fdump);

        fclose(fdump);
      }
      mlog("...done", MLOG_CLOSE);
    }

    free(fof_FirstHalo_os);
    free(halo_NextHaloInFOFGroup_os);
    free(halo_FOFGroup_os);
  }

  MPI_Allreduce(MPI_IN_PLACE, &n_halos_kept, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &n_fof_groups_kept, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mlog("Read %d halos in %d fof_groups.", MLOG_MESG, n_halos_kept, n_fof_groups_kept);

  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  return trees_info;
}


void initialize_halo_storage()
{
  int           *n_store_snapshots     = &(run_globals.NStoreSnapshots);
  halo_t      ***snapshot_halo         = &(run_globals.SnapshotHalo);
  fof_group_t ***snapshot_fof_group    = &(run_globals.SnapshotFOFGroup);
  int         ***snapshot_index_lookup = &(run_globals.SnapshotIndexLookup);
  trees_info_t **snapshot_trees_info   = &(run_globals.SnapshotTreesInfo);

  int            last_snap             = 0;

  mlog("Initializing halo storage arrays...", MLOG_OPEN);

  // Find what the last requested output snapshot is
  for (int ii = 0; ii < run_globals.NOutputSnaps; ii++)
    if (run_globals.ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals.ListOutputSnaps[ii];

  // Allocate an array of last_snap halo array pointers
  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC)
    *n_store_snapshots = last_snap + 1;
  else
    *n_store_snapshots = 1;

  *snapshot_halo         = (halo_t**)calloc(*n_store_snapshots, sizeof(halo_t *));
  *snapshot_fof_group    = (fof_group_t**)calloc(*n_store_snapshots, sizeof(fof_group_t *));
  *snapshot_index_lookup = (int**)calloc(*n_store_snapshots, sizeof(int *));
  *snapshot_trees_info   = (trees_info_t*)calloc(*n_store_snapshots, sizeof(trees_info_t));

  for (int ii = 0; ii < *n_store_snapshots; ii++)
  {
    (*snapshot_trees_info)[ii].n_halos = -1;
    (*snapshot_index_lookup)[ii]       = NULL;
  }

  // loop through and read all snapshots
  if (run_globals.params.FlagInteractive || run_globals.params.FlagMCMC)
  {
    mlog("Preloading input trees and halos...", MLOG_OPEN);
    for (int i_snap = 0; i_snap <= last_snap; i_snap++)
      read_halos(i_snap, &((*snapshot_halo)[i_snap]), &((*snapshot_fof_group)[i_snap]), &((*snapshot_index_lookup)[i_snap]), *snapshot_trees_info);
    mlog("...done", MLOG_CLOSE);
  }

  mlog("...done", MLOG_CLOSE);
}


void free_halo_storage()
{
  int           n_store_snapshots     = run_globals.NStoreSnapshots;
  halo_t      **snapshot_halo         = run_globals.SnapshotHalo;
  fof_group_t **snapshot_fof_group    = run_globals.SnapshotFOFGroup;
  int         **snapshot_index_lookup = run_globals.SnapshotIndexLookup;
  trees_info_t *snapshot_trees_info   = run_globals.SnapshotTreesInfo;

  // Free all of the remaining allocated galaxies, halos and fof groups
  mlog("Freeing FOF groups and halos...", MLOG_MESG);
  for (int ii = 0; ii < n_store_snapshots; ii++)
  {
    free(snapshot_halo[ii]);
    free(snapshot_fof_group[ii]);
    free(snapshot_index_lookup[ii]);
  }
  free(snapshot_halo);
  free(snapshot_fof_group);
  free(snapshot_index_lookup);
  free(snapshot_trees_info);
}
