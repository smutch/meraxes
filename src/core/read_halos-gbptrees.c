#include "meraxes.h"
#include <assert.h>
#include <hdf5_hl.h>
#include <math.h>
#include "tree_flags.h"

static void halo_catalog_filename(
    char* simulation_dir,
    char* catalog_file_prefix,
    int snapshot,
    char* group_type,
    int sub,
    int* i_layout,
    char* fname)
{

    // if we need to determine the filename structure...
    if (*i_layout == -1) {
        for (*i_layout = 0; *i_layout < 4; (*i_layout)++) {
            if (*i_layout == 0)
                sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
            else if (*i_layout == 1)
                sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
            else if (*i_layout == 2)
                sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
            else if (*i_layout == 3)
                sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);

            FILE* fin;
            if ((fin = fopen(fname, "rb")) != NULL) {
                fclose(fin);
                break;
            }
        }
    }

    // ensure we have a valid i_layout value.
    if (*i_layout < 0 || *i_layout > 3) {
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
    FILE* fin,
    int* i_file,
    int* n_files,
    int* n_halos_file,
    int* n_halos_total)
{
    fread(i_file, sizeof(int), 1, fin);
    fread(n_files, sizeof(int), 1, fin);
    fread(n_halos_file, sizeof(int), 1, fin);
    fread(n_halos_total, sizeof(int), 1, fin);
}

//! This is the structure for a halo in the catalog files
typedef struct catalog_halo_t {
    long long id_MBP; //!< ID of most bound particle in structure
    double M_vir; //!< Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
    int n_particles; //!< Number of particles in the structure
    float position_COM[3]; //!< Centre-of-mass position      [Mpc/h]
    float position_MBP[3]; //!< Most bound particle position [Mpc/h]
    float velocity_COM[3]; //!< Centre-of-mass velocity      [km/s]
    float velocity_MBP[3]; //!< Most bound particle velocity [km/s]
    float R_vir; //!< Virial radius [Mpc/h]
    float R_halo; //!< Distance of last halo particle from MBP [Mpc/h]
    float R_max; //!< Radius of maximum circular velocity     [Mpc/h]
    float V_max; //!< Maximum circular velocity               [km/s]
    float sigma_v; //!< Total 3D velocity dispersion            [km/s]
    float ang_mom[3]; //!< Specific angular momentum vector        [Mpc/h*km/s]
    float q_triaxial; //!< Triaxial shape parameter q=b/a
    float s_triaxial; //!< Triaxial shape parameter s=c/a
    float shape_eigen_vectors[3][3]; //!< Normalized triaxial shape eigenvectors
    char padding[8]; //!< Alignment padding
} catalog_halo_t;

static void read_catalog_halos(
    FILE** fin,
    char* simulation_dir,
    char* catalog_file_prefix,
    int snapshot,
    int* flayout_switch,
    int* i_file,
    int* n_halos_file,
    int* i_halo_in_file,
    int* i_halo,
    catalog_halo_t* halo,
    int n_to_read,
    int type_flag)
{
    char fname[STRLEN];
    char halo_type[10];
    int dummy;
    int n_from_this_file;

    switch (type_flag) {
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
    if ((*fin) == NULL) {
        halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, halo_type, *i_file, flayout_switch, fname);
        *fin = fopen(fname, "rb");
        if (*fin == NULL) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(34494);
        }
        read_catalogs_header(*fin, &dummy, &dummy, n_halos_file, &dummy);
    }

    // Have we already read all the halos in this file?
    if ((*i_halo_in_file) >= (*n_halos_file)) {
        fclose(*fin);
        (*i_file)++;
        (*i_halo_in_file) = 0;
        halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, halo_type, *i_file, flayout_switch, fname);
        *fin = fopen(fname, "rb");
        if (*fin == NULL) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(34494);
        }
        read_catalogs_header(*fin, &dummy, &dummy, n_halos_file, &dummy);
    }

    // Read in as many halos as we can from this file
    if ((*i_halo_in_file + n_to_read) <= *n_halos_file) {
        fread(&(halo[*i_halo]), sizeof(catalog_halo_t), n_to_read, *fin);
        *i_halo += n_to_read;
        *i_halo_in_file += n_to_read;
    }
    else {
        // read in as many as we can from this file and then get the rest from the next file
        n_from_this_file = (*n_halos_file) - *i_halo_in_file;

        fread(&(halo[*i_halo]), sizeof(catalog_halo_t), n_from_this_file, *fin);
        *i_halo += n_from_this_file;
        *i_halo_in_file += n_from_this_file;
        n_to_read -= n_from_this_file;
        read_catalog_halos(fin, simulation_dir, catalog_file_prefix, snapshot, flayout_switch, i_file, n_halos_file, i_halo_in_file, i_halo, halo, n_to_read, type_flag);
    }
}

static void inline convert_input_virial_props(double* Mvir, double* Rvir, double* Vvir, double* FOFMvirModifier, int len, int snapshot, bool flag_highres)
{
    if (len >= 0) {
        // Update the virial properties for subhalos
        *Mvir = calculate_Mvir(*Mvir, len, flag_highres);
        *Rvir = calculate_Rvir(*Mvir, snapshot);
    }
    else {
        // Convert the mass unit for FoFs
        *Mvir /= 1.0e10;
        if (run_globals.RequestedMassRatioModifier == 1) {
            // Modifier the FoF mass and update the virial radius
            *FOFMvirModifier = interpolate_modifier(run_globals.mass_ratio_modifier, log10(*Mvir / run_globals.params.Hubble_h) + 10.0);
            *Mvir *= *FOFMvirModifier;
            *Rvir = calculate_Rvir(*Mvir, snapshot);
        }
    }
    *Vvir = calculate_Vvir(*Mvir, *Rvir);
}

//! Buffered read of hdf5 trees into halo structures
void read_trees__gbptrees(
    int snapshot,
    halo_t* halo,
    int n_halos,
    fof_group_t* fof_group,
    int n_fof_groups,
    int n_requested_forests,
    int* n_halos_kept,
    int* n_fof_groups_kept,
    int* index_lookup)
{
    // I guess this should ideally be equal to the chunk size of the input hdf5 file...
    int buffer_size = 5000;
    int group_buffer_size = 0; // This will be calculated below
    int n_read = 0;
    int n_to_read = 0;
    bool keep_flag;

    FILE* fin_catalogs = NULL;
    int flayout_switch = -1;
    int i_catalog_file = 0;
    int n_halos_in_catalog_file = 0;
    int i_halo_in_catalog_file = 0;
    int i_halo = 0;
    int Len;

    FILE* fin_groups = NULL;
    int group_flayout_switch = -1;
    int i_group_file = 0;
    int n_groups_in_catalog_file = 0;
    int i_group_in_catalog_file = 0;
    int i_group = 0;
    int n_groups = 0;
    int first_group_index = 0;
    int last_group_index = 1; // N.B. Must be init to different value from first_group_index

    char catalog_file_prefix[50];
    char simulation_dir[STRLEN];
    char fname[STRLEN];
    hid_t fd;
    catalog_halo_t* catalog_buffer;
    catalog_halo_t* group_buffer;

    //! Tree entry struct
    typedef struct tree_entry_t {
        int id;
        int flags;
        int desc_id;
        int tree_id;
        int file_offset;
        int desc_index;
        int n_particle_peak;
        int central_index;
        int forest_id;
        int group_index;
    } tree_entry_t;

    mlog("Doing read...", MLOG_OPEN);

    if (run_globals.mpi_rank == 0) {
        // open the tree file
        sprintf(fname, "%s/trees/horizontal_trees_%03d.hdf5", run_globals.params.SimulationDir, snapshot);
        if ((fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(EXIT_FAILURE);
        }
    }

    sprintf(catalog_file_prefix, "%s", run_globals.params.CatalogFilePrefix);
    sprintf(simulation_dir, "%s", run_globals.params.SimulationDir);

    tree_entry_t* tree_buffer;
    tree_buffer = malloc(sizeof(tree_entry_t) * buffer_size);
    catalog_buffer = malloc(sizeof(catalog_halo_t) * buffer_size);

    *n_halos_kept = 0;
    *n_fof_groups_kept = 0;

    size_t dst_size = sizeof(tree_entry_t);
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
        while (n_read < n_halos) {
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
    MPI_Bcast(&group_buffer_size, 1, MPI_INT, 0, run_globals.mpi_comm);
    mlog("Using group buffer size = %d", MLOG_MESG, group_buffer_size);
    group_buffer = malloc(sizeof(catalog_halo_t) * group_buffer_size);

    // Now actually do the buffered read...
    n_read = 0;
    n_to_read = 0;
    keep_flag = true;
    while (n_read < n_halos) {
        if ((n_halos - n_read) >= buffer_size)
            n_to_read = buffer_size;
        else
            n_to_read = n_halos - n_read;

        // read in a tree_buffer of the trees
        if (run_globals.mpi_rank == 0)
            H5TBread_records(fd, "trees", n_read, (hsize_t)n_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);
        MPI_Bcast(tree_buffer, n_to_read * sizeof(tree_entry_t), MPI_BYTE, 0, run_globals.mpi_comm);

        first_group_index = tree_buffer[0].group_index;
        if (first_group_index == last_group_index) {
            i_group = 1;
            memcpy(&(group_buffer[0]), &(group_buffer[n_groups - 1]), sizeof(catalog_halo_t));
        }
        else
            i_group = 0;
        last_group_index = tree_buffer[n_to_read - 1].group_index;
        n_groups = last_group_index - first_group_index + 1;

        // read in the corresponding catalog entrys
        if (run_globals.mpi_rank == 0) {
            i_halo = 0;
            read_catalog_halos(&fin_catalogs, simulation_dir, catalog_file_prefix, snapshot,
                &flayout_switch, &i_catalog_file, &n_halos_in_catalog_file,
                &i_halo_in_catalog_file, &i_halo, catalog_buffer, n_to_read, 1);
            read_catalog_halos(&fin_groups, simulation_dir, catalog_file_prefix, snapshot,
                &group_flayout_switch, &i_group_file, &n_groups_in_catalog_file,
                &i_group_in_catalog_file, &i_group, group_buffer, n_groups - i_group, 0);
        }
        MPI_Bcast(catalog_buffer, n_to_read * sizeof(catalog_halo_t), MPI_BYTE, 0, run_globals.mpi_comm);
        MPI_Bcast(group_buffer, n_groups * sizeof(catalog_halo_t), MPI_BYTE, 0, run_globals.mpi_comm);

        // paste the data into the halo structures
        for (int jj = 0; jj < n_to_read; jj++) {
            if (run_globals.RequestedForestId != NULL) {
                if (bsearch(&(tree_buffer[jj].forest_id), run_globals.RequestedForestId,
                        (size_t)n_requested_forests, sizeof(int), compare_ints)
                    != NULL)
                    keep_flag = true;
                else
                    keep_flag = false;
            }

            if (keep_flag) {
                halo_t* cur_halo = &(halo[*n_halos_kept]);
                catalog_halo_t* cur_cat_halo = &(catalog_buffer[jj]);
                tree_entry_t* cur_tree_entry = &(tree_buffer[jj]);

                cur_halo->ID = cur_tree_entry->id;
                cur_halo->TreeFlags = cur_tree_entry->flags;
                cur_halo->SnapOffset = cur_tree_entry->file_offset;
                cur_halo->DescIndex = cur_tree_entry->desc_index;
                cur_halo->NextHaloInFOFGroup = NULL;

                if (index_lookup)
                    index_lookup[*n_halos_kept] = n_read + jj;

                if (n_read + jj == tree_buffer[jj].central_index) {
                    cur_halo->Type = 0;

                    assert((*n_fof_groups_kept) < run_globals.NFOFGroupsMax);
                    assert((tree_buffer[jj].group_index - first_group_index) < n_groups);
                    catalog_halo_t* cur_cat_group = &(group_buffer[tree_buffer[jj].group_index - first_group_index]);
                    fof_group_t* cur_group = &(fof_group[*n_fof_groups_kept]);

                    cur_group->Mvir = cur_cat_group->M_vir;
                    cur_group->Rvir = cur_cat_group->R_vir;
                    cur_group->FOFMvirModifier = 1.0;

                    convert_input_virial_props(&(cur_group->Mvir),
                        &(cur_group->Rvir),
                        &(cur_group->Vvir),
                        &(cur_group->FOFMvirModifier),
                        -1,
                        snapshot,
                        false);

                    fof_group[(*n_fof_groups_kept)++].FirstHalo = &(halo[*n_halos_kept]);
                }
                else {
                    cur_halo->Type = 1;
                    halo[(*n_halos_kept) - 1].NextHaloInFOFGroup = &(halo[*n_halos_kept]);
                }

                cur_halo->FOFGroup = &(fof_group[(*n_fof_groups_kept) - 1]);

                // paste in the halo properties
                cur_halo->Len = cur_cat_halo->n_particles;
                cur_halo->Pos[0] = cur_cat_halo->position_MBP[0];
                cur_halo->Pos[1] = cur_cat_halo->position_MBP[1];
                cur_halo->Pos[2] = cur_cat_halo->position_MBP[2];
                cur_halo->Vel[0] = cur_cat_halo->velocity_COM[0];
                cur_halo->Vel[1] = cur_cat_halo->velocity_COM[1];
                cur_halo->Vel[2] = cur_cat_halo->velocity_COM[2];
                cur_halo->Rvir = cur_cat_halo->R_vir;
                cur_halo->Vmax = cur_cat_halo->V_max;
                cur_halo->AngMom[0] = cur_cat_halo->ang_mom[0];
                cur_halo->AngMom[1] = cur_cat_halo->ang_mom[1];
                cur_halo->AngMom[2] = cur_cat_halo->ang_mom[2];
                cur_halo->Galaxy = NULL;
                cur_halo->Mvir = cur_cat_halo->M_vir;

                // double check that PBC conditions are met!
                cur_halo->Pos[0] = apply_pbc_pos(cur_halo->Pos[0]);
                cur_halo->Pos[1] = apply_pbc_pos(cur_halo->Pos[1]);
                cur_halo->Pos[2] = apply_pbc_pos(cur_halo->Pos[2]);

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
                    snapshot,
                    false);

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

    // close the catalogs and trees files
    if (run_globals.mpi_rank == 0) {
        H5Fclose(fd);

        if (fin_catalogs)
            fclose(fin_catalogs);
    }

    // free the buffers
    free(group_buffer);
    free(catalog_buffer);
    free(tree_buffer);

    mlog(" ...done", MLOG_CLOSE);
}

trees_info_t read_trees_info__gbptrees(int snapshot)
{
    trees_info_t trees_info;

    if (run_globals.mpi_rank == 0) {
        // open the tree file
        char fname[STRLEN];
        hid_t fd;

        sprintf(fname, "%s/trees/horizontal_trees_%03d.hdf5", run_globals.params.SimulationDir, snapshot);
        if ((fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(EXIT_FAILURE);
        }

        // read the info attributes
        H5LTget_attribute_int(fd, "trees", "n_halos", &(trees_info.n_halos));
        H5LTget_attribute_int(fd, "trees", "n_halos_max", &(trees_info.n_halos_max));
        H5LTget_attribute_int(fd, "trees", "max_tree_id", &(trees_info.max_tree_id));
        H5LTget_attribute_int(fd, "trees", "n_fof_groups", &(trees_info.n_fof_groups));
        H5LTget_attribute_int(fd, "trees", "n_fof_groups_max", &(trees_info.n_fof_groups_max));

        H5Fclose(fd);
    }

    // broadcast the tree file info
    MPI_Bcast(&trees_info, sizeof(trees_info_t), MPI_BYTE, 0, run_globals.mpi_comm);

    return trees_info;
}

//! Buffered read of hdf5 trees into halo structures (read the extended trees)
void read_trees__extended_gbptrees(
    int snapshot,
    halo_t* halo,
    int n_halos,
    fof_group_t* fof_group,
    int n_fof_groups,
    int n_requested_forests,
    int* n_halos_kept,
    int* n_fof_groups_kept,
    int* index_lookup)
{
    // I guess this should ideally be equal to the chunk size of the input hdf5 file...
    int buffer_size = 5000;
    int group_buffer_size = 0; // This will be calculated below
    int halo_buffer_size = 0; // This will be calculated below
    int n_read = 0;
    int n_to_read = 0;
    bool keep_flag;

    FILE* fin_catalogs = NULL;
    int flayout_switch = -1;
    int i_catalog_file = 0;
    int n_halos_in_catalog_file = 0;
    int i_halo_in_catalog_file = 0;
    int i_halo = 0;
    int n_catalogs = 0;
    int first_halo_index = 0;
    int last_halo_index = 1;
    int Len;

    FILE* fin_groups = NULL;
    int group_flayout_switch = -1;
    int i_group_file = 0;
    int n_groups_in_catalog_file = 0;
    int i_group_in_catalog_file = 0;
    int i_group = 0;
    int n_groups = 0;
    int first_group_index = 0;
    int last_group_index = 1; // N.B. Must be init to different value from first_group_index

    char catalog_file_prefix[50];
    char simulation_dir[STRLEN];
    char fname[STRLEN];
    hid_t fd;
    catalog_halo_t* catalog_buffer;
    catalog_halo_t* group_buffer;

    //! Tree entry struct
    typedef struct tree_entry_t {
        int id;
        int flags;
        int desc_id;
        int tree_id;
        int file_offset;
        int desc_index;
        int n_particle_peak;
        int central_index;
        int forest_id;
        int group_index;
        int original_index;
        float position_offset[3];
    } tree_entry_t;

    mlog("Doing read...", MLOG_OPEN);

    if (run_globals.mpi_rank == 0) {
        // open the tree file
        sprintf(fname, "%s/trees/horizontal_trees_%03d.hdf5", run_globals.params.SimulationDir, snapshot);
        if ((fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
            mlog("Failed to open file %s", MLOG_MESG, fname);
            ABORT(EXIT_FAILURE);
        }
    }

    sprintf(catalog_file_prefix, "%s", run_globals.params.CatalogFilePrefix);
    sprintf(simulation_dir, "%s", run_globals.params.SimulationDir);

    tree_entry_t* tree_buffer;
    tree_buffer = malloc(sizeof(tree_entry_t) * buffer_size);

    *n_halos_kept = 0;
    *n_fof_groups_kept = 0;

    size_t dst_size = sizeof(tree_entry_t);
    size_t dst_offsets[12] = {
        HOFFSET(tree_entry_t, id),
        HOFFSET(tree_entry_t, flags),
        HOFFSET(tree_entry_t, desc_id),
        HOFFSET(tree_entry_t, tree_id),
        HOFFSET(tree_entry_t, file_offset),
        HOFFSET(tree_entry_t, desc_index),
        HOFFSET(tree_entry_t, n_particle_peak),
        HOFFSET(tree_entry_t, central_index),
        HOFFSET(tree_entry_t, forest_id),
        HOFFSET(tree_entry_t, group_index),
        HOFFSET(tree_entry_t, original_index),
        HOFFSET(tree_entry_t, position_offset)
    };
    size_t dst_sizes[12] = {
        sizeof(tree_buffer[0].id),
        sizeof(tree_buffer[0].flags),
        sizeof(tree_buffer[0].desc_id),
        sizeof(tree_buffer[0].tree_id),
        sizeof(tree_buffer[0].file_offset),
        sizeof(tree_buffer[0].desc_index),
        sizeof(tree_buffer[0].n_particle_peak),
        sizeof(tree_buffer[0].central_index),
        sizeof(tree_buffer[0].forest_id),
        sizeof(tree_buffer[0].group_index),
        sizeof(tree_buffer[0].original_index),
        sizeof(tree_buffer[0].position_offset)
    };

    // Calculate the maximum catalog buffer size we require to hold a single
    // buffer worth of subgroups.  This is necessary as the extended trees might 
    // use same halos in the catalog file
    // Calculate the maximum group buffer size we require to hold a single
    // buffer worth of subgroups.  This is necessary as subfind can produce
    // many groups with no subgroups in them and so the group buffer may need
    // to be larger than the subgroup buffer.
    if (run_globals.mpi_rank == 0){
        int tmp_size;
        while (n_read < n_halos) {
            if ((n_halos - n_read) >= buffer_size)
                n_to_read = buffer_size;
            else
                n_to_read = n_halos - n_read;

            H5TBread_records(fd, "trees", n_read, (hsize_t)n_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);
            n_read += n_to_read;

            tmp_size = tree_buffer[n_to_read - 1].original_index - tree_buffer[0].original_index + 1;
            if (tmp_size > halo_buffer_size)
                halo_buffer_size = tmp_size;
            tmp_size = tree_buffer[n_to_read - 1].group_index - tree_buffer[0].group_index + 1;
            if (tmp_size > group_buffer_size)
                group_buffer_size = tmp_size;
        }
    }
    MPI_Bcast(&halo_buffer_size, 1, MPI_INT, 0, run_globals.mpi_comm);
    mlog("Using catalog buffer size = %d", MLOG_MESG, halo_buffer_size);
    catalog_buffer = malloc(sizeof(catalog_halo_t) * halo_buffer_size);

    MPI_Bcast(&group_buffer_size, 1, MPI_INT, 0, run_globals.mpi_comm);
    mlog("Using group buffer size = %d", MLOG_MESG, group_buffer_size);
    group_buffer = malloc(sizeof(catalog_halo_t) * group_buffer_size);

    // Now actually do the buffered read...
    n_read = 0;
    n_to_read = 0;
    keep_flag = true;
    while (n_read < n_halos) {
        if ((n_halos - n_read) >= buffer_size)
            n_to_read = buffer_size;
        else
            n_to_read = n_halos - n_read;

        // read in a tree_buffer of the trees
        if (run_globals.mpi_rank == 0)
            H5TBread_records(fd, "trees", n_read, (hsize_t)n_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);
        MPI_Bcast(tree_buffer, n_to_read * sizeof(tree_entry_t), MPI_BYTE, 0, run_globals.mpi_comm);

        first_halo_index = tree_buffer[0].original_index;
        if (first_halo_index == last_halo_index) {
            i_halo = 1;
            memcpy(&(catalog_buffer[0]), &(catalog_buffer[n_catalogs - 1]), sizeof(catalog_halo_t));
        }
        else
            i_halo = 0;
        last_halo_index = tree_buffer[n_to_read - 1].original_index;
        n_catalogs = last_halo_index - first_halo_index + 1;

        first_group_index = tree_buffer[0].group_index;
        if (first_group_index == last_group_index) {
            i_group = 1;
            memcpy(&(group_buffer[0]), &(group_buffer[n_groups - 1]), sizeof(catalog_halo_t));
        }
        else
            i_group = 0;
        last_group_index = tree_buffer[n_to_read - 1].group_index;
        n_groups = last_group_index - first_group_index + 1;

        // read in the corresponding catalog entrys
        if (run_globals.mpi_rank == 0) {
            i_halo = 0;
            read_catalog_halos(&fin_catalogs, simulation_dir, catalog_file_prefix, snapshot,
                &flayout_switch, &i_catalog_file, &n_halos_in_catalog_file,
                &i_halo_in_catalog_file, &i_halo, catalog_buffer, n_catalogs - i_halo, 1);
            read_catalog_halos(&fin_groups, simulation_dir, catalog_file_prefix, snapshot,
                &group_flayout_switch, &i_group_file, &n_groups_in_catalog_file,
                &i_group_in_catalog_file, &i_group, group_buffer, n_groups - i_group, 0);
        }
        MPI_Bcast(catalog_buffer, n_catalogs * sizeof(catalog_halo_t), MPI_BYTE, 0, run_globals.mpi_comm);
        MPI_Bcast(group_buffer, n_groups * sizeof(catalog_halo_t), MPI_BYTE, 0, run_globals.mpi_comm);

        // paste the data into the halo structures
        for (int jj = 0; jj < n_to_read; jj++) {
            if (run_globals.RequestedForestId != NULL) {
                if (bsearch(&(tree_buffer[jj].forest_id), run_globals.RequestedForestId,
                        (size_t)n_requested_forests, sizeof(int), compare_ints)
                    != NULL)
                    keep_flag = true;
                else
                    keep_flag = false;
            }

            if (keep_flag) {
                tree_entry_t* cur_tree_entry = &(tree_buffer[jj]);

                assert((*n_halos_kept) < run_globals.NHalosMax);
                assert((tree_buffer[jj].original_index - first_halo_index) < n_catalogs);
                catalog_halo_t* cur_cat_halo = &(catalog_buffer[tree_buffer[jj].original_index - first_halo_index]);
                halo_t* cur_halo = &(halo[*n_halos_kept]);

                cur_halo->ID = cur_tree_entry->id;
                cur_halo->TreeFlags = cur_tree_entry->flags;
                cur_halo->SnapOffset = cur_tree_entry->file_offset;
                cur_halo->DescIndex = cur_tree_entry->desc_index;
                cur_halo->NextHaloInFOFGroup = NULL;

                if (index_lookup)
                    index_lookup[*n_halos_kept] = n_read + jj;

                if (n_read + jj == tree_buffer[jj].central_index) {
                    cur_halo->Type = 0;

                    assert((*n_fof_groups_kept) < run_globals.NFOFGroupsMax);
                    assert((tree_buffer[jj].group_index - first_group_index) < n_groups);
                    catalog_halo_t* cur_cat_group = &(group_buffer[tree_buffer[jj].group_index - first_group_index]);
                    fof_group_t* cur_group = &(fof_group[*n_fof_groups_kept]);

                    cur_group->Mvir = cur_cat_group->M_vir;
                    cur_group->Rvir = cur_cat_group->R_vir;
                    cur_group->FOFMvirModifier = 1.0;

                    convert_input_virial_props(&(cur_group->Mvir),
                        &(cur_group->Rvir),
                        &(cur_group->Vvir),
                        &(cur_group->FOFMvirModifier),
                        -1,
                        snapshot,
                        true);

                    fof_group[(*n_fof_groups_kept)++].FirstHalo = &(halo[*n_halos_kept]);
                }
                else {
                    cur_halo->Type = 1;
                    halo[(*n_halos_kept) - 1].NextHaloInFOFGroup = &(halo[*n_halos_kept]);
                }

                cur_halo->FOFGroup = &(fof_group[(*n_fof_groups_kept) - 1]);

                // paste in the halo properties
                cur_halo->Len = cur_cat_halo->n_particles;
                cur_halo->Pos[0] = cur_cat_halo->position_MBP[0] + cur_tree_entry->position_offset[0];
                cur_halo->Pos[1] = cur_cat_halo->position_MBP[1] + cur_tree_entry->position_offset[1];
                cur_halo->Pos[2] = cur_cat_halo->position_MBP[2] + cur_tree_entry->position_offset[2];
                cur_halo->Vel[0] = cur_cat_halo->velocity_COM[0];
                cur_halo->Vel[1] = cur_cat_halo->velocity_COM[1];
                cur_halo->Vel[2] = cur_cat_halo->velocity_COM[2];
                cur_halo->Rvir = cur_cat_halo->R_vir;
                cur_halo->Vmax = cur_cat_halo->V_max;
                cur_halo->AngMom[0] = cur_cat_halo->ang_mom[0];
                cur_halo->AngMom[1] = cur_cat_halo->ang_mom[1];
                cur_halo->AngMom[2] = cur_cat_halo->ang_mom[2];
                cur_halo->Galaxy = NULL;
                cur_halo->Mvir = cur_cat_halo->M_vir;

                // double check that PBC conditions are met!
                cur_halo->Pos[0] = apply_pbc_pos(cur_halo->Pos[0]);
                cur_halo->Pos[1] = apply_pbc_pos(cur_halo->Pos[1]);
                cur_halo->Pos[2] = apply_pbc_pos(cur_halo->Pos[2]);

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
                    snapshot,
                    check_for_flag(TREE_CASE_HIGHRES & TREE_CASE_TRAITOR, cur_halo->TreeFlags));

                (*n_halos_kept)++;
            }
        }

        n_read += n_to_read;
    }

    // close the catalogs and trees files
    if (run_globals.mpi_rank == 0) {
        H5Fclose(fd);

        if (fin_catalogs)
            fclose(fin_catalogs);
    }

    // free the buffers
    free(group_buffer);
    free(catalog_buffer);
    free(tree_buffer);

    mlog(" ...done", MLOG_CLOSE);
}
