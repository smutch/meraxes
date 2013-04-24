#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>


void calc_hdf5_props(run_globals_struct *run_globals)
{

  /*
   * Prepare an HDF5 to receive the output galaxy data.
   * Here we store the data in an hdf5 table for easily appending new data.
   */

  hdf5_output_struct     *h5props = &(run_globals->hdf5props);
  galaxy_output_struct   galout;
  hid_t                  array_nmag_f_tid;

	int                    i;  // dummy

  // If we are calculating any magnitudes then increment the number of
  // output properties appropriately.
  h5props->n_props = 18;  // not inc. magnitudes

  // Size of a single galaxy entry.
  h5props->dst_size = sizeof(galaxy_output_struct);

  // Create datatypes for different size arrays
  hid_t array3f_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){3});

  // Calculate the offsets of our struct members in memory
  h5props->dst_offsets     = SID_malloc(sizeof(size_t)*h5props->n_props);
  // Calculate the sizes of our struct members in memory.
  h5props->dst_sizes       = SID_malloc(sizeof(size_t)*h5props->n_props);
  // Give each galaxy property a field name in the table
  h5props->field_names     = SID_malloc(sizeof(const char*)*h5props->n_props);
  // Assign a type to each galaxy property field in the table.
  h5props->field_types     = SID_malloc(sizeof(hid_t)*h5props->n_props);

  i=0;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Type);
  h5props->dst_sizes[i]    = sizeof(galout.Type);
  h5props->field_names[i]  = "Type";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, HaloIndex);
  h5props->dst_sizes[i]    = sizeof(galout.HaloIndex);
  h5props->field_names[i]  = "HaloIndex";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, CentralGal);
  h5props->dst_sizes[i]    = sizeof(galout.CentralGal);
  h5props->field_names[i]  = "CentralGal";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Pos);
  h5props->dst_sizes[i]    = sizeof(galout.Pos);
  h5props->field_names[i]  = "Pos";
  h5props->field_types[i++]  = array3f_tid;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Vel);
  h5props->dst_sizes[i]    = sizeof(galout.Vel);
  h5props->field_names[i]  = "Vel";
  h5props->field_types[i++]  = array3f_tid;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Spin);
  h5props->dst_sizes[i]    = sizeof(galout.Spin);
  h5props->field_names[i]  = "Spin";
  h5props->field_types[i++]  = array3f_tid;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Len);
  h5props->dst_sizes[i]   = sizeof(galout.Len);
  h5props->field_names[i] = "Len";
  h5props->field_types[i++] = H5T_NATIVE_INT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Mvir);
  h5props->dst_sizes[i]   = sizeof(galout.Mvir);
  h5props->field_names[i] = "Mvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, dM);
  h5props->dst_sizes[i]   = sizeof(galout.dM);
  h5props->field_names[i] = "dM";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, dMdt);
  h5props->dst_sizes[i]   = sizeof(galout.dMdt);
  h5props->field_names[i] = "dMdt";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Rvir);
  h5props->dst_sizes[i]   = sizeof(galout.Rvir);
  h5props->field_names[i] = "Rvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Vvir);
  h5props->dst_sizes[i]   = sizeof(galout.Vvir);
  h5props->field_names[i] = "Vvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Vmax);
  h5props->dst_sizes[i]   = sizeof(galout.Vmax);
  h5props->field_names[i] = "Vmax";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, VelDisp);
  h5props->dst_sizes[i]   = sizeof(galout.VelDisp);
  h5props->field_names[i] = "VelDisp";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, StellarMass);
  h5props->dst_sizes[i]   = sizeof(galout.StellarMass);
  h5props->field_names[i] = "StellarMass";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Sfr);
  h5props->dst_sizes[i]   = sizeof(galout.Sfr);
  h5props->field_names[i] = "Sfr";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Cos_Inc);
  h5props->dst_sizes[i]   = sizeof(galout.Cos_Inc);
  h5props->field_names[i] = "Cos_Inc";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, MergTime);
  h5props->dst_sizes[i]   = sizeof(galout.MergTime);
  h5props->field_names[i] = "MergTime";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  // DEBUG
  if(i != h5props->n_props){
    SID_log_error("Incorrect number of galaxy properties in HDF5 file.");
    ABORT(EXIT_FAILURE);
  }

}


void prep_hdf5_file(run_globals_struct *run_globals, char fname[STRLEN])
{

  hsize_t  chunk_size        = 10;
  int     *fill_data         = NULL;
  hid_t    file_id;
  hid_t    file_group_id;
  hid_t    snap_group_id;
  char     target_group[100];
  int      filenr;                   // dummy
  hid_t    status;
  hdf5_output_struct h5props = run_globals->hdf5props;

  // Create a new file
  file_id = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  // Create a group for each output file and snapshot
  for (filenr = run_globals->params.FirstFile; filenr <= run_globals->params.LastFile; filenr++) {

    for (int i_snap=0; i_snap<NOUT; i_snap++) {

      sprintf(target_group, "Snap%03d", run_globals->ListOutputSnaps[i_snap]);
      snap_group_id = H5Gcreate(file_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Make the table
      status=H5TBmake_table( "Galaxy Table", snap_group_id, "Galaxies",
          h5props.n_props,0, h5props.dst_size, h5props.field_names,
          h5props.dst_offsets, h5props.field_types, chunk_size, fill_data, 0,
          NULL );

      H5Gclose(snap_group_id);

    }

  }

  // Close the HDF5 file.
  status = H5Fclose(file_id);

}


void write_galaxy(run_globals_struct *run_globals, galaxy_struct *gal, int i_out, char fname[STRLEN])
{

  /*
   * Write a single galaxy to the hdf5 file table.
   */

  herr_t status;
  hid_t  file_id, group_id;
  char   target_group[100];
  galaxy_output_struct galout;
  hdf5_output_struct h5props = run_globals->hdf5props;

  // Open the file.
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  // Open the relevant group.
  sprintf(target_group, "Snap%03d", run_globals->ListOutputSnaps[i_out]);
  group_id = H5Gopen(file_id, target_group, H5P_DEFAULT);

  // Write the galaxy.
  prepare_galaxy_for_output(gal, &galout, i_out);
  status = H5TBappend_records(group_id, "Galaxies", 1, h5props.dst_size, h5props.dst_offsets, h5props.dst_sizes,
      &galout);

  // Close the group
  status = H5Gclose(group_id);

  // Close the file.
  status = H5Fclose(file_id);

}


void write_snapshot(run_globals_struct *run_globals, int NGal, int i_out, char fname[STRLEN])
{

  /*
   * Write a batch of galaxies to the output HDF5 table.
   */

  herr_t status;
  hid_t  file_id;
  hid_t  group_id;
  char   target_group[100];
  galaxy_output_struct galout;
  galaxy_struct *gal = NULL;
  hdf5_output_struct h5props = run_globals->hdf5props;

  // Open the file.
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  // Open the relevant group.
  sprintf(target_group, "Snap%03d", run_globals->ListOutputSnaps[i_out]);
  group_id = H5Gopen(file_id, target_group, H5P_DEFAULT);

  // Write the galaxies.
  gal = run_globals->FirstGal; 
  do {
    prepare_galaxy_for_output(gal, &galout, i_out);

    status=H5TBappend_records(group_id, "Galaxies", 1, h5props.dst_size,
        h5props.dst_offsets, h5props.dst_sizes, &galout);

    gal = gal->Next;
  } while (gal!=NULL);

  // Close the group.
  status = H5Gclose(group_id);

  // Close the file.
  status = H5Fclose(file_id);

}


void write_hdf5_attrs(run_globals_struct *run_globals, int NGal, int i_out, char fname[STRLEN])
{

  /*
   * Write the HDF5 file attributes.
   */

  herr_t  status;
  hid_t   file_id;
  hid_t   dataset_id;
  hid_t   attribute_id;
  hid_t   dataspace_id;
  hid_t   group_id;
  hsize_t dims;
  char    target_group[100];

  // Open the output file and galaxy dataset.
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  // Open the relevant group.
  sprintf(target_group, "Snap%03d", run_globals->ListOutputSnaps[i_out]);
  group_id = H5Gopen(file_id, target_group, H5P_DEFAULT);

  dataset_id = H5Dopen(group_id, "Galaxies", H5P_DEFAULT);

  // Create the data space for the attributes.
  dims = 1;
  dataspace_id = H5Screate_simple(1, &dims, NULL);

  // Write the total number of galaxies.
  attribute_id = H5Acreate(dataset_id, "NGalaxies", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_INT, &NGal);
  status = H5Aclose(attribute_id);

  // Close the dataspace.
  status = H5Sclose(dataspace_id);

  // Close to the dataset.
  status = H5Dclose(dataset_id);

  // Close the group.
  status = H5Gclose(group_id);

  // Close the file.
  status = H5Fclose(file_id);

}


