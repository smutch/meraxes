#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>


void calc_hdf5_props()
{

  /*
   * Prepare an HDF5 to receive the output galaxy data.
   * Here we store the data in an hdf5 table for easily appending new data.
   */

  galaxy_output_struct   galout;
  hid_t                  array_nmag_f_tid;

	int                    i;  // dummy

  // If we are calculating any magnitudes then increment the number of
  // output properties appropriately.
  HDF5_n_props = 18;  // not inc. magnitudes

  // Size of a single galaxy entry.
  HDF5_dst_size = sizeof(galaxy_output_struct);

  // Create datatypes for different size arrays
  hid_t array3f_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){3});

  // Calculate the offsets of our struct members in memory
  HDF5_dst_offsets     = mymalloc(sizeof(size_t)*HDF5_n_props);
  // Calculate the sizes of our struct members in memory.
  HDF5_dst_sizes       = mymalloc(sizeof(size_t)*HDF5_n_props);
  // Give each galaxy property a field name in the table
  HDF5_field_names     = mymalloc(sizeof(const char*)*HDF5_n_props);
  // Assign a type to each galaxy property field in the table.
  HDF5_field_types     = mymalloc(sizeof(hid_t)*HDF5_n_props);

  i=0;

  HDF5_dst_offsets[i]  = HOFFSET(struct GALAXY_OUTPUT, Type);
  HDF5_dst_sizes[i]    = sizeof(galout.Type);
  HDF5_field_names[i]  = "Type";
  HDF5_field_types[i++]  = H5T_NATIVE_INT;

  HDF5_dst_offsets[i]  = HOFFSET(struct GALAXY_OUTPUT, HaloIndex);
  HDF5_dst_sizes[i]    = sizeof(galout.HaloIndex);
  HDF5_field_names[i]  = "HaloIndex";
  HDF5_field_types[i++]  = H5T_NATIVE_INT;

  HDF5_dst_offsets[i]  = HOFFSET(struct GALAXY_OUTPUT, CentralGal);
  HDF5_dst_sizes[i]    = sizeof(galout.CentralGal);
  HDF5_field_names[i]  = "CentralGal";
  HDF5_field_types[i++]  = H5T_NATIVE_INT;

  HDF5_dst_offsets[i]  = HOFFSET(struct GALAXY_OUTPUT, Pos);
  HDF5_dst_sizes[i]    = sizeof(galout.Pos);
  HDF5_field_names[i]  = "Pos";
  HDF5_field_types[i++]  = array3f_tid;

  HDF5_dst_offsets[i]  = HOFFSET(struct GALAXY_OUTPUT, Vel);
  HDF5_dst_sizes[i]    = sizeof(galout.Vel);
  HDF5_field_names[i]  = "Vel";
  HDF5_field_types[i++]  = array3f_tid;

  HDF5_dst_offsets[i]  = HOFFSET(struct GALAXY_OUTPUT, Spin);
  HDF5_dst_sizes[i]    = sizeof(galout.Spin);
  HDF5_field_names[i]  = "Spin";
  HDF5_field_types[i++]  = array3f_tid;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Len);
  HDF5_dst_sizes[i]   = sizeof(galout.Len);
  HDF5_field_names[i] = "Len";
  HDF5_field_types[i++] = H5T_NATIVE_INT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Mvir);
  HDF5_dst_sizes[i]   = sizeof(galout.Mvir);
  HDF5_field_names[i] = "Mvir";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, dM);
  HDF5_dst_sizes[i]   = sizeof(galout.dM);
  HDF5_field_names[i] = "dM";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, dMdt);
  HDF5_dst_sizes[i]   = sizeof(galout.dMdt);
  HDF5_field_names[i] = "dMdt";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Rvir);
  HDF5_dst_sizes[i]   = sizeof(galout.Rvir);
  HDF5_field_names[i] = "Rvir";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Vvir);
  HDF5_dst_sizes[i]   = sizeof(galout.Vvir);
  HDF5_field_names[i] = "Vvir";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Vmax);
  HDF5_dst_sizes[i]   = sizeof(galout.Vmax);
  HDF5_field_names[i] = "Vmax";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, VelDisp);
  HDF5_dst_sizes[i]   = sizeof(galout.VelDisp);
  HDF5_field_names[i] = "VelDisp";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, StellarMass);
  HDF5_dst_sizes[i]   = sizeof(galout.StellarMass);
  HDF5_field_names[i] = "StellarMass";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Sfr);
  HDF5_dst_sizes[i]   = sizeof(galout.Sfr);
  HDF5_field_names[i] = "Sfr";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, Cos_Inc);
  HDF5_dst_sizes[i]   = sizeof(galout.Cos_Inc);
  HDF5_field_names[i] = "Cos_Inc";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  HDF5_dst_offsets[i] = HOFFSET(struct GALAXY_OUTPUT, MergTime);
  HDF5_dst_sizes[i]   = sizeof(galout.MergTime);
  HDF5_field_names[i] = "MergTime";
  HDF5_field_types[i++] = H5T_NATIVE_FLOAT;

  // DEBUG
  if(i != HDF5_n_props){
    SID_log_error("Incorrect number of galaxy properties in HDF5 file.");
    ABORT(EXIT_FAILURE);
  }

}


void prep_hdf5_file(char *fname)
{

  hsize_t  chunk_size        = 10;
  int     *fill_data         = NULL;
  hid_t    file_id;
  hid_t    file_group_id;
  hid_t    snap_group_id;
  char     target_group[100];
  int      filenr;                   // dummy
  hid_t    status;

  // Create a new file
  file_id = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  // Create a group for each output file and snapshot
  for (filenr = FirstFile; filenr <= LastFile; filenr++) {

    for (int i_snap=0; i_snap<NOUT; i_snap++) {

      sprintf(target_group, "Snap%03d", ListOutputSnaps[i_snap]);
      snap_group_id = H5Gcreate(file_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Make the table
      status=H5TBmake_table( "Galaxy Table", snap_group_id, "Galaxies", HDF5_n_props,0,
          HDF5_dst_size,HDF5_field_names, HDF5_dst_offsets, HDF5_field_types,
          chunk_size, fill_data, 0, NULL );

      H5Gclose(snap_group_id);

    }

  }

  // Close the HDF5 file.
  status = H5Fclose(file_id);

}


void write_galaxy(galaxy_struct *gal, int i_out)
{

  /*
   * Write a single galaxy to the hdf5 file table.
   */

  herr_t status;
  hid_t  file_id, group_id;
  char   target_group[100];
  char   fname[1000];
  galaxy_output_struct galout;

  // Generate the filename to be opened.
  sprintf(fname, "%s/%s.hdf5", OutputDir, FileNameGalaxies);

  // Open the file.
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  // Open the relevant group.
  sprintf(target_group, "Snap%03d", ListOutputSnaps[i_out]);
  group_id = H5Gopen(file_id, target_group, H5P_DEFAULT);

  // Write the galaxy.
  prepare_galaxy_for_output(gal, &galout, i_out);
  status = H5TBappend_records(group_id, "Galaxies", 1, HDF5_dst_size, HDF5_dst_offsets, HDF5_dst_sizes,
      galout);

  // Close the group
  status = H5Gclose(group_id);

  // Close the file.
  status = H5Fclose(file_id);

}


void write_snapshot(run_globals_struct *run_globals, int NGal, int i_out)
{

  /*
   * Write a batch of galaxies to the output HDF5 table.
   */

  herr_t status;
  hid_t  file_id;
  hid_t  group_id;
  char   target_group[100];
  char   fname[1000];
  galaxy_output_struct galout;
  galaxy_struct *gal = NULL;

  // Generate the filename to be opened.
  sprintf(fname, "%s/%s.hdf5", OutputDir, FileNameGalaxies, filenr);

  // Open the file.
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  // Open the relevant group.
  sprintf(target_group, "Snap%03d", ListOutputSnaps[i_out]);
  group_id = H5Gopen(file_id, target_group, H5P_DEFAULT);

  // Write the galaxies.
  gal = run_globals->FirstGal; 
  do {
    prepare_galaxy_for_output(gal, &galout, i_out);

    status=H5TBappend_records(group_id, "Galaxies", 1, HDF5_dst_size,
        HDF5_dst_offsets, HDF5_dst_sizes, galout);

    gal = gal->Next;
  } while (gal!=NULL);

  // Close the group.
  status = H5Gclose(group_id);

  // Close the file.
  status = H5Fclose(file_id);

}


void write_hdf5_attrs(int NGal, int i_out, int filenr)
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
  char    fname[1000];

  // Generate the filename to be opened.
  sprintf(fname, "%s/%s.hdf5", OutputDir, FileNameGalaxies, filenr);

  // Open the output file and galaxy dataset.
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

  // Open the relevant group.
  sprintf(target_group, "Snap%03d", ListOutputSnaps[i_out]);
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


