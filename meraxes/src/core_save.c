#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>

void prepare_galaxy_for_output(
  run_globals_struct   *run_globals,
  galaxy_struct         gal,        
  galaxy_output_struct *galout,     
  int                   i_snap)     
{

  double Hubble_h = run_globals->params.Hubble_h;
  run_units_struct *units = &(run_globals->units);

  galout->Type = (int)gal.Type;

  for(int ii=0; ii<3; ii++)
  {
    galout->Pos[ii]  = (float)(gal.Pos[ii] / Hubble_h);
    galout->Vel[ii]  = (float)(gal.Vel[ii]);
  }

  galout->Len         = (int)(gal.Len);
  galout->Mvir        = (float)(gal.Mvir / Hubble_h);
  galout->dM          = (float)(gal.dM / Hubble_h);
  galout->dMdt        = (float)(gal.dMdt * units->UnitMass_in_g / units->UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  galout->Rvir        = (float)(gal.Rvir / Hubble_h);
  galout->Vvir        = (float)(gal.Vvir);
  galout->Vmax        = (float)(gal.Vmax);
  galout->StellarMass = (float)(gal.StellarMass / 1.e10 / Hubble_h);
  galout->Sfr         = (float)(gal.Sfr[i_snap] * units->UnitMass_in_g / units->UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  galout->Cos_Inc     = (float)(gal.Cos_Inc);
  galout->MergTime    = (float)(gal.MergTime * units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR / Hubble_h);

}

void calc_hdf5_props(run_globals_struct *run_globals)
{

  /*
   * Prepare an HDF5 to receive the output galaxy data.
   * Here we store the data in an hdf5 table for easily appending new data.
   */

  hdf5_output_struct     *h5props = &(run_globals->hdf5props);
  galaxy_output_struct   galout;

	int                    i;  // dummy

  // If we are calculating any magnitudes then increment the number of
  // output properties appropriately.
  h5props->n_props = 14;  // not inc. magnitudes

  // Size of a single galaxy entry.
  h5props->dst_size = sizeof(galaxy_output_struct);

  // Create datatypes for different size arrays
  hid_t array3f_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){3});

  // Calculate the offsets of our struct members in memory
  h5props->dst_offsets     = SID_malloc(sizeof(size_t)*h5props->n_props);
  // Calculate the sizes of our struct members in memory.
  h5props->dst_field_sizes = SID_malloc(sizeof(size_t)*h5props->n_props);
  // Give each galaxy property a field name in the table
  h5props->field_names     = SID_malloc(sizeof(const char*)*h5props->n_props);
  // Assign a type to each galaxy property field in the table.
  h5props->field_types     = SID_malloc(sizeof(hid_t)*h5props->n_props);

  i=0;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Type);
  h5props->dst_field_sizes[i]    = sizeof(galout.Type);
  h5props->field_names[i]  = "Type";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Pos);
  h5props->dst_field_sizes[i]    = sizeof(galout.Pos);
  h5props->field_names[i]  = "Pos";
  h5props->field_types[i++]  = array3f_tid;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_struct, Vel);
  h5props->dst_field_sizes[i]    = sizeof(galout.Vel);
  h5props->field_names[i]  = "Vel";
  h5props->field_types[i++]  = array3f_tid;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Len);
  h5props->dst_field_sizes[i]   = sizeof(galout.Len);
  h5props->field_names[i] = "Len";
  h5props->field_types[i++] = H5T_NATIVE_INT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Mvir);
  h5props->dst_field_sizes[i]   = sizeof(galout.Mvir);
  h5props->field_names[i] = "Mvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, dM);
  h5props->dst_field_sizes[i]   = sizeof(galout.dM);
  h5props->field_names[i] = "dM";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, dMdt);
  h5props->dst_field_sizes[i]   = sizeof(galout.dMdt);
  h5props->field_names[i] = "dMdt";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Rvir);
  h5props->dst_field_sizes[i]   = sizeof(galout.Rvir);
  h5props->field_names[i] = "Rvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Vvir);
  h5props->dst_field_sizes[i]   = sizeof(galout.Vvir);
  h5props->field_names[i] = "Vvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Vmax);
  h5props->dst_field_sizes[i]   = sizeof(galout.Vmax);
  h5props->field_names[i] = "Vmax";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, StellarMass);
  h5props->dst_field_sizes[i]   = sizeof(galout.StellarMass);
  h5props->field_names[i] = "StellarMass";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Sfr);
  h5props->dst_field_sizes[i]   = sizeof(galout.Sfr);
  h5props->field_names[i] = "Sfr";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, Cos_Inc);
  h5props->dst_field_sizes[i]   = sizeof(galout.Cos_Inc);
  h5props->field_names[i] = "Cos_Inc";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_struct, MergTime);
  h5props->dst_field_sizes[i]   = sizeof(galout.MergTime);
  h5props->field_names[i] = "MergTime";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  // DEBUG
  if(i != h5props->n_props){
    SID_log_error("Incorrect number of galaxy properties in HDF5 file.");
    ABORT(EXIT_FAILURE);
  }

}


void prep_hdf5_file(run_globals_struct *run_globals)
{

  hid_t    file_id;
  hid_t    status;

  // Create a new file
  file_id = H5Fcreate(run_globals->FNameOut, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  // Close the HDF5 file.
  status = H5Fclose(file_id);

}


void write_snapshot(run_globals_struct *run_globals, int n_write, int i_out)
{

  /*
   * Write a batch of galaxies to the output HDF5 table.
   */

  herr_t status;
  hid_t  file_id;
  hid_t  group_id;
  hsize_t  chunk_size        = 10;
  int     *fill_data         = NULL;
  char   target_group[20];
  galaxy_output_struct galout;
  galaxy_struct *gal = NULL;
  hdf5_output_struct h5props = run_globals->hdf5props;
  int gal_count = 0;

  // Create the file.
  file_id = H5Fopen(run_globals->FNameOut, H5F_ACC_RDWR, H5P_DEFAULT);

  // Create the relevant group.
  sprintf(target_group, "Snap%03d", (run_globals->ListOutputSnaps)[i_out]);
  group_id = H5Gcreate(file_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Make the table
  status=H5TBmake_table("Galaxies", group_id, "Galaxies",
      h5props.n_props, n_write, h5props.dst_size, h5props.field_names,
      h5props.dst_offsets, h5props.field_types, chunk_size, fill_data, 0,
      NULL );

  // // Write the galaxies.
  gal = run_globals->FirstGal; 
  while (gal!=NULL) {
    // Don't output galaxies which merged at this timestep
    if (gal->Type < 3)
    {
      prepare_galaxy_for_output(run_globals, *gal, &galout, i_out);
      status=H5TBwrite_records(group_id, "Galaxies", gal_count, 1, h5props.dst_size,
          h5props.dst_offsets, h5props.dst_field_sizes, &galout);
      gal_count++;
    }

    gal = gal->Next;
  }
 
  // Close the group.
  status = H5Gclose(group_id);

  // Close the file.
  status = H5Fclose(file_id);

}


