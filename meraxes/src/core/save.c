#include "meraxes.h"
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>

static void inline h5_write_attribute(hid_t loc, const char *name, hid_t datatype, hid_t dataset_id, void *data)
{
  herr_t status;
  hid_t attr_id;

  attr_id = H5Acreate(loc, name, datatype, dataset_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_id, datatype, data);
  if (status<0)
  {
    SID_log("Error writing attribute '%s'", SID_LOG_COMMENT, name);
    ABORT(EXIT_FAILURE);
  }
  status = H5Aclose(attr_id);
  if (status<0)
  {
    SID_log("Error closing attribute '%s'", SID_LOG_COMMENT, name);
    ABORT(EXIT_FAILURE);
  }
}

void prepare_galaxy_for_output(
  run_globals_t   *run_globals,
  galaxy_t         gal,
  galaxy_output_t *galout,
  int                   i_snap)
{

  double Hubble_h = run_globals->params.Hubble_h;
  run_units_t *units = &(run_globals->units);

  galout->id_MBP = (long long)gal.id_MBP;
  galout->ID = (int)gal.ID;
  galout->Type = (int)gal.Type;
  if((!gal.ghost_flag) && (gal.Halo->FOFGroup->FirstHalo->Galaxy!=NULL))
    galout->CentralGal = (int)gal.Halo->FOFGroup->FirstHalo->Galaxy->output_index;
  else
    galout->CentralGal = -1;
  galout->GhostFlag = (int)gal.ghost_flag;

  for(int ii=0; ii<3; ii++)
  {
    galout->Pos[ii]  = (float)(gal.Pos[ii] / Hubble_h);
    galout->Vel[ii]  = (float)(gal.Vel[ii]);
  }

  galout->Len         = (int)(gal.Len);
  galout->Mvir        = (float)(gal.Mvir / Hubble_h);
  galout->dM          = (float)(gal.dM / Hubble_h);
  galout->dMdt        = (float)(gal.dM/gal.dt * units->UnitMass_in_g / units->UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  galout->Rvir        = (float)(gal.Rvir / Hubble_h);
  galout->Vvir        = (float)(gal.Vvir);
  galout->Vmax        = (float)(gal.Vmax);
  galout->StellarMass = (float)(gal.StellarMass / Hubble_h);
  galout->Sfr         = (float)(gal.Sfr * units->UnitMass_in_g / units->UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  galout->Cos_Inc     = (float)(gal.Cos_Inc);
  galout->MergTime    = (float)(gal.MergTime * units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR / Hubble_h);
  galout->LTTime      = (float)(gal.LTTime * units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR / Hubble_h);

  prepare_magnitudes_for_output(run_globals, gal, galout, i_snap);

}

void calc_hdf5_props(run_globals_t *run_globals)
{

  /*
   * Prepare an HDF5 to receive the output galaxy data.
   * Here we store the data in an hdf5 table for easily appending new data.
   */

  hdf5_output_t   *h5props = &(run_globals->hdf5props);
  galaxy_output_t  galout;
  int              i;                                   // dummy

  // If we are calculating any magnitudes then increment the number of
  // output properties appropriately.
  h5props->n_props = 19;

#ifdef CALC_MAGS
  int n_photo_bands = run_globals->photo.NBands;
  h5props->n_props +=2;
  h5props->array_nmag_f_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){n_photo_bands});
#endif

  // Size of a single galaxy entry.
  h5props->dst_size = sizeof(galaxy_output_t);

  // Create datatypes for different size arrays
  h5props->array3f_tid      = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){3});

  // Calculate the offsets of our struct members in memory
  h5props->dst_offsets     = SID_malloc(sizeof(size_t)*h5props->n_props);
  // Calculate the sizes of our struct members in memory.
  h5props->dst_field_sizes = SID_malloc(sizeof(size_t)*h5props->n_props);
  // Give each galaxy property a field name in the table
  h5props->field_names     = SID_malloc(sizeof(const char*)*h5props->n_props);
  // Assign a type to each galaxy property field in the table.
  h5props->field_types     = SID_malloc(sizeof(hid_t)*h5props->n_props);

  i=0;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_t, id_MBP);
  h5props->dst_field_sizes[i]    = sizeof(galout.id_MBP);
  h5props->field_names[i]  = "id_MBP";
  h5props->field_types[i++]  = H5T_NATIVE_LLONG;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_t, ID);
  h5props->dst_field_sizes[i]    = sizeof(galout.ID);
  h5props->field_names[i]  = "ID";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_t, Type);
  h5props->dst_field_sizes[i]    = sizeof(galout.Type);
  h5props->field_names[i]  = "Type";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_t, CentralGal);
  h5props->dst_field_sizes[i]    = sizeof(galout.CentralGal);
  h5props->field_names[i]  = "CentralGal";
  h5props->field_types[i++]  = H5T_NATIVE_INT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, GhostFlag);
  h5props->dst_field_sizes[i]   = sizeof(galout.GhostFlag);
  h5props->field_names[i] = "GhostFlag";
  h5props->field_types[i++] = H5T_NATIVE_INT;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_t, Pos);
  h5props->dst_field_sizes[i]    = sizeof(galout.Pos);
  h5props->field_names[i]  = "Pos";
  h5props->field_types[i++]  = h5props->array3f_tid;

  h5props->dst_offsets[i]  = HOFFSET(galaxy_output_t, Vel);
  h5props->dst_field_sizes[i]    = sizeof(galout.Vel);
  h5props->field_names[i]  = "Vel";
  h5props->field_types[i++]  = h5props->array3f_tid;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Len);
  h5props->dst_field_sizes[i]   = sizeof(galout.Len);
  h5props->field_names[i] = "Len";
  h5props->field_types[i++] = H5T_NATIVE_INT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Mvir);
  h5props->dst_field_sizes[i]   = sizeof(galout.Mvir);
  h5props->field_names[i] = "Mvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, dM);
  h5props->dst_field_sizes[i]   = sizeof(galout.dM);
  h5props->field_names[i] = "dM";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, dMdt);
  h5props->dst_field_sizes[i]   = sizeof(galout.dMdt);
  h5props->field_names[i] = "dMdt";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Rvir);
  h5props->dst_field_sizes[i]   = sizeof(galout.Rvir);
  h5props->field_names[i] = "Rvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Vvir);
  h5props->dst_field_sizes[i]   = sizeof(galout.Vvir);
  h5props->field_names[i] = "Vvir";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Vmax);
  h5props->dst_field_sizes[i]   = sizeof(galout.Vmax);
  h5props->field_names[i] = "Vmax";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, StellarMass);
  h5props->dst_field_sizes[i]   = sizeof(galout.StellarMass);
  h5props->field_names[i] = "StellarMass";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Sfr);
  h5props->dst_field_sizes[i]   = sizeof(galout.Sfr);
  h5props->field_names[i] = "Sfr";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Cos_Inc);
  h5props->dst_field_sizes[i]   = sizeof(galout.Cos_Inc);
  h5props->field_names[i] = "Cos_Inc";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, MergTime);
  h5props->dst_field_sizes[i]   = sizeof(galout.MergTime);
  h5props->field_names[i] = "MergTime";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, LTTime);
  h5props->dst_field_sizes[i]   = sizeof(galout.LTTime);
  h5props->field_names[i] = "LTTime";
  h5props->field_types[i++] = H5T_NATIVE_FLOAT;

#ifdef CALC_MAGS
  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, Mag);
  h5props->dst_field_sizes[i]   = sizeof(float) * n_photo_bands;
  h5props->field_names[i] = "Mag";
  h5props->field_types[i++] = h5props->array_nmag_f_tid;

  h5props->dst_offsets[i] = HOFFSET(galaxy_output_t, MagDust);
  h5props->dst_field_sizes[i]   = sizeof(float) * n_photo_bands;
  h5props->field_names[i] = "MagDust";
  h5props->field_types[i++] = h5props->array_nmag_f_tid;
#endif

  // DEBUG
  if(i != h5props->n_props){
    SID_log_error("Incorrect number of galaxy properties in HDF5 file.");
    ABORT(EXIT_FAILURE);
  }

}


void prep_hdf5_file(run_globals_t *run_globals)
{

  hid_t    file_id, str_t, ds_id, group_id;
  hsize_t  dims = 1;
  char     names[50][STRLEN];
  void    *addresses[50];
  int ii;

  // Create a new file
  file_id = H5Fcreate(run_globals->FNameOut, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  // Set up reusable dataspaces and types
  ds_id = H5Screate_simple(1, &dims, NULL);
  str_t = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_t, STRLEN);

  // Open the group
  group_id = H5Gcreate(file_id, "InputParams", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  ii=0;
  addresses[ii] = &(run_globals->params.OutputDir);
  sprintf(names[ii++],  "OutputDir");
  addresses[ii] = &(run_globals->params.FileNameGalaxies);
  sprintf(names[ii++],  "FileNameGalaxies");
  addresses[ii] = &(run_globals->params.SimName);
  sprintf(names[ii++],  "SimName");
  addresses[ii] = &(run_globals->params.SimulationDir);
  sprintf(names[ii++],  "SimulationDir");
  addresses[ii] = &(run_globals->params.FileWithOutputSnaps);
  sprintf(names[ii++],  "FileWithOutputSnaps");
  if(run_globals->NRequestedForests > 0)
  {
    addresses[ii] = &(run_globals->params.ForestIDFile);
    sprintf(names[ii++], "ForestIDFile");
  }

  for(int jj=0; jj<ii; jj++)
    h5_write_attribute(group_id, names[jj], str_t, ds_id, addresses[jj]);

  ii=0;
  addresses[ii] = &(run_globals->params.NEverySnap);
  sprintf(names[ii++],  "NEverySnap");
  addresses[ii] = &(run_globals->params.NScanSnap);
  sprintf(names[ii++],  "NScanSnap");
  addresses[ii] = &(run_globals->params.TotalSimSnaps);
  sprintf(names[ii++],  "TotalSimSnaps");
  addresses[ii] = &(run_globals->params.LastSnapshotNr);
  sprintf(names[ii++],  "LastSnapshotNr");
  addresses[ii] = &(run_globals->params.FirstFile);
  sprintf(names[ii++],  "FirstFile");
  addresses[ii] = &(run_globals->params.LastFile);
  sprintf(names[ii++],  "LastFile");
  addresses[ii] = &(run_globals->params.SnaplistLength);
  sprintf(names[ii++],  "SnaplistLength");
  addresses[ii] = &(run_globals->NRequestedForests);
  sprintf(names[ii++], "NRequestedForests");
  if(SID.n_proc>0)
  {
    addresses[ii] = &(SID.n_proc);
    sprintf(names[ii++], "NOutputFiles");
  }

  for(int jj=0; jj<ii; jj++)
    h5_write_attribute(group_id, (const char *)(names[jj]), H5T_NATIVE_INT, ds_id, addresses[jj]);

  ii=0;
  addresses[ii] = &(run_globals->params.BoxSize);
  sprintf(names[ii++],  "BoxSize");
  addresses[ii] = &(run_globals->params.VolumeFactor);
  sprintf(names[ii++],  "VolumeFactor");
  addresses[ii] = &(run_globals->params.ThreshMajorMerger);
  sprintf(names[ii++],  "ThreshMajorMerger");
  addresses[ii] = &(run_globals->params.RecycleFraction);
  sprintf(names[ii++],  "RecycleFraction");
  addresses[ii] = &(run_globals->params.Hubble_h);
  sprintf(names[ii++],  "Hubble_h");
  addresses[ii] = &(run_globals->params.BaryonFrac);
  sprintf(names[ii++],  "BaryonFrac");
  addresses[ii] = &(run_globals->params.OmegaM);
  sprintf(names[ii++],  "OmegaM");
  addresses[ii] = &(run_globals->params.OmegaK);
  sprintf(names[ii++],  "OmegaK");
  addresses[ii] = &(run_globals->params.OmegaLambda);
  sprintf(names[ii++],  "OmegaLambda");
  addresses[ii] = &(run_globals->params.PartMass);
  sprintf(names[ii++],  "PartMass");
  addresses[ii] = &(run_globals->params.MergerTimeFactor);
  sprintf(names[ii++],  "MergerTimeFactor");

  for(int jj=0; jj<ii; jj++)
    h5_write_attribute(group_id, names[jj], H5T_NATIVE_DOUBLE, ds_id, addresses[jj]);

  // Close the group
  H5Gclose(group_id);

#ifdef CALC_MAGS
  // Open the group
  group_id = H5Gcreate(file_id, "InputParams/photo", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  ii=0;
  addresses[ii] = &(run_globals->params.SSPModel);
  sprintf(names[ii++],  "SSPModel");
  addresses[ii] = &(run_globals->params.IMF);
  sprintf(names[ii++],  "IMF");
  addresses[ii] = &(run_globals->params.MagSystem);
  sprintf(names[ii++],  "MagSystem");
  addresses[ii] = &(run_globals->params.MagBands);
  sprintf(names[ii++],  "MagBands");

  for(int jj=0; jj<ii; jj++)
    h5_write_attribute(group_id, names[jj], str_t, ds_id, addresses[jj]);

  // Close the group
  H5Gclose(group_id);
#endif

  // Open the group
  group_id = H5Gcreate(file_id, "InputParams/physics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  ii=0;
  addresses[ii] = &(run_globals->params.physics.peak);
  sprintf(names[ii++],  "peak");
  addresses[ii] = &(run_globals->params.physics.sigma);
  sprintf(names[ii++],  "sigma");
  addresses[ii] = &(run_globals->params.physics.stellarfrac);
  sprintf(names[ii++],  "stellarfrac");
  addresses[ii] = &(run_globals->params.physics.peak_evo);
  sprintf(names[ii++],  "peak_evo");
  addresses[ii] = &(run_globals->params.physics.sigma_evo);
  sprintf(names[ii++],  "sigma_evo");
  addresses[ii] = &(run_globals->params.physics.stellarfrac_evo);
  sprintf(names[ii++],  "stellarfrac_evo");
  addresses[ii] = &(run_globals->params.physics.bhgrowthfactor);
  sprintf(names[ii++],  "bhgrowthfactor");

  for(int jj=0; jj<ii; jj++)
    h5_write_attribute(group_id, names[jj], H5T_NATIVE_DOUBLE, ds_id, addresses[jj]);

  ii=0;
  addresses[ii] = &(run_globals->params.physics.funcprop);
  sprintf(names[ii++],  "funcprop");

  h5_write_attribute(group_id, names[0], H5T_NATIVE_INT, ds_id, addresses[0]);

  // Close the group
  H5Gclose(group_id);

#ifdef GITREF_STR
  // Save the git ref if requested
  char tempstr[45];
  hid_t attr_id;

  sprintf(tempstr, GITREF_STR);
  attr_id = H5Acreate(file_id, "GitRef", str_t, ds_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_id, str_t, tempstr);
  H5Aclose(attr_id);
#endif

#ifdef USE_TOCF
  if(run_globals->params.TOCF_Flag)
  {
    // Open the group
    group_id = H5Gcreate(file_id, "InputParams/tocf", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    ii=0;
    addresses[ii] = &(tocf_params.dim);
    sprintf(names[ii++],  "dim");
    addresses[ii] = &(tocf_params.HII_dim);
    sprintf(names[ii++],  "HII_dim");
    addresses[ii] = &(tocf_params.numcores);
    sprintf(names[ii++],  "numcores");
    addresses[ii] = &(tocf_params.HII_filter);
    sprintf(names[ii++],  "HII_filter");
    addresses[ii] = &(tocf_params.uvb_feedback);
    sprintf(names[ii++],  "uvb_feedback");
    addresses[ii] = &(tocf_params.compute_mfp);
    sprintf(names[ii++],  "compute_mfp");

    for(int jj=0; jj<ii; jj++)
      h5_write_attribute(group_id, names[jj], H5T_NATIVE_INT, ds_id, addresses[jj]);

    ii=0;
    addresses[ii] = &(tocf_params.ram);
    sprintf(names[ii++],  "ram");
    addresses[ii] = &(tocf_params.HII_eff_factor);
    sprintf(names[ii++],  "HII_eff_factor");
    addresses[ii] = &(tocf_params.r_bubble_min);
    sprintf(names[ii++],  "r_bubble_min");
    addresses[ii] = &(tocf_params.r_bubble_max);
    sprintf(names[ii++],  "r_bubble_max");
    addresses[ii] = &(tocf_params.gamma_halo_bias);
    sprintf(names[ii++],  "gamma_halo_bias");
    addresses[ii] = &(tocf_params.delta_r_HII_factor);
    sprintf(names[ii++],  "delta_r_HII_factor");
    addresses[ii] = &(tocf_params.m_0_sm);
    sprintf(names[ii++],  "m_0_sm");
    addresses[ii] = &(tocf_params.a_sm);
    sprintf(names[ii++],  "a_sm");
    addresses[ii] = &(tocf_params.b_sm);
    sprintf(names[ii++],  "b_sm");
    addresses[ii] = &(tocf_params.c_sm);
    sprintf(names[ii++],  "c_sm");
    addresses[ii] = &(tocf_params.d_sm);
    sprintf(names[ii++],  "d_sm");

    for(int jj=0; jj<ii; jj++)
      h5_write_attribute(group_id, names[jj], H5T_NATIVE_FLOAT, ds_id, addresses[jj]);

    ii=0;
    addresses[ii] = &(tocf_params.ion_tvir_min);
    sprintf(names[ii++],  "ion_tvir_min");

    for(int jj=0; jj<ii; jj++)
      h5_write_attribute(group_id, names[jj], H5T_NATIVE_DOUBLE, ds_id, addresses[jj]);

    // Close the group
    H5Gclose(group_id);
  }
#endif

  // Close the HDF5 file.
  H5Fclose(file_id);

  H5Sclose(ds_id);
  H5Tclose(str_t);

}

static void inline save_walk_indices(
  run_globals_t *run_globals,
  hid_t               file_id,
  int                 i_out,
  int                 prev_i_out,
  int                *descendant_index,
  int                *first_progenitor_index,
  int                *next_progenitor_index,
  int                 old_count,
  int                 n_write)
{

  hsize_t dim[1];
  char target[50];

  dim[0] = (hsize_t)old_count;
  sprintf(target, "Snap%03d/DescendantIndices", (run_globals->ListOutputSnaps)[prev_i_out]);
  H5LTmake_dataset(file_id,target,1,dim,H5T_NATIVE_INT,descendant_index);
  sprintf(target, "Snap%03d/NextProgenitorIndices", (run_globals->ListOutputSnaps)[prev_i_out]);
  H5LTmake_dataset(file_id,target,1,dim,H5T_NATIVE_INT,next_progenitor_index);

  dim[0] = (hsize_t)n_write;
  sprintf(target, "Snap%03d/FirstProgenitorIndices", (run_globals->ListOutputSnaps)[i_out]);
  H5LTmake_dataset(file_id,target,1,dim,H5T_NATIVE_INT,first_progenitor_index);

}


void write_snapshot(run_globals_t *run_globals, int n_write, int i_out, int *last_n_write)
{

  /*
   * Write a batch of galaxies to the output HDF5 table.
   */

  hid_t                 file_id;
  hid_t                 group_id;
  hid_t                 ds_id;
  hsize_t               chunk_size             = 10000;
  galaxy_output_t *output_buffer          = NULL;
  int                  *fill_data              = NULL;
  char                  target_group[20];
  galaxy_t        *gal                    = NULL;
  hdf5_output_t    h5props                = run_globals->hdf5props;
  int                   gal_count              = 0;
  int                   old_count              = 0;
  hsize_t               dims                   = 1;
  double                temp;
  int                  *descendant_index       = NULL;
  int                  *first_progenitor_index = NULL;
  int                  *next_progenitor_index  = NULL;
  int                   calc_descendants_i_out = -1;
  int                   prev_snapshot          = -1;
  int                   index                  = -1;
  int                   corrected_snapshot;

  SID_log("rank %d: Writing output file...", SID_LOG_OPEN|SID_LOG_TIMER|SID_LOG_ALLRANKS, SID.My_rank);

  // Create the file.
  file_id = H5Fopen(run_globals->FNameOut, H5F_ACC_RDWR, H5P_DEFAULT);

  // Create the relevant group.
  sprintf(target_group, "Snap%03d", (run_globals->ListOutputSnaps)[i_out]);
  group_id = H5Gcreate(file_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Make the table
  H5TBmake_table("Galaxies", group_id, "Galaxies",
      h5props.n_props, n_write, h5props.dst_size, h5props.field_names,
      h5props.dst_offsets, h5props.field_types, chunk_size, fill_data, 0,
      NULL );

  // If the immediately preceeding snapshot was also written, then save the
  // descendent indices
  prev_snapshot = run_globals->ListOutputSnaps[i_out]-1;
  if (i_out > 0)
    for (int ii=0; ii<NOUT; ii++)
      if (run_globals->ListOutputSnaps[ii] == prev_snapshot)
      {
        calc_descendants_i_out = ii;
        break;
      }

  // Assign the write order indices to each galaxy and store the old indices if required
  gal_count = 0;
  old_count = 0;
  if (calc_descendants_i_out>-1)
  {
    descendant_index       = SID_malloc(sizeof(int)* (*last_n_write));
    next_progenitor_index  = SID_malloc(sizeof(int)* (*last_n_write));
    first_progenitor_index = SID_malloc(sizeof(int)* n_write);

    for (int ii=0; ii<*last_n_write; ii++)
    {
      descendant_index[ii] = -1;
      next_progenitor_index[ii] = -1;
    }
    for (int ii=0; ii<n_write; ii++)
      first_progenitor_index[ii] = -1;

    gal = run_globals->FirstGal;
    while (gal!=NULL) {
      if (gal->Type < 3)
      {
        if (gal->output_index > -1)
        {
          first_progenitor_index[gal_count] = gal->output_index;
          descendant_index[gal->output_index] = gal_count;
          old_count++;
        }
        gal->output_index = gal_count++;
      }
      gal=gal->Next;
    }

    // Here we want to walk the progenitor indices to tag on galaxies which
    // have merged in this timestep and also set their descendant_index.
    gal = run_globals->FirstGal;
    while (gal!=NULL) {
      if (gal->Type == 3)
      {
        descendant_index[gal->output_index] = gal->MergerTarget->output_index;
        old_count++;
        index = first_progenitor_index[gal->MergerTarget->output_index];
        if (index>-1)
        {
          while(next_progenitor_index[index]>-1)
            index = next_progenitor_index[index];
          next_progenitor_index[index] = gal->output_index;
        }
      }
      gal = gal->Next;
    }

    save_walk_indices(run_globals, file_id, i_out, calc_descendants_i_out,
        descendant_index, first_progenitor_index, next_progenitor_index,
        *last_n_write, n_write);

    // Free the allocated arrays
    SID_free(SID_FARG first_progenitor_index);
    SID_free(SID_FARG next_progenitor_index);
    SID_free(SID_FARG descendant_index);

  } else
  {
    gal = run_globals->FirstGal;
    while (gal!=NULL) {
      if (gal->Type < 3)
        gal->output_index = gal_count++;
      gal = gal->Next;
    }
  }

  if (n_write != gal_count)
  {
    SID_log("We don't have the expected number of galaxies in save...", SID_LOG_COMMENT);
    SID_log("gal_count=%d, n_write=%d", SID_LOG_COMMENT, gal_count, n_write);
    ABORT(EXIT_FAILURE);
  }

  // Write the galaxies.
  // In order to speed things up, we will chunk our write.
  // This can cause significant memory overhead if `chunk_size` is large.
  gal_count = 0;
  gal = run_globals->FirstGal;
  output_buffer = SID_malloc(sizeof(galaxy_output_t)*chunk_size);
  int buffer_count = 0;
  while (gal!=NULL) {
    // Don't output galaxies which merged at this timestep
    if (gal->Type < 3)
    {
      prepare_galaxy_for_output(run_globals, *gal, &(output_buffer[buffer_count]), i_out);
      buffer_count++;
    }
    if(buffer_count==chunk_size)
    {
      H5TBwrite_records(group_id, "Galaxies", gal_count, buffer_count, h5props.dst_size,
          h5props.dst_offsets, h5props.dst_field_sizes, output_buffer);
      gal_count += buffer_count;
      buffer_count = 0;
    }
    gal = gal->Next;
  }

  // Write any remaining galaxies in the buffer
  if(buffer_count>0)
  {
    H5TBwrite_records(group_id, "Galaxies", gal_count, buffer_count, h5props.dst_size,
        h5props.dst_offsets, h5props.dst_field_sizes, output_buffer);
    gal_count += buffer_count;
  }

  if (n_write != gal_count)
  {
    SID_log("We don't have the expected number of galaxies in save...", SID_LOG_COMMENT);
    SID_log("gal_count=%d, n_write=%d", SID_LOG_COMMENT, gal_count, n_write);
    ABORT(EXIT_FAILURE);
  }

  // Free the output buffer
  SID_free(SID_FARG output_buffer);

  // Save a few useful attributes
  ds_id = H5Screate_simple(1, &dims, NULL);

  corrected_snapshot = get_corrected_snapshot(run_globals, run_globals->ListOutputSnaps[i_out]);
  h5_write_attribute(group_id, "Redshift", H5T_NATIVE_DOUBLE, ds_id, &(run_globals->ZZ[run_globals->ListOutputSnaps[i_out]]));
  h5_write_attribute(group_id, "CorrectedSnap", H5T_NATIVE_INT, ds_id, &corrected_snapshot);

  temp = run_globals->LTTime[run_globals->ListOutputSnaps[i_out]] * run_globals->units.UnitLength_in_cm / run_globals->units.UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR / run_globals->params.Hubble_h;
  h5_write_attribute(group_id, "LTTime", H5T_NATIVE_DOUBLE, ds_id, &temp);

  temp = global_ionizing_emmisivity(run_globals);
  temp *= pow(run_globals->params.Hubble_h, 3);  // Factor out hubble constants
  h5_write_attribute(group_id, "GlobalIonizingEmissivity", H5T_NATIVE_DOUBLE, ds_id, &temp);

  H5Sclose(ds_id);

#ifdef USE_TOCF
  if(run_globals->params.TOCF_Flag)
    save_tocf_grids(run_globals, group_id, run_globals->ListOutputSnaps[i_out]);
#endif

  // Close the group.
  H5Gclose(group_id);

  // Close the file.
  H5Fclose(file_id);

  // Update the value of last_n_write
  *last_n_write = n_write;

  SID_log("rank %d: ...done", SID_LOG_CLOSE|SID_LOG_ALLRANKS, SID.My_rank);

}
