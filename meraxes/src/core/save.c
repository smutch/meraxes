#include <unistd.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>
#include "meraxes.h"
#include "parse_paramfile.h"
#include "git.h"

float current_mwmsa(galaxy_t *gal, int i_snap)
{
  double *LTTime = run_globals.LTTime;
  double mwmsa_num = gal->mwmsa_num;
  double mwmsa_denom = gal->mwmsa_denom;
  int snapshot = run_globals.ListOutputSnaps[i_snap];

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
  {
    mwmsa_num += gal->NewStars[ii] * LTTime[snapshot-ii];
    mwmsa_denom += gal->NewStars[ii];
  }

  return (float)((mwmsa_num/mwmsa_denom) - LTTime[snapshot]);
}


void prepare_galaxy_for_output(
  galaxy_t         gal,
  galaxy_output_t *galout,
  int              i_snap)
{
  run_units_t *units = &(run_globals.units);

  galout->id_MBP = (long long)gal.id_MBP;
  galout->ID     = (int)gal.ID;
  galout->Type   = (int)gal.Type;
  if (!gal.ghost_flag)
  {
    galout->CentralGal = (int)(gal.Halo->FOFGroup->FirstOccupiedHalo->Galaxy->output_index);
    galout->FOFMvir    = (float)(gal.Halo->FOFGroup->Mvir);
  }
  else
  {
    galout->CentralGal = -1;
    galout->FOFMvir    = -1.0;
  }
  galout->GhostFlag = (int)gal.ghost_flag;

  for (int ii = 0; ii < 3; ii++)
  {
    galout->Pos[ii] = (float)(gal.Pos[ii]);
    galout->Vel[ii] = (float)(gal.Vel[ii]);
  }

  galout->Len                = (int)(gal.Len);
  galout->MaxLen             = (int)(gal.MaxLen);
  galout->Mvir               = (float)(gal.Mvir);
  galout->Rvir               = (float)(gal.Rvir);
  galout->Vvir               = (float)(gal.Vvir);
  galout->Vmax               = (float)(gal.Vmax);
  galout->Spin               = (float)(gal.Spin);
  galout->HotGas             = (float)(gal.HotGas);
  galout->MetalsHotGas       = (float)(gal.MetalsHotGas);
  galout->ColdGas            = (float)(gal.ColdGas);
  galout->MetalsColdGas      = (float)(gal.MetalsColdGas);
  galout->Mcool              = (float)(gal.Mcool);
  galout->StellarMass        = (float)(gal.StellarMass);
  galout->GrossStellarMass   = (float)(gal.GrossStellarMass);
  galout->BlackHoleMass      = (float)(gal.BlackHoleMass);
  galout->DiskScaleLength    = (float)(gal.DiskScaleLength);
  galout->MetalsStellarMass  = (float)(gal.MetalsStellarMass);
  galout->Sfr                = (float)(gal.Sfr * units->UnitMass_in_g / units->UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  galout->EjectedGas         = (float)(gal.EjectedGas);
  galout->MetalsEjectedGas   = (float)(gal.MetalsEjectedGas);
  galout->Rcool              = (float)(gal.Rcool);
  galout->Cos_Inc            = (float)(gal.Cos_Inc);
  galout->BaryonFracModifier = (float)(gal.BaryonFracModifier);
  galout->MvirCrit           = (float)(gal.MvirCrit);
  galout->MergTime           = (float)(gal.MergTime * units->UnitLength_in_cm / units->UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR);
  galout->MergerStartRadius  = (float)(gal.MergerStartRadius);
  galout->MWMSA              = current_mwmsa(&gal, i_snap);

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
    galout->NewStars[ii] = (float)(gal.NewStars[ii]);

  prepare_magnitudes_for_output(gal, galout, i_snap);
}

void calc_hdf5_props()
{
  /*
   * Prepare an HDF5 to receive the output galaxy data.
   * Here we store the data in an hdf5 table for easily appending new data.
   */

  if (!run_globals.params.FlagMCMC)
  {
    hdf5_output_t *h5props = &(run_globals.hdf5props);
    galaxy_output_t galout;
    int i;                                                // dummy

    h5props->n_props = 39;

#ifdef CALC_MAGS
    // If we are calculating any magnitudes then increment the number of
    // output properties appropriately.
    int n_photo_bands = run_globals.photo->NBands;
    h5props->n_props         += 2;
    h5props->array_nmag_f_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){ n_photo_bands });
#endif

    // Size of a single galaxy entry.
    h5props->dst_size = sizeof(galaxy_output_t);

    // Create datatypes for different size arrays
    h5props->array3f_tid       = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){ 3 });
    h5props->array_nhist_f_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, (hsize_t[]){ N_HISTORY_SNAPS });

    // Calculate the offsets of our struct members in memory
    h5props->dst_offsets = SID_malloc(sizeof(size_t) * h5props->n_props);
    // Calculate the sizes of our struct members in memory.
    h5props->dst_field_sizes = SID_malloc(sizeof(size_t) * h5props->n_props);
    // Give each galaxy property a field name in the table
    h5props->field_names = SID_malloc(sizeof(const char*) * h5props->n_props);
    // Assign a type to each galaxy property field in the table.
    h5props->field_types = SID_malloc(sizeof(hid_t) * h5props->n_props);
    // Store the **output** units of each property for writing to the master file.
    // Units should be compatible with the python astropy.units module.
    h5props->field_units = SID_malloc(sizeof(const char*) * h5props->n_props);
    // Store the **output** h conversion for each property.  The string will be
    // parsed by python eval(), substituting h for the appropriate value at read
    // time and v for the property value.
    h5props->field_h_conv = SID_malloc(sizeof(const char*) * h5props->n_props);

    i = 0;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, id_MBP);
    h5props->dst_field_sizes[i] = sizeof(galout.id_MBP);
    h5props->field_names[i]     = "id_MBP";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_LLONG;

#ifdef CALC_MAGS
    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Mag);
    h5props->dst_field_sizes[i] = sizeof(float) * n_photo_bands;
    h5props->field_names[i]     = "Mag";
    h5props->field_units[i]     = "mag";
    h5props->field_h_conv[i]    = "v - 5*log10(h)";
    h5props->field_types[i++]   = h5props->array_nmag_f_tid;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MagDust);
    h5props->dst_field_sizes[i] = sizeof(float) * n_photo_bands;
    h5props->field_names[i]     = "MagDust";
    h5props->field_units[i]     = "mag";
    h5props->field_h_conv[i]     = "v - 5*log10(h)";
    h5props->field_types[i++]   = h5props->array_nmag_f_tid;
#endif

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, ID);
    h5props->dst_field_sizes[i] = sizeof(galout.ID);
    h5props->field_names[i]     = "ID";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_LLONG;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Type);
    h5props->dst_field_sizes[i] = sizeof(galout.Type);
    h5props->field_names[i]     = "Type";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_INT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, CentralGal);
    h5props->dst_field_sizes[i] = sizeof(galout.CentralGal);
    h5props->field_names[i]     = "CentralGal";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_INT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, GhostFlag);
    h5props->dst_field_sizes[i] = sizeof(galout.GhostFlag);
    h5props->field_names[i]     = "GhostFlag";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_INT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Len);
    h5props->dst_field_sizes[i] = sizeof(galout.Len);
    h5props->field_names[i]     = "Len";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_INT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MaxLen);
    h5props->dst_field_sizes[i] = sizeof(galout.MaxLen);
    h5props->field_names[i]     = "MaxLen";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_INT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Pos);
    h5props->dst_field_sizes[i] = sizeof(galout.Pos);
    h5props->field_names[i]     = "Pos";
    h5props->field_units[i]     = "Mpc";  // comoving
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = h5props->array3f_tid;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Vel);
    h5props->dst_field_sizes[i] = sizeof(galout.Vel);
    h5props->field_names[i]     = "Vel";
    h5props->field_units[i]     = "km/s";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = h5props->array3f_tid;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Spin);
    h5props->dst_field_sizes[i] = sizeof(galout.Spin);
    h5props->field_names[i]     = "Spin";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Mvir);
    h5props->dst_field_sizes[i] = sizeof(galout.Mvir);
    h5props->field_names[i]     = "Mvir";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Rvir);
    h5props->dst_field_sizes[i] = sizeof(galout.Rvir);
    h5props->field_names[i]     = "Rvir";
    h5props->field_units[i]     = "Mpc";  // physical
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Vvir);
    h5props->dst_field_sizes[i] = sizeof(galout.Vvir);
    h5props->field_names[i]     = "Vvir";
    h5props->field_units[i]     = "km/s";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Vmax);
    h5props->dst_field_sizes[i] = sizeof(galout.Vmax);
    h5props->field_names[i]     = "Vmax";
    h5props->field_units[i]     = "km/s";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, FOFMvir);
    h5props->dst_field_sizes[i] = sizeof(galout.FOFMvir);
    h5props->field_names[i]     = "FOFMvir";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, HotGas);
    h5props->dst_field_sizes[i] = sizeof(galout.HotGas);
    h5props->field_names[i]     = "HotGas";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MetalsHotGas);
    h5props->dst_field_sizes[i] = sizeof(galout.MetalsHotGas);
    h5props->field_names[i]     = "MetalsHotGas";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, ColdGas);
    h5props->dst_field_sizes[i] = sizeof(galout.ColdGas);
    h5props->field_names[i]     = "ColdGas";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MetalsColdGas);
    h5props->dst_field_sizes[i] = sizeof(galout.MetalsColdGas);
    h5props->field_names[i]     = "MetalsColdGas";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Mcool);
    h5props->dst_field_sizes[i] = sizeof(galout.Mcool);
    h5props->field_names[i]     = "Mcool";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, DiskScaleLength);
    h5props->dst_field_sizes[i] = sizeof(galout.DiskScaleLength);
    h5props->field_names[i]     = "DiskScaleLength";
    h5props->field_units[i]     = "Mpc";  // real
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, StellarMass);
    h5props->dst_field_sizes[i] = sizeof(galout.StellarMass);
    h5props->field_names[i]     = "StellarMass";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, GrossStellarMass);
    h5props->dst_field_sizes[i] = sizeof(galout.GrossStellarMass);
    h5props->field_names[i]     = "GrossStellarMass";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MetalsStellarMass);
    h5props->dst_field_sizes[i] = sizeof(galout.MetalsStellarMass);
    h5props->field_names[i]     = "MetalsStellarMass";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Sfr);
    h5props->dst_field_sizes[i] = sizeof(galout.Sfr);
    h5props->field_names[i]     = "Sfr";
    h5props->field_units[i]     = "solMass/yr";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, EjectedGas);
    h5props->dst_field_sizes[i] = sizeof(galout.EjectedGas);
    h5props->field_names[i]     = "EjectedGas";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MetalsEjectedGas);
    h5props->dst_field_sizes[i] = sizeof(galout.MetalsEjectedGas);
    h5props->field_names[i]     = "MetalsEjectedGas";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, BlackHoleMass);
    h5props->dst_field_sizes[i] = sizeof(galout.BlackHoleMass);
    h5props->field_names[i]     = "BlackHoleMass";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Rcool);
    h5props->dst_field_sizes[i] = sizeof(galout.Rcool);
    h5props->field_names[i]     = "Rcool";
    h5props->field_units[i]     = "Mpc";  // real
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Cos_Inc);
    h5props->dst_field_sizes[i] = sizeof(galout.Cos_Inc);
    h5props->field_names[i]     = "Cos_Inc";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MergTime);
    h5props->dst_field_sizes[i] = sizeof(galout.MergTime);
    h5props->field_names[i]     = "MergTime";
    h5props->field_units[i]     = "Myr";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MergerStartRadius);
    h5props->dst_field_sizes[i] = sizeof(galout.MergerStartRadius);
    h5props->field_names[i]     = "MergerStartRadius";
    h5props->field_units[i]     = "Mpc";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, BaryonFracModifier);
    h5props->dst_field_sizes[i] = sizeof(galout.BaryonFracModifier);
    h5props->field_names[i]     = "BaryonFracModifier";
    h5props->field_units[i]     = "None";
    h5props->field_h_conv[i]    = "None";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MvirCrit);
    h5props->dst_field_sizes[i] = sizeof(galout.MvirCrit);
    h5props->field_names[i]     = "MvirCrit";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, MWMSA);
    h5props->dst_field_sizes[i] = sizeof(galout.MWMSA);
    h5props->field_names[i]     = "MWMSA";
    h5props->field_units[i]     = "Myr";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = H5T_NATIVE_FLOAT;

    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, NewStars);
    h5props->dst_field_sizes[i] = sizeof(galout.NewStars);
    h5props->field_names[i]     = "NewStars";
    h5props->field_units[i]     = "1e10 solMass";
    h5props->field_h_conv[i]    = "v/h";
    h5props->field_types[i++]   = h5props->array_nhist_f_tid;

    // DEBUG
    if (i != h5props->n_props)
    {
      SID_log_error("Incorrect number of galaxy properties in HDF5 file.");
      ABORT(EXIT_FAILURE);
    }
  }
}

void prep_hdf5_file()
{
  hid_t file_id;

  // create a new file
  if (access(run_globals.FNameOut, F_OK) != -1)
    remove(run_globals.FNameOut);
  file_id = H5Fcreate(run_globals.FNameOut, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // store the file number and total number of cores
  H5LTset_attribute_int(file_id, "/", "iCore", &(SID.My_rank), 1);
  H5LTset_attribute_int(file_id, "/", "NCores", &(SID.n_proc), 1);

  // close the file
  H5Fclose(file_id);
}


void create_master_file()
{
  hid_t file_id, group_id;
  char fname[STRLEN];
  hdf5_output_t *h5props = &(run_globals.hdf5props);
  char **params_tag = h5props->params_tag;
  void **params_addr = h5props->params_addr;
  int   *params_type = h5props->params_type;
  int    params_count = h5props->params_count;

  SID_log("Creating master file...", SID_LOG_OPEN | SID_LOG_TIMER);

  // Create a new file
  sprintf(fname, "%s/%s.hdf5", run_globals.params.OutputDir, run_globals.params.FileNameGalaxies);
  if (access(fname, F_OK) != -1)
    remove(fname);
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Open the group
  {
    const char *group_name = {"InputParams"};
    group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Save all of the input params
    for (int ii = 0; (ii < params_count) && (params_type[ii] != PARAM_TYPE_UNUSED); ii++)
    {
      switch (params_type[ii])
      {
        case PARAM_TYPE_STRING:
          H5LTset_attribute_string(file_id, group_name, params_tag[ii], params_addr[ii]);
          break;
        case PARAM_TYPE_INT:
          H5LTset_attribute_int(file_id, group_name, params_tag[ii], params_addr[ii], 1);
          break;
        case PARAM_TYPE_DOUBLE:
          H5LTset_attribute_double(file_id, group_name, params_tag[ii], params_addr[ii], 1);
          break;
        case PARAM_TYPE_FLOAT:
          H5LTset_attribute_float(file_id, group_name, params_tag[ii], params_addr[ii], 1);
          break;
        case PARAM_TYPE_LONGLONG:
          H5LTset_attribute_long_long(file_id, group_name, params_tag[ii], params_addr[ii], 1);
          break;
        default:
          ABORT(EXIT_FAILURE);
          break;
      }
    }

    // Close the group
    H5Gclose(group_id);
  }

  // save the units of each galaxy property and grid
  {
    const char *group_name[2] = {"Units", "HubbleConversions"};
    for(int ii=0; ii<2; ii++)
    {
      group_id = H5Gcreate(file_id, group_name[ii], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group_id);
    }

    for (int ii = 0; ii < h5props->n_props; ii++)
    {
      H5LTset_attribute_string(file_id, group_name[0], h5props->field_names[ii], h5props->field_units[ii]);
      H5LTset_attribute_string(file_id, group_name[1], h5props->field_names[ii], h5props->field_h_conv[ii]);
    }
  }

  if (run_globals.params.Flag_PatchyReion)
  {
    {
      const char *group_name = {"Units/Grids"};
      group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5LTset_attribute_string(file_id, group_name, "xH", "None");
      H5LTset_attribute_string(file_id, group_name, "J_21", "10e-21 erg/s/Hz/cm/cm/sr");
      H5LTset_attribute_string(file_id, group_name, "J_21_at_ionization", "10e-21 erg/s/Hz/cm/cm/sr");
      H5LTset_attribute_string(file_id, group_name, "z_at_ionization", "None");
      H5LTset_attribute_string(file_id, group_name, "Mvir_crit", "1e10 solMass");
      H5LTset_attribute_string(file_id, group_name, "StellarMass", "1e10 solMass");
      H5LTset_attribute_string(file_id, group_name, "Sfr", "solMass/yr");
      H5LTset_attribute_string(file_id, group_name, "deltax", "None");
      H5Gclose(group_id);
    }

    {
      const char *group_name = {"HubbleConversions/Grids"};
      group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5LTset_attribute_string(file_id, group_name, "xH", "None");
      H5LTset_attribute_string(file_id, group_name, "J_21", "v*(h**2)");
      H5LTset_attribute_string(file_id, group_name, "J_21_at_ionization", "v*(h**2)");
      H5LTset_attribute_string(file_id, group_name, "z_at_ionization", "None");
      H5LTset_attribute_string(file_id, group_name, "Mvir_crit", "v/h");
      H5LTset_attribute_string(file_id, group_name, "StellarMass", "v/h");
      H5LTset_attribute_string(file_id, group_name, "Sfr", "None");
      H5LTset_attribute_string(file_id, group_name, "deltax", "None");
      H5Gclose(group_id);
    }
  }


#ifdef GITREF_STR
  // Save the git ref and diff if requested
  H5LTmake_dataset_string(file_id, "gitdiff", GITDIFF_STR);
  H5LTset_attribute_string(file_id, "gitdiff", "gitref", GITREF_STR);
#endif

  // save the number of cores used in this run
  H5LTset_attribute_int(file_id, "/", "NCores", &(SID.n_proc), 1);


  char target_group[50];
  char source_ds[50];
  char source_group[50];
  char source_file[STRLEN];
  char relative_source_file[50];
  hid_t snap_group_id;
  hid_t source_file_id;
  hid_t source_group_id;
  hsize_t core_n_gals;
  double temp;
  int unsampled_snapshot;

  // Now create soft links to all of the files and datasets that make up this run
  for (int i_out = 0, snap_n_gals = 0; i_out < run_globals.NOutputSnaps; i_out++, snap_n_gals = 0, temp = 0)
  {
    sprintf(target_group, "Snap%03d", run_globals.ListOutputSnaps[i_out]);
    snap_group_id = H5Gcreate(file_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (int i_core = 0; i_core < SID.n_proc; i_core++)
    {
      sprintf(target_group, "Core%d", i_core);
      group_id = H5Gcreate(snap_group_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      sprintf(source_file, "%s/%s_%d.hdf5", run_globals.params.OutputDir, run_globals.params.FileNameGalaxies, i_core);
      sprintf(relative_source_file, "%s_%d.hdf5", run_globals.params.FileNameGalaxies, i_core);
      sprintf(source_ds, "Snap%03d/Galaxies", run_globals.ListOutputSnaps[i_out]);
      H5Lcreate_external(relative_source_file, source_ds, group_id, "Galaxies", H5P_DEFAULT, H5P_DEFAULT);

      source_file_id = H5Fopen(source_file, H5F_ACC_RDONLY, H5P_DEFAULT);
      H5TBget_table_info(source_file_id, source_ds, NULL, &core_n_gals);
      snap_n_gals                += (int)core_n_gals;

      if (i_core == 0)
        H5LTget_attribute_int(source_file_id, source_ds, "unsampled_snapshot", &unsampled_snapshot);

      // if they exists, then also create a link to walk indices
      sprintf(source_group, "Snap%03d", run_globals.ListOutputSnaps[i_out]);
      source_group_id = H5Gopen(source_file_id, source_group, H5P_DEFAULT);
      if (H5LTfind_dataset(source_group_id, "FirstProgenitorIndices"))
      {
        sprintf(source_ds, "Snap%03d/FirstProgenitorIndices", run_globals.ListOutputSnaps[i_out]);
        H5Lcreate_external(relative_source_file, source_ds, group_id, "FirstProgenitorIndices", H5P_DEFAULT, H5P_DEFAULT);
      }
      if (H5LTfind_dataset(source_group_id, "NextProgenitorIndices"))
      {
        sprintf(source_ds, "Snap%03d/NextProgenitorIndices", run_globals.ListOutputSnaps[i_out]);
        H5Lcreate_external(relative_source_file, source_ds, group_id, "NextProgenitorIndices", H5P_DEFAULT, H5P_DEFAULT);
      }
      if (H5LTfind_dataset(source_group_id, "DescendantIndices"))
      {
        sprintf(source_ds, "Snap%03d/DescendantIndices", run_globals.ListOutputSnaps[i_out]);
        H5Lcreate_external(relative_source_file, source_ds, group_id, "DescendantIndices", H5P_DEFAULT, H5P_DEFAULT);
      }

      H5Gclose(source_group_id);
      H5Gclose(group_id);
      H5Fclose(source_file_id);

      if((i_core == 0) && (run_globals.params.Flag_PatchyReion))
      {
        // create links to the 21cmFAST grids that exist
        gen_grids_fname(run_globals.ListOutputSnaps[i_out], relative_source_file, true);
        gen_grids_fname(run_globals.ListOutputSnaps[i_out], source_file, false);
        if (access(source_file, F_OK) != -1)
        {
          source_file_id = H5Fopen(source_file, H5F_ACC_RDONLY, H5P_DEFAULT);
          H5Lcreate_external(relative_source_file, "/", snap_group_id, "Grids", H5P_DEFAULT, H5P_DEFAULT);
          H5Fclose(source_file_id);
        }

        // sprintf(source_ds, "Snap%03d/PowerSpectrum", run_globals.ListOutputSnaps[i_out]);
        // if ((H5LTpath_valid(source_file_id, source_ds, FALSE)))
        // {
        //   sprintf(target_ds, "PowerSpectrum");
        //   H5Lcreate_external(relative_source_file, source_ds, snap_group_id, target_ds, H5P_DEFAULT, H5P_DEFAULT);
        // }

        // sprintf(source_ds, "Snap%03d/RegionSizeDist", run_globals.ListOutputSnaps[i_out]);
        // if ((H5LTpath_valid(source_file_id, source_ds, FALSE)))
        // {
        //   sprintf(target_ds, "RegionSizeDist");
        //   H5Lcreate_external(relative_source_file, source_ds, snap_group_id, target_ds, H5P_DEFAULT, H5P_DEFAULT);
        // }
      }

    }

    // Save a few useful attributes
    sprintf(target_group, "Snap%03d", run_globals.ListOutputSnaps[i_out]);

    // save the total number of galaxies at this snapshot
    H5LTset_attribute_int(file_id, target_group, "NGalaxies", &snap_n_gals, 1);

    H5LTset_attribute_double(file_id, target_group, "Redshift", &(run_globals.ZZ[run_globals.ListOutputSnaps[i_out]]), 1);
    H5LTset_attribute_int(file_id, target_group, "UnsampledSnapshot", &unsampled_snapshot, 1);

    temp = run_globals.LTTime[run_globals.ListOutputSnaps[i_out]] * run_globals.units.UnitLength_in_cm / run_globals.units.UnitVelocity_in_cm_per_s / SEC_PER_MEGAYEAR;
    H5LTset_attribute_double(file_id, target_group, "LTTime", &temp, 1);

    H5Gclose(snap_group_id);
  }

  // Close the HDF5 file.
  H5Fclose(file_id);

  SID_log(" ...done", SID_LOG_CLOSE);
}

static void inline save_walk_indices(
  hid_t          file_id,
  int            i_out,
  int            prev_i_out,
  int           *descendant_index,
  int           *first_progenitor_index,
  int           *next_progenitor_index,
  int            old_count,
  int            n_write)
{
  hsize_t dim[1];
  char target[50];

  dim[0] = (hsize_t)old_count;
  sprintf(target, "Snap%03d/DescendantIndices", (run_globals.ListOutputSnaps)[prev_i_out]);
  H5LTmake_dataset(file_id, target, 1, dim, H5T_NATIVE_INT, descendant_index);
  sprintf(target, "Snap%03d/NextProgenitorIndices", (run_globals.ListOutputSnaps)[prev_i_out]);
  H5LTmake_dataset(file_id, target, 1, dim, H5T_NATIVE_INT, next_progenitor_index);

  dim[0] = (hsize_t)n_write;
  sprintf(target, "Snap%03d/FirstProgenitorIndices", (run_globals.ListOutputSnaps)[i_out]);
  H5LTmake_dataset(file_id, target, 1, dim, H5T_NATIVE_INT, first_progenitor_index);
}


static inline bool pass_write_check(galaxy_t *gal, bool flag_merger)
{
  if (
    // Test for non-merger galaxy to be written in the current snap
    (!flag_merger && (gal->Type < 3)) // && ((gal->output_index > -1) || (gal->StellarMass >= 1e-10)))
    // and this is the test for a merger to be accounted for in descendant / progenitor arrays
    || (flag_merger && (gal->Type == 3) && (gal->output_index > -1))
    )
    return true;
  else
    return false;
}


void write_snapshot(
  int            n_write,
  int            i_out,
  int           *last_n_write,
  trees_info_t  *trees_info)
{
  /*
   * Write a batch of galaxies to the output HDF5 table.
   */

  hid_t file_id;
  hid_t group_id;
  hsize_t chunk_size             = 10000;
  galaxy_output_t *output_buffer = NULL;
  int *fill_data                 = NULL;
  char target_group[20];
  galaxy_t *gal               = NULL;
  hdf5_output_t h5props       = run_globals.hdf5props;
  int gal_count               = 0;
  int old_count               = 0;
  int *first_progenitor_index = NULL;
  int *next_progenitor_index  = NULL;
  int calc_descendants_i_out  = -1;
  int prev_snapshot           = -1;
  int write_count             = 0;

  SID_log("Writing output file (n_write = %d)...", SID_LOG_OPEN | SID_LOG_TIMER, n_write);

  // We aren't going to write any galaxies that have zero stellar mass, so
  // modify n_write appropriately...
  gal = run_globals.FirstGal;
  while (gal != NULL)
  {
    if (pass_write_check(gal, false))
      write_count++;
    gal = gal->Next;
  }

  if (n_write != write_count)
  {
    SID_log("Excluding %d ~zero mass galaxies...", SID_LOG_COMMENT, n_write - write_count);
    SID_log("New write count = %d", SID_LOG_COMMENT, write_count);
    n_write = write_count;
  }

  // Create the file.
  file_id = H5Fopen(run_globals.FNameOut, H5F_ACC_RDWR, H5P_DEFAULT);

  // Create the relevant group.
  sprintf(target_group, "Snap%03d", (run_globals.ListOutputSnaps)[i_out]);
  group_id = H5Gcreate(file_id, target_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Make the table
  H5TBmake_table("Galaxies", group_id, "Galaxies",
                 h5props.n_props, n_write, h5props.dst_size, h5props.field_names,
                 h5props.dst_offsets, h5props.field_types, chunk_size, fill_data, 0,
                 NULL);

  // Store the unsampled snapshot number
  H5LTset_attribute_int(group_id, "Galaxies", "unsampled_snapshot", &(trees_info->unsampled_snapshot), 1);

  // If the immediately preceeding snapshot was also written, then save the
  // descendent indices
  prev_snapshot = run_globals.ListOutputSnaps[i_out] - 1;
  if (i_out > 0)
    for (int ii = 0; ii < run_globals.NOutputSnaps; ii++)
      if (run_globals.ListOutputSnaps[ii] == prev_snapshot)
      {
        calc_descendants_i_out = ii;
        break;
      }

  // Assign the write order indices to each galaxy and store the old indices if required
  gal_count = 0;
  old_count = 0;
  if (calc_descendants_i_out > -1)
  {
    // malloc the arrays
    int *descendant_index       = SID_malloc(sizeof(int) * (*last_n_write));
    next_progenitor_index  = SID_malloc(sizeof(int) * (*last_n_write));
    first_progenitor_index = SID_malloc(sizeof(int) * n_write);

    // initialise all entries to -1
    for (int ii = 0; ii < *last_n_write; ii++)
    {
      descendant_index[ii]      = -1;
      next_progenitor_index[ii] = -1;
    }
    for (int ii = 0; ii < n_write; ii++)
      first_progenitor_index[ii] = -1;

    // loop through the current galaxies and save their first progenitor
    // indices as their previous output_index, and the descendent indices of
    // the last snapshot to what will be the output index when the current
    // galaxy is written.
    gal = run_globals.FirstGal;
    while (gal != NULL)
    {
      if (pass_write_check(gal, false))
      {
        if (gal->output_index > -1)
        {
          first_progenitor_index[gal_count] = gal->output_index;

          assert(gal->output_index < *last_n_write);

          descendant_index[gal->output_index] = gal_count;
          old_count++;
        }
        gal->output_index = gal_count++;
      }
      gal = gal->Next;
    }

    // Here we want to walk the progenitor indices to tag on galaxies which
    // have merged in this timestep and also set their descendant_index.
    gal = run_globals.FirstGal;
    while (gal != NULL)
    {
      if (pass_write_check(gal, true))
      {
        assert((gal->output_index < *last_n_write) && (gal->output_index >= 0));
        descendant_index[gal->output_index] = gal->MergerTarget->output_index;
        old_count++;

        assert(gal->MergerTarget->output_index < n_write);
        if (gal->MergerTarget->output_index >= 0)
        {
          int index = first_progenitor_index[gal->MergerTarget->output_index];
          if (index > -1)
          {
            while (next_progenitor_index[index] > -1)
              index = next_progenitor_index[index];
            next_progenitor_index[index] = gal->output_index;
          }
        }
      }
      gal = gal->Next;
    }

    save_walk_indices(file_id, i_out, calc_descendants_i_out,
                      descendant_index, first_progenitor_index, next_progenitor_index,
                      *last_n_write, n_write);

    // Free the allocated arrays
    SID_free(SID_FARG first_progenitor_index);
    SID_free(SID_FARG next_progenitor_index);
    SID_free(SID_FARG descendant_index);
  }
  else
  {
    gal = run_globals.FirstGal;
    while (gal != NULL)
    {
      if (pass_write_check(gal, false))
        gal->output_index = gal_count++;
      gal = gal->Next;
    }
  }

  if (n_write != gal_count)
  {
    fprintf(stderr, "We don't have the expected number of galaxies in save...");
    fprintf(stderr, "gal_count=%d, n_write=%d", gal_count, n_write);
    ABORT(EXIT_FAILURE);
  }

  // Write the galaxies.
  // In order to speed things up, we will chunk our write.
  // This can cause significant memory overhead if `chunk_size` is large.
  gal_count     = 0;
  gal           = run_globals.FirstGal;
  output_buffer = SID_calloc(sizeof(galaxy_output_t) * (int)chunk_size);
  int buffer_count = 0;
  while (gal != NULL)
  {
    // Don't output galaxies which merged at this timestep
    if (pass_write_check(gal, false))
    {
      prepare_galaxy_for_output(*gal, &(output_buffer[buffer_count]), i_out);
      buffer_count++;
    }
    if (buffer_count == (int)chunk_size)
    {
      H5TBwrite_records(group_id, "Galaxies", gal_count, buffer_count, h5props.dst_size,
                        h5props.dst_offsets, h5props.dst_field_sizes, output_buffer);
      gal_count   += buffer_count;
      buffer_count = 0;
    }
    gal = gal->Next;
  }

  // Write any remaining galaxies in the buffer
  if (buffer_count > 0)
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

  if (run_globals.params.Flag_PatchyReion && !check_if_reionization_complete())
    save_reion_output_grids(run_globals.ListOutputSnaps[i_out]);

  // Close the group.
  H5Gclose(group_id);

  // Close the file.
  H5Fclose(file_id);

  // Update the value of last_n_write
  *last_n_write = n_write;

  SID_log("...done", SID_LOG_CLOSE);
}

