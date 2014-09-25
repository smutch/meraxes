#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "meraxes.h"

void init_luminosities(run_globals_t *run_globals, galaxy_t *gal)
{
#ifdef CALC_MAGS
  for (int ii = 0; ii < NOUT; ii++)
    for (int jj = 0; jj < run_globals->photo.NBands; jj++)
      gal->Lum[jj][ii] = 0.0;
#else
  return;
#endif
}

void sum_luminosities(run_globals_t *run_globals, galaxy_t *parent, galaxy_t *gal, int outputbin)
{
#ifdef CALC_MAGS
  int n_bands = run_globals->photo.NBands;
  for (int ii = 0; ii < n_bands; ii++)
    parent->Lum[ii][outputbin] += gal->Lum[ii][outputbin];
#else
  return;
#endif
}

void prepare_magnitudes_for_output(run_globals_t *run_globals, galaxy_t gal, galaxy_output_t *galout, int i_snap)
{
#ifdef CALC_MAGS
  int n_bands = run_globals->photo.NBands;
  double LumDust[n_bands];
  double Hubble_h = run_globals->params.Hubble_h;
  for (int ii = 0; ii < n_bands; ii++)
  {
    galout->Mag[ii] = (float)(lum_to_mag(gal.Lum[ii][i_snap])) - 5.0 * log10(Hubble_h);
    apply_dust(n_bands, gal, LumDust, i_snap);
    galout->MagDust[ii] = (float)(lum_to_mag(LumDust[ii])) - 5.0 * log10(Hubble_h);
  }
#else
  return;
#endif
}

#ifdef CALC_MAGS
static int inline phototab_index(
  phototabs_t *photo,
  int          i_band,
  int          i_metal,
  int          i_age)
{
  return (int)((i_age) + photo->NAges * ((i_metal) + photo->NMetals * (i_band)));
}

static void init_jump_index(run_globals_t *run_globals)
{
  // This function precomputes a jump table that allows us to quickly jump to
  // the nearest AgeTab index for any given age.  The larger the NJumps, the
  // more accurate we can be...

  float age;
  int idx;
  float jumpfac;
  phototabs_t *photo = &(run_globals->photo);
  int *jumptab       = photo->JumpTable;
  float *AgeTab      = photo->Ages;

  jumpfac = N_PHOTO_JUMPS / (AgeTab[photo->NAges - 1] - AgeTab[1]);

  for (int ii = 0; ii < N_PHOTO_JUMPS; ii++)
  {
    age = AgeTab[1] + ii / jumpfac;
    idx = 1;
    while (AgeTab[idx + 1] < age)
      idx++;
    jumptab[ii] = idx;
  }

  photo->JumpFactor = jumpfac;
}
#endif

#if defined(DEBUG) && defined(CALC_MAGS)
static void print_phototab(run_globals_t *run_globals, int i_metal)
{
  phototabs_t *photo = &(run_globals->photo);
  float *phototab    = photo->Table;
  int n_ages         = photo->NAges;
  int n_bands        = photo->NBands;

  printf("--------------------\n");
  printf("PHOTOTAB i_metal=%d\n", i_metal);
  for (int i_age = 0; i_age < n_ages; i_age++)
  {
    for (int i_band = 0; i_band < n_bands; i_band++)
      printf("    %.2f", phototab[phototab_index(photo, i_band, i_metal, i_age)]);
    printf("\n");
  }
  printf("--------------------\n");
}
#endif


void read_photometric_tables(run_globals_t *run_globals)
{
#ifdef CALC_MAGS
  run_params_t *run_params = &(run_globals->params);

  // malloc the phototabs struct (see cleanup_mags() for free)
  run_globals->photo       = SID_malloc(sizeof(phototabs_t));
  phototabs_t *photo       = &(run_globals->photo);

  float **Metals           = &(photo->Metals);
  float **AgeTab           = &(photo->Ages);
  char (**MagBands)[5] = &(photo->MagBands);
  float **PhotoTab    = &(photo->Table);
  int n_table_entries = 0;


  if (SID.My_rank == 0)
  {
    hid_t fin;
    hid_t group;
    hsize_t dims[1];
    char name[STRLEN];
    char *bp;
    int i_group = 0;
    float *table_ds;
    int start_ind                = 0;
    double UnitTime_in_Megayears = run_globals->units.UnitTime_in_Megayears;
    double Hubble_h              = run_params->Hubble_h;
    char temp[STRLEN];
    H5E_auto2_t old_func;
    void *old_client_data;
    herr_t error_stack = 0;


    SID_log("Reading photometric tables...", SID_LOG_OPEN);

    // Open the file
    sprintf(name, "%s/%s/%s/%s.hdf5", run_params->PhotometricTablesDir, run_params->SSPModel, run_params->IMF, run_params->MagSystem);
    SID_log("Using tables in file %s", SID_LOG_COMMENT, name);
    fin = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Read the list of metallicities
    SID_log("Reading metallicities...", SID_LOG_COMMENT);
    H5LTget_dataset_info(fin, "metallicities", dims, NULL, NULL);
    photo->NMetals = (int)(dims[0]);
    *Metals        = (float*)SID_malloc(sizeof(float) * (size_t)(dims[0]));
    H5LTread_dataset_float(fin, "metallicities", *Metals);

    // Read the list of ages
    SID_log("Reading ages...", SID_LOG_COMMENT);
    H5LTget_dataset_info(fin, "age", dims, NULL, NULL);
    photo->NAges = (int)(dims[0]);
    *AgeTab      = (float*)SID_malloc(sizeof(float) * (size_t)(dims[0]));
    H5LTread_dataset_float(fin, "age", *AgeTab);

    // Convert the ages from Myr to log10(internal time units)
    for (int ii = 0; ii < photo->NAges; ii++)
      (*AgeTab)[ii] = log10((*AgeTab)[ii] / UnitTime_in_Megayears * Hubble_h);

    // Temporarily turn off hdf5 errors
    H5Eget_auto(error_stack, &old_func, &old_client_data);
    H5Eset_auto(error_stack, NULL, NULL);

    // Parse the requested magnitude bands string and count the number of bands
    // we are going to use
    i_group = 0;
    sprintf(temp, "%s", run_params->MagBands);
    bp = strtok(temp, ",");
    while (bp != NULL)
    {
      group = H5Gopen(fin, bp, H5P_DEFAULT);
      if (group > 0)
      {
        i_group++;
        H5Gclose(group);
      }
      else
        SID_log_warning("Requested magnitude band `%s' not preset in input photometric tables - skipping...", SID_LOG_COMMENT, bp);

      bp = strtok(NULL, " ,\n");
    }
    photo->NBands = i_group;

    if (i_group > MAX_PHOTO_NBANDS)
    {
      SID_log_error("Requested number of valid magnitude bands exceeds maximum (%d > %d)!", SID_LOG_COMMENT, i_group, MAX_PHOTO_NBANDS);
      ABORT(EXIT_FAILURE);
    }
    else if (i_group == 0)
    {
      SID_log("No valid magnitude bands requested!", SID_LOG_COMMENT);
      SID_log("Exiting... Please recompile with CALC_MAGS=0...", SID_LOG_COMMENT);
      ABORT(EXIT_FAILURE);
    }

    // Now we have the number of magnitude bands, metallicities and ages so we
    // can malloc the photometric table itself.
    n_table_entries = photo->NBands * photo->NMetals * photo->NAges;
    SID_log("N table entries = %d", SID_LOG_COMMENT, n_table_entries);
    *PhotoTab = (float*)SID_malloc(sizeof(float) * (size_t)n_table_entries);

    // Finally - loop through the requested bands string one more time, save the
    // band name and store the photometric table data
    *MagBands = SID_malloc(sizeof(char[5]) * (size_t)i_group);
    table_ds  = SID_malloc(sizeof(float) * (size_t)(photo->NAges));
    bp        = strtok(run_params->MagBands, ",");
    SID_log("Reading in photometric table...", SID_LOG_COMMENT);
    i_group = 0;
    while (bp != NULL)
    {
      group = H5Gopen(fin, bp, H5P_DEFAULT);
      if (group > 0)
      {
        sprintf(&((*MagBands)[0][i_group]), "%s", bp);
        for (int i_metal = 0; i_metal < (photo->NMetals); i_metal++)
        {
          sprintf(name, "%0.4f", (*Metals)[i_metal]);
          H5LTread_dataset_float(group, name, table_ds);

          start_ind = phototab_index(photo, i_group, i_metal, 0);
          for (int ii = 0; ii < (photo->NAges); ii++)
            (*PhotoTab)[ii + start_ind] = table_ds[ii];
        }
        H5Gclose(group);
        i_group++;
      }
      bp = strtok(NULL, " ,\n");
    }

    // Restore hdf5 error handling
    H5Eset_auto(error_stack, old_func, old_client_data);

    // Deal with any zeros
    for (int ii = 0; ii < n_table_entries; ii++)
      if ((*PhotoTab)[ii] == 0)
        (*PhotoTab)[ii] = 70.0;

    // Convert metallicities to log10 for interpolation
    for (int ii = 0; ii < photo->NMetals; ii++)
      (*Metals)[ii] = log10((*Metals)[ii]);

    // free temp arrays
    SID_free(SID_FARG table_ds);

    // Close the file
    H5Fclose(fin);
  }

  // if using MPI, broadcast the necessary data to all ranks
  SID_Bcast(&(photo->NMetals), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    *Metals = (float*)SID_malloc(sizeof(float) * photo->NMetals);
  SID_Bcast(*Metals, sizeof(float) * photo->NMetals, 0, SID.COMM_WORLD);

  SID_Bcast(&(photo->NAges), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    *AgeTab = (float*)SID_malloc(sizeof(float) * photo->NAges);
  SID_Bcast(*AgeTab, sizeof(float) * photo->NAges, 0, SID.COMM_WORLD);

  SID_Bcast(&(photo->NBands), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    *MagBands = (float*)SID_malloc(sizeof(char[5]) * photo->NBands);
  SID_Bcast(*MagBands, sizeof(char[5]) * photo->NBands, 0, SID.COMM_WORLD);

  SID_Bcast(&(n_table_entries), sizeof(int), 0, SID.COMM_WORLD);
  if (SID.My_rank > 0)
    *PhotoTab = (float*)SID_malloc(sizeof(float) * n_table_entries);
  SID_Bcast(*PhotoTab, sizeof(float) * n_table_entries 0, SID.COMM_WORLD);

#ifdef DEBUG
  print_phototab(run_globals, 0);
#endif

  init_jump_index(run_globals);

  SID_log(" ...done", SID_LOG_CLOSE);

#else
  return;
#endif // CALC_MAGS
}

#ifdef CALC_MAGS
static int inline get_jump_index(double age, float *AgeTab, int *jumptab, float jumpfac)
{
  return jumptab[(int)((age - AgeTab[1]) * jumpfac)];
}

static void find_interpolated_lum(
  run_globals_t *run_globals,
  double         timenow,
  double         timetarget,
  double         metallicity,
  int           *metals_ind,
  int           *age_ind,
  double        *fage1,
  double        *fage2,
  double        *fmet1,
  double        *fmet2)
{
  // TODO: There is a lot of float/double calculations here. I should tidy this up...

  int k, i, idx;
  float age, frac;
  float fa1, fa2, fm1, fm2;

  phototabs_t *photo = &(run_globals->photo);
  float *Metals      = photo->Metals;
  float *AgeTab      = photo->Ages;
  int *JumpTable     = photo->JumpTable;
  float JumpFactor   = photo->JumpFactor;
  int n_ages         = photo->NAges;
  int n_metals       = photo->NMetals;

  age = (float)(timenow - timetarget);

  if (age > 0)
  {
    age = log10(age);

    if (age > AgeTab[n_ages - 1])         // beyond table, take latest entry
    {
      k   = n_ages - 2;
      fa1 = 0;
      fa2 = 1;
    }
    else if (age < AgeTab[1])     // age younger than 1st enty, take 1st entry
    {
      k   = 0;
      fa1 = 0;
      fa2 = 1;
    }
    else
    {
      idx = get_jump_index(age, AgeTab, JumpTable, JumpFactor);
      while (AgeTab[idx + 1] < age)
        idx++;
      k    = idx;
      frac = (age - AgeTab[idx]) / (AgeTab[idx + 1] - AgeTab[idx]);
      fa1  = 1 - frac;
      fa2  = frac;
    }
  }
  else                           // this lies in the past
  {
    k   = 0;
    fa1 = 0;
    fa2 = 0;
  }

  // Now interpolate also for the metallicity
  metallicity = log10(metallicity);

  if (metallicity > Metals[n_metals - 1])         // beyond table, take latest entry
  {
    i   = n_metals - 2;
    fm1 = 0;
    fm2 = 1;
  }
  else if (metallicity < Metals[0])       // age younger than 1st enty, take 1st entry
  {
    i   = 0;
    fm1 = 1;
    fm2 = 0;
  }
  else
  {
    idx = 0;
    while (Metals[idx + 1] < metallicity)
      idx++;
    i    = idx;
    frac = (metallicity - Metals[idx]) / (Metals[idx + 1] - Metals[idx]);
    fm1  = 1 - frac;
    fm2  = frac;
  }

  *metals_ind = i;
  *age_ind    = k;

  *fage1 = (double)fa1;
  *fage2 = (double)fa2;
  *fmet1 = (double)fm1;
  *fmet2 = (double)fm2;
}
#endif


void add_to_luminosities(
  run_globals_t *run_globals,
  galaxy_t      *gal,
  double         burst_mass,
  double         metallicity,
  double         burst_time)
{
#ifdef CALC_MAGS
  phototabs_t *photo = &(run_globals->photo);
  double Hubble_h    = run_globals->params.Hubble_h;
  float *PhotoTab    = run_globals->photo.Table;
  int n_bands        = run_globals->photo.NBands;
  int metals_ind;
  int age_ind;
  double X1, X2;
  double f1, f2, fmet1, fmet2;

  // Convert burst_mass into 1e6 Msol/(h=Hubble_h) units
  X1 = burst_mass * 1e4 / Hubble_h;
  X2 = -0.4 * M_LN10;

  for (int outputbin = 0; outputbin < NOUT; outputbin++)
  {
    find_interpolated_lum(run_globals, burst_time, run_globals->LTTime[run_globals->ListOutputSnaps[outputbin]], metallicity,
                          &metals_ind, &age_ind, &f1, &f2, &fmet1, &fmet2);

    for (int i_band = 0; i_band < n_bands; i_band++)
      gal->Lum[i_band][outputbin] +=
        X1 * exp(X2 * (fmet1 * (f1 *
                                PhotoTab[phototab_index(photo, i_band, metals_ind, age_ind)]
                                + f2 *
                                PhotoTab[phototab_index(photo, i_band, metals_ind, age_ind + 1)]
                                ) + fmet2 * (f1 *
                                             PhotoTab[phototab_index(photo, i_band, metals_ind + 1, age_ind)]
                                             + f2 *
                                             PhotoTab[phototab_index(photo, i_band, metals_ind + 1, age_ind + 1)]
                                             )));
  }

#else
  return;
#endif // CALC_MAGS
}

double lum_to_mag(double lum)
{
  if (lum > 0)
    return -2.5 * log10(lum);
  else
    return 99.0;
}

void cleanup_mags(run_globals_t *run_globals)
{
#ifdef CALC_MAGS
  SID_free(SID_FARG run_globals->photo.Table);
  SID_free(SID_FARG run_globals->photo.MagBands);
  SID_free(SID_FARG run_globals->photo.Ages);
  SID_free(SID_FARG run_globals->photo.Metals);
  SID_free(SID_FARG run_globals->photo);
  H5Tclose(run_globals->hdf5props.array_nmag_f_tid);
#else
  return;
#endif
}
