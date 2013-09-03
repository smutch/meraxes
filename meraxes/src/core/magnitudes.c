#include <math.h>
#include "meraxes.h"

void init_luminosities(galaxy_struct *gal)
{
#ifdef CALC_MAGS
  for(int ii=0; ii<NOUT; ii++)
    for(int jj=0; jj<N_PHOTO_BANDS; jj++)
      gal->Lum[jj][ii]     = 0.0;
#else
  return;
#endif
}

void sum_luminosities(galaxy_struct *parent, galaxy_struct *gal, int outputbin)
{
#ifdef CALC_MAGS
  for(int ii=0; ii < N_PHOTO_BANDS; ii++)
    parent->Lum[ii][outputbin]     += gal->Lum[ii][outputbin];
#else
  return;
#endif
}

void prepare_magnitudes_for_output(galaxy_struct gal, galaxy_output_struct *galout, int i_snap)
{
#ifdef CALC_MAGS
  double LumDust[N_PHOTO_BANDS];
  for(int ii=0; ii<N_PHOTO_BANDS; ii++)
  {
    galout->Mag[ii] = (float)(lum_to_mag(gal.Lum[ii][i_snap]));
    apply_dust(gal, LumDust, i_snap);
    galout->MagDust[ii] = (float)(lum_to_mag(LumDust[ii]));
  }
#else
  return;
#endif
}

static int inline phototab_index(
  phototabs_struct   *photo,
  int                 i_band,     
  int                 i_metal,    
  int                 i_age)      
{
  return (int)((i_age)+photo->NAges*((i_metal)+photo->NMetals*(i_band)));
}

static void init_jump_index(run_globals_struct *run_globals)
{

  // This function precomputes a jump table that allows us to quickly jump to
  // the nearest AgeTab index for any given age.  The larger the NJumps, the
  // more accurate we can be...

  float  age;
  int    idx;
  int   *jumptab = run_globals->photo.JumpTable;
  float  jumpfac;
  float *AgeTab  = run_globals->photo.Ages;

  jumpfac = N_PHOTO_JUMPS / (AgeTab[N_PHOTO_AGES - 1] - AgeTab[1]);

  for(int ii = 0; ii < N_PHOTO_JUMPS; ii++)
  {
    age = AgeTab[1] + ii / jumpfac;
    idx = 1;
    while(AgeTab[idx + 1] < age)
      idx++;
    jumptab[ii] = idx;
  }

  run_globals->photo.JumpFactor = jumpfac;

}


void read_photometric_tables(run_globals_struct *run_globals)
{
#ifdef CALC_MAGS

  h5id_t              fin;
  h5id_t              group;
  hsize_t             dims[1];
  char                name[STRLEN];
  run_params_struct  *run_params   = &(run_globals->params);
  phototabs_struct   *photo        = &(run_globals->photo);
  float              *Metals       = photo->Metals;
  float              *AgeTab       = photo->Ages;
  char              *(MagBands[5]    ) = photo->MagBands;
  float              *PhotoTab     = photo->Table;
  char               *bp;
  int                 i_group      = 0;
  float              *table_ds;
  int                 start_ind    = 0;
  int                 n_table_entries = 0;


  SID_log("Reading photometric tables...", SID_LOG_OPEN);

  // Open the file
  sprintf(name, "%s/%s/%s/%s.hdf5", run_params->PhotometricTablesDir, run_params->SSPModel, run_params->IMF, run_params->MagSystem);
  fin = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read the list of metallicities
  H5LTget_dataset_info(fin, "metallicities", dims, NULL, NULL);
  photo->NMetals = (int)(dim[0]);
  Metals = (float *)SID_malloc(sizeof(float) * (size_t)(dims[0]));
  H5LTread_dataset_float(fin, "metallicities", Metals);

  // Read the list of ages
  bp = strtok(run_params->MagBands, ",");
  group = H5Gopen(fin, bp, H5P_DEFAULT);
  sprintf(name, "%0.4f", Metals[0]);
  H5LTget_dataset_info(group, name, dims, NULL, NULL);
  H5Gclose(group);
  photo->NAges = (int)(dim[0]);
  H5LTread_dataset_float(fin, "ages", AgeTab);

  // Convert the ages from Myr to log10(internal time units)
  for(int ii=0; ii<photo->NAges; ii++)
    AgeTab[ii] = log10(AgeTab[ii] / UnitTime_in_Megayears * Hubble_h);

  // Parse the requested magnitude bands string and count the number of bands
  // we are going to use
  i_group = 0;
  bp = strtok(run_params->MagBands, ",");
  while (bp!=NULL){

    if(group = H5Gopen(fin, bp, H5P_DEFAULT))
    {
      i_group++;
      H5Gclose(group);
    } else
      SID_log_warning("Requested magnitude band %s not preset in input photometric tables - skipping...", SID_LOG_COMMENT);

    bp = strtok(NULL, " ,\n");
  }
  photo->NMagBands = i_group;

  // Now we have the number of magnitude bands, metallicities and ages so we
  // can malloc the photometric table itself.
  n_table_entries = photo->NMagBands * photo->NMetals * photo->NAges;
  PhotoTab = (float *)SID_malloc(sizeof(float) * (size_t)n_entries);

  // Finally - loop through the requested bands string one more time, save the
  // band name and store the potometric table data
  MagBands = (char[5] *)SID_malloc(sizeof(char[5]), (size_t)i_group);
  table_ds = SID_malloc(sizeof(float) * (size_t)(photo->NAges));
  bp = strtok(run_params->MagBands, ",");
  while (bp!=NULL){
    if(group = H5Gopen(fin, bp, H5P_DEFAULT))
    {

      sprintf(&(MagBands[i_group]), "%s", bp);
      for(int i_metal=0; i_metal<photo->NMetals; i_metal++)
      {
      sprintf(name, "%0.4f", Metals[i_metal]);
      H5LTread_dataset_float(group, name, table_ds);

      start_ind = phototab_index(photo, i_group, i_metal, 0);
      for(int ii=0; ii<(photo->NAges); ii++)
        PhotoTab[ii] = table_ds[ii+start_ind];

      }
      H5Gclose(group);
      i_group++;
    }
    bp = strtok(NULL, " ,\n");
  }
  SID_free(SID_FARG table_ds);

  // Deal with any zeros
  for(int ii=0; ii<n_table_entries; ii++)
    if(PhotoTab[ii] == 0)
      PhotoTab[ii] = 70.0;

  // Convert metallicities to log10 for interpolation  
  for(int ii=0; ii<photo->NMetals; ii++)
    Metals[ii] = log10(Metals[ii]);

  // Close the file
  H5Fclose(fin);

  SID_log(" ...done", SID_LOG_CLOSE);

#else
  return;
#endif // CALC_MAGS
}

static int inline get_jump_index(double age, float *AgeTab, int *jumptab, float jumpfac)
{
  return jumptab[(int) ((age - AgeTab[1]) * jumpfac)];
}

static void find_interpolated_lum(
  run_globals_struct *run_globals,
  double              timenow,    
  double              timetarget, 
  double              metallicity,
  int                *metals_ind, 
  int                *age_ind,    
  double             *fage1,      
  double             *fage2,      
  double             *fmet1,      
  double             *fmet2)      
{

  // TODO: There is a lot of float/double caclulations here. I should tidy this up...

  int k, i, idx;
  float age, frac;
  float fa1, fa2, fm1, fm2;

  float *Metals     = run_globals->photo.Metals;
  float *AgeTab     = run_globals->photo.Ages;
  int   *JumpTable  = run_globals->photo.JumpTable;
  float  JumpFactor = run_globals->photo.JumpFactor;

  age = (float)(timenow - timetarget);

  if(age > 0)
  {
    age = log10(age);

    if(age > AgeTab[N_PHOTO_AGES - 1])	 // beyond table, take latest entry 
    {
      k = N_PHOTO_AGES - 2;
      fa1 = 0;
      fa2 = 1;
    }
    else if(age < AgeTab[1])	 // age younger than 1st enty, take 1st entry 
    {
      k = 0;
      fa1 = 0;
      fa2 = 1;
    }
    else
    {
      idx = get_jump_index(age, AgeTab, JumpTable, JumpFactor);
      while(AgeTab[idx + 1] < age)
        idx++;
      k = idx;
      frac = (age - AgeTab[idx]) / (AgeTab[idx + 1] - AgeTab[idx]);
      fa1 = 1 - frac;
      fa2 = frac;
    }
  }
  else				 // this lies in the past 
  {
    k = 0;
    fa1 = 0;
    fa2 = 0;
  }

  // Now interpolate also for the metallicity 
  metallicity = log10(metallicity);

  if(metallicity > Metals[N_PHOTO_METALS - 1])	 // beyond table, take latest entry 
  {
    i = N_PHOTO_METALS - 2;
    fm1 = 0;
    fm2 = 1;
  }
  else if(metallicity < Metals[0])	 // age younger than 1st enty, take 1st entry 
  {
    i = 0;
    fm1 = 1;
    fm2 = 0;
  }
  else
  {
    idx = 0;
    while(Metals[idx + 1] < metallicity)
      idx++;
    i = idx;
    frac = (metallicity - Metals[idx]) / (Metals[idx + 1] - Metals[idx]);
    fm1 = 1 - frac;
    fm2 = frac;
  }

  *metals_ind = i;
  *age_ind = k;

  *fage1 = (double)fa1;
  *fage2 = (double)fa2;
  *fmet1 = (double)fm1;
  *fmet2 = (double)fm2;
}


void add_to_luminosities(
  run_globals_struct *run_globals,
  galaxy_struct      *gal,        
  double              burst_mass, 
  double              metallicity,
  double              burst_time) 
{
#ifdef CALC_MAGS
  double  X1             , X2;
  double  Hubble_h     = run_globals->params.Hubble_h;
  float  *PhotoTab     = run_globals->photo.Table;
  int     metals_ind;
  int     age_ind;
  double  f1, f2, fmet1, fmet2;

  // Convert burst_mass into 1e11 Msol/(h=Hubble_h) units
  X1 = burst_mass * 0.1 / Hubble_h;
  X2 = -0.4 * M_LN10;

  for(int outputbin = 0; outputbin < NOUT; outputbin++)
  {

    find_interpolated_lum(run_globals, burst_time, run_globals->LTTime[run_globals->ListOutputSnaps[outputbin]], metallicity,
        &metals_ind, &age_ind, &f1, &f2, &fmet1, &fmet2);

    // NB - tables give luminosities for a 1.0e^11 M_sun burst 
    for(int i_band = 0; i_band < N_PHOTO_BANDS; i_band++)
      gal->Lum[i_band][outputbin] +=
        X1 * exp(X2 * (fmet1 * (f1 *  PhotoTab[phototab_index(i_band,
                  metals_ind, age_ind)]+ f2 * PhotoTab[phototab_index(i_band,
                    metals_ind, age_ind+1)]) + fmet2 * (f1 *
                    PhotoTab[phototab_index(i_band, metals_ind+1, age_ind)] +
                    f2 * PhotoTab[phototab_index(i_band, metals_ind+1,
                      age_ind+1)])));
  }

#else
  return;
#endif // CALC_MAGS
}

double lum_to_mag(double lum)
{
  
  if(lum > 0)
    return -2.5 * log10(lum);
  else
    return 99.0;
}

void cleanup_mags(run_globals_struct *run_globals)
{
#ifdef CALC_MAGS
  SID_free(SID_FARG run_globals->photo.PhotoTab);
  SID_free(SID_FARG run_globals->photo.MagBands);
  SID_free(SID_FARG run_globals->photo.Ages);
  SID_free(SID_FARG run_globals->photo.Metals);
#else
  return
#endif
}
