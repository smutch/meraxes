#include <math.h>
#include "meraxes.h"

static int inline phototab_index(
  int                 i_band,       
  int                 i_metal,
  int                 i_age)        
{
  return (int)((i_age)+N_PHOTO_METALS*((i_metal)+N_PHOTO_BANDS*(i_band)));
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

  // DEBUG
  // SID_log("=================", SID_LOG_COMMENT);
  // SID_log("Jump index table:", SID_LOG_COMMENT);
  // for(int ii=0; ii<N_PHOTO_JUMPS; ii++)
  // {
  //   age = AgeTab[1] + ii / jumpfac;
  //   idx = 1;
  //   while(AgeTab[idx + 1] < age)
  //     idx++;
  //   SID_log("%.3f  ->  %d (%.3f)", SID_LOG_COMMENT, age, idx, AgeTab[idx]);
  // }
  // SID_log("=================", SID_LOG_COMMENT);
}

void read_photometric_tables(run_globals_struct *run_globals)
{
  FILE  *fin;
  char   filename[STRLEN];
  int    i_age;
  float *Metals                = run_globals->photo.Metals;
  float *PhotoTab              = run_globals->photo.Table;
  float *AgeTab                = run_globals->photo.Ages;
  float  dummy;
  float  Hubble_h              = run_globals->params.Hubble_h;
  float  UnitTime_in_Megayears = run_globals->units.UnitTime_in_Megayears;

  SID_log("Reading photometric tables...", SID_LOG_OPEN);

  Metals[0] = 0.0001;  // TODO: This shouldn't be hard coded
  Metals[1] = 0.0004;  // TODO: This shouldn't be hard coded
  Metals[2] = 0.004;  // TODO: This shouldn't be hard coded
  Metals[3] = 0.008;  // TODO: This shouldn't be hard coded
  Metals[4] = 0.02;  // TODO: This shouldn't be hard coded
  Metals[5] = 0.05;  // TODO: This shouldn't be hard coded

  for (int i_Z=0; i_Z<N_PHOTO_METALS; i_Z++)
  {
    sprintf(filename, "%s/Z%.4f_salp.bc03", run_globals->params.PhotometricTablesDir, Metals[i_Z]);

    if(!(fin = fopen(filename, "r")))
    {
      SID_log("file `%s' not found.", SID_LOG_COMMENT, filename);
      ABORT(EXIT_FAILURE);
    }

    i_age = 0;
    while(fscanf(fin, " %f %f %f %f %f %f %f ", 
          &dummy,
          &PhotoTab[phototab_index(0, i_Z, i_age)],
          &PhotoTab[phototab_index(1, i_Z, i_age)],
          &PhotoTab[phototab_index(2, i_Z, i_age)],
          &PhotoTab[phototab_index(3, i_Z, i_age)],
          &PhotoTab[phototab_index(4, i_Z, i_age)],
          &AgeTab[i_age])
        !=EOF)
    {

      // convert AgeTab from log10(Gyr) to log10(internal time units) 
      AgeTab[i_age] = pow(10.0, AgeTab[i_age]) / 1.0e6 /
        UnitTime_in_Megayears * Hubble_h;
      AgeTab[i_age] = log10(AgeTab[i_age]);

      for(int ii=0; ii<N_PHOTO_BANDS; ii++)
        if(PhotoTab[phototab_index(ii, i_Z, i_age)] == 0)
          PhotoTab[phototab_index(ii, i_Z, i_age)] = 70.0;

      i_age++;
    }
    SID_log("Read %d ages from %s", SID_LOG_COMMENT, i_age, filename);
    fclose(fin);
  }

  // DEBUG
  // SID_log("=============================", SID_LOG_COMMENT);
  // SID_log("Checking input table (Z0.02):", SID_LOG_COMMENT);
  // for(int ii=0; ii<N_PHOTO_AGES; ii++)
  //   SID_log("%.3f  %.3f  %.3f  %.3f  %.3f :: %.3e", SID_LOG_COMMENT,
  //         PhotoTab[phototab_index(0, 4, ii)],
  //         PhotoTab[phototab_index(1, 4, ii)],
  //         PhotoTab[phototab_index(2, 4, ii)],
  //         PhotoTab[phototab_index(3, 4, ii)],
  //         PhotoTab[phototab_index(4, 4, ii)],
  //         pow(10.0, AgeTab[ii])/Hubble_h * UnitTime_in_Megayears);
  // SID_log("=============================", SID_LOG_COMMENT);

  init_jump_index(run_globals);

  // Lastly - convert the metallicities into log10 values for interpolation purposes
  for(int ii=0; ii<N_PHOTO_METALS; ii++)
    Metals[ii] = log10(Metals[ii]);

  SID_log(" ...done", SID_LOG_CLOSE);
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
  double  X1             , X2;
  double  Hubble_h     = run_globals->params.Hubble_h;
  float  *PhotoTab     = run_globals->photo.Table;
  int     metals_ind;
  int     age_ind;
  double  f1, f2, fmet1, fmet2;

  // Convert burst_mass into 1e11 Msol/(h=Hubble_h) units
  X1 = burst_mass * 0.1 / Hubble_h;
  X2 = -0.4 * M_LN10;

  // DEBUG
  // float  *AgeTab     = run_globals->photo.Ages;
  // float  *Metals     = run_globals->photo.Metals;
  // SID_log("=============================", SID_LOG_COMMENT);
  // find_interpolated_lum(run_globals, burst_time, run_globals->LTTime[run_globals->ListOutputSnaps[0]], metallicity,
  //     &metals_ind, &age_ind, &f1, &f2, &fmet1, &fmet2);
  // SID_log("lum[0] before = %.3e", SID_LOG_COMMENT, gal->Lum[0][0]);
  // SID_log("lum[1] before = %.3e", SID_LOG_COMMENT, gal->Lum[1][0]);
  // SID_log("target_time = %.3e", SID_LOG_COMMENT, run_globals->LTTime[run_globals->ListOutputSnaps[0]]);
  // SID_log("burst_time = %.3e", SID_LOG_COMMENT, burst_time);
  // SID_log("burst_mass = %.3e", SID_LOG_COMMENT, burst_mass);
  // SID_log("metallicity = %.3e", SID_LOG_COMMENT, metallicity);
  // SID_log("age = %.3e", SID_LOG_COMMENT, (burst_time-run_globals->LTTime[run_globals->ListOutputSnaps[0]])/Hubble_h*run_globals->units.UnitTime_in_Megayears);
  // SID_log("Checking input table (Z%.5f):", SID_LOG_COMMENT, pow(10., Metals[metals_ind]));
  // for(int ii=0; ii<N_PHOTO_AGES; ii++)
  //   SID_log("%.3f  %.3f  %.3f  %.3f  %.3f :: %.3e", SID_LOG_COMMENT,
  //         PhotoTab[phototab_index(0, metals_ind, ii)],
  //         PhotoTab[phototab_index(1, metals_ind, ii)],
  //         PhotoTab[phototab_index(2, metals_ind, ii)],
  //         PhotoTab[phototab_index(3, metals_ind, ii)],
  //         PhotoTab[phototab_index(4, metals_ind, ii)],
  //         pow(10.0, AgeTab[ii])/Hubble_h * run_globals->units.UnitTime_in_Megayears);
  // SID_log("=============================", SID_LOG_COMMENT);
  
  for(int outputbin = 0; outputbin < NOUT; outputbin++)
  {

    find_interpolated_lum(run_globals, burst_time, run_globals->LTTime[run_globals->ListOutputSnaps[outputbin]], metallicity,
        &metals_ind, &age_ind, &f1, &f2, &fmet1, &fmet2);

    // DEBUG
    // SID_log("burst_time = %.2f; burst_mass = %.2f; metals_ind=%d; age_ind=%d; f1=%.2f; f2=%.2f; fmet1=%.2f; fmet2=%.2f", SID_LOG_COMMENT, 
    //     burst_time, burst_mass, metals_ind, age_ind, f1, f2, fmet1, fmet2);

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

  // DEBUG
  // SID_log("lum[0] after = %.3e (mag = %.2f)", SID_LOG_COMMENT, gal->Lum[0][0], lum_to_mag(gal->Lum[0][0]));
  // SID_log("lum[1] after = %.3e (mag = %.2f)", SID_LOG_COMMENT, gal->Lum[1][0], lum_to_mag(gal->Lum[1][0]));
  // if(gal->ID==0)
  //   ABORT(EXIT_SUCCESS);

}

double lum_to_mag(double lum)
{
  
  if(lum > 0)
    return -2.5 * log10(lum);
  else
    return 99.0;
}
