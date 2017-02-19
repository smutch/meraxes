#include "meraxes.h"
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>

// Qin, Y. et al., 2017. Dark-ages Reionization and Galaxy Formation Simulation VIII.
// Suppressed growth of dark matter halos during the Epoch of Reionization.
// Monthly Notices of the Royal Astronomical Society, 9(November), p.stx083.
// Available at: http://mnras.oxfordjournals.org/lookup/doi/10.1093/mnras/stx083.
// DO NOT CHANGE THESE UNLESS YOU KNOW WHAT THESE MEAN
#define NFIELDS 8
#define N_START 0
#define N_LOGMS 31
#define DELTA_M 1.0
#define MIN_LOGM 7.5
#define MAX_LOGM 11.5
#define M_OFFSET 0.5

void read_mass_ratio_modifiers(int snapshot)
{
  if (strlen(run_globals.params.MassRatioModifier) == 0)
  {
    run_globals.RequestedMassRatioModifier = -1;
    SID_log("No Mass Ratio Modifier :(", SID_LOG_COMMENT);
  }
  else{
    run_globals.mass_ratio_modifier = SID_malloc(sizeof(Modifier) * N_LOGMS);
    const size_t dst_size           = sizeof(Modifier);
    const size_t dst_sizes[NFIELDS] = {
      sizeof(run_globals.mass_ratio_modifier[0].logMmin),
      sizeof(run_globals.mass_ratio_modifier[0].logMmax),
      sizeof(run_globals.mass_ratio_modifier[0].mass_mean),
      sizeof(run_globals.mass_ratio_modifier[0].mass_errl),
      sizeof(run_globals.mass_ratio_modifier[0].mass_erru),
      sizeof(run_globals.mass_ratio_modifier[0].ratio),
      sizeof(run_globals.mass_ratio_modifier[0].ratio_errl),
      sizeof(run_globals.mass_ratio_modifier[0].ratio_erru)
    };

    const size_t dst_offset[NFIELDS] = {
      HOFFSET(Modifier, logMmin),
      HOFFSET(Modifier, logMmax),
      HOFFSET(Modifier, mass_mean),
      HOFFSET(Modifier, mass_errl),
      HOFFSET(Modifier, mass_erru),
      HOFFSET(Modifier, ratio),
      HOFFSET(Modifier, ratio_errl),
      HOFFSET(Modifier, ratio_erru)
    };

    if (SID.My_rank == 0)
    {
      hid_t fd;
      char  fname[STRLEN];
      char  tablename[STRLEN];

      sprintf(fname, "%s", run_globals.params.MassRatioModifier);
      fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      sprintf(tablename, "%03d", snapshot);

      H5TBread_fields_name(fd, tablename,"log10(mass_lower),log10(mass_upper),mass_mean,mass_errl,mass_erru,ratio_mean,ratio_errl,ratio_erru",
                           N_START, N_LOGMS, dst_size, dst_offset, dst_sizes, run_globals.mass_ratio_modifier);
      H5Fclose(fd);
    }
    SID_Bcast(run_globals.mass_ratio_modifier, sizeof(run_globals.mass_ratio_modifier), 0, SID.COMM_WORLD);
  }
}


void read_baryon_frac_modifiers(int snapshot)
{
  if (strlen(run_globals.params.BaryonFracModifier) == 0)
  {
    run_globals.RequestedBaryonFracModifier = -1;
    SID_log("No Baryon Fraction Modifier :(", SID_LOG_COMMENT);
  }
  else{
    run_globals.baryon_frac_modifier = SID_malloc(sizeof(Modifier) * N_LOGMS);
    const size_t dst_size           = sizeof(Modifier);
    const size_t dst_sizes[NFIELDS] = {
      sizeof(run_globals.baryon_frac_modifier[0].logMmin),
      sizeof(run_globals.baryon_frac_modifier[0].logMmax),
      sizeof(run_globals.baryon_frac_modifier[0].mass_mean),
      sizeof(run_globals.baryon_frac_modifier[0].mass_errl),
      sizeof(run_globals.baryon_frac_modifier[0].mass_erru),
      sizeof(run_globals.baryon_frac_modifier[0].ratio),
      sizeof(run_globals.baryon_frac_modifier[0].ratio_errl),
      sizeof(run_globals.baryon_frac_modifier[0].ratio_erru)
    };

    const size_t dst_offset[NFIELDS] = {
      HOFFSET(Modifier, logMmin),
      HOFFSET(Modifier, logMmax),
      HOFFSET(Modifier, mass_mean),
      HOFFSET(Modifier, mass_errl),
      HOFFSET(Modifier, mass_erru),
      HOFFSET(Modifier, ratio),
      HOFFSET(Modifier, ratio_errl),
      HOFFSET(Modifier, ratio_erru)
    };

    if (SID.My_rank == 0)
    {
      hid_t fd;
      char  fname[STRLEN];
      char  tablename[STRLEN];

      sprintf(fname, "%s", run_globals.params.BaryonFracModifier);
      fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      sprintf(tablename, "%03d", snapshot);

      H5TBread_fields_name(fd, tablename,"log10(mass_lower),log10(mass_upper),mass_mean,mass_errl,mass_erru,fb_mean,fb_errl,fb_erru",
                           N_START, N_LOGMS, dst_size, dst_offset, dst_sizes, run_globals.baryon_frac_modifier);
      H5Fclose(fd);
    }
    SID_Bcast(run_globals.baryon_frac_modifier, sizeof(run_globals.baryon_frac_modifier), 0, SID.COMM_WORLD);
  }
}


double interpolate_modifier(Modifier *modifier_data, double logM)
{
  if (logM < modifier_data[0].logMmin)
    return modifier_data[0].ratio;

  if (logM > modifier_data[N_LOGMS - 1].logMmin)
    return modifier_data[N_LOGMS - 1].ratio;

  double logM_below, ratio_below, ratio_above, ratio;
  int    i;

  i = 0;
  while (logM > modifier_data[i].logMmin)
    i++;

  logM_below  = modifier_data[i].logMmin;
  ratio_below = modifier_data[i].ratio;
  ratio_above = modifier_data[i + 1].ratio;

  ratio       = ratio_below + (ratio_above - ratio_below) / DELTA_M * (logM - logM_below);

  return ratio;
}