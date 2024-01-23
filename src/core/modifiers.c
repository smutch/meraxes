#include <hdf5_hl.h>
#include <string.h>

#include "meraxes.h"
#include "modifiers.h"

void read_mass_ratio_modifiers(int snapshot)
{
  if (strlen(run_globals.params.MassRatioModifier) == 0) {
    run_globals.RequestedMassRatioModifier = -1;
    mlog("No Mass Ratio Modifier.", MLOG_MESG);
  } else {
    run_globals.mass_ratio_modifier = malloc(sizeof(Modifier) * N_LOGMS);
    const size_t dst_size = sizeof(Modifier);
    const size_t dst_sizes[NFIELDS] = {
      sizeof(run_globals.mass_ratio_modifier[0].logMmin),    sizeof(run_globals.mass_ratio_modifier[0].logMmax),
      sizeof(run_globals.mass_ratio_modifier[0].mass_mean),  sizeof(run_globals.mass_ratio_modifier[0].mass_errl),
      sizeof(run_globals.mass_ratio_modifier[0].mass_erru),  sizeof(run_globals.mass_ratio_modifier[0].ratio),
      sizeof(run_globals.mass_ratio_modifier[0].ratio_errl), sizeof(run_globals.mass_ratio_modifier[0].ratio_erru)
    };

    const size_t dst_offset[NFIELDS] = { HOFFSET(Modifier, logMmin),    HOFFSET(Modifier, logMmax),
                                         HOFFSET(Modifier, mass_mean),  HOFFSET(Modifier, mass_errl),
                                         HOFFSET(Modifier, mass_erru),  HOFFSET(Modifier, ratio),
                                         HOFFSET(Modifier, ratio_errl), HOFFSET(Modifier, ratio_erru) };

    if (run_globals.mpi_rank == 0) {
      hid_t fd;
      char fname[STRLEN];
      char tablename[STRLEN];

      sprintf(fname, "%s", run_globals.params.MassRatioModifier);
      fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      sprintf(tablename, "%03d", snapshot);

      H5TBread_fields_name(
        fd,
        tablename,
        "log10(mass_lower),log10(mass_upper),mass_mean,mass_errl,mass_erru,ratio_mean,ratio_errl,ratio_erru",
        N_START,
        N_LOGMS,
        dst_size,
        dst_offset,
        dst_sizes,
        run_globals.mass_ratio_modifier);
      H5Fclose(fd);
    }
    MPI_Bcast(
      run_globals.mass_ratio_modifier, sizeof(run_globals.mass_ratio_modifier), MPI_BYTE, 0, run_globals.mpi_comm);
  }
}

void read_baryon_frac_modifiers(int snapshot)
{
  if (strlen(run_globals.params.BaryonFracModifier) == 0) {
    run_globals.RequestedBaryonFracModifier = -1;
    mlog("No Baryon Fraction Modifier.", MLOG_MESG);
  } else {
    run_globals.baryon_frac_modifier = malloc(sizeof(Modifier) * N_LOGMS);
    const size_t dst_size = sizeof(Modifier);
    const size_t dst_sizes[NFIELDS] = {
      sizeof(run_globals.baryon_frac_modifier[0].logMmin),    sizeof(run_globals.baryon_frac_modifier[0].logMmax),
      sizeof(run_globals.baryon_frac_modifier[0].mass_mean),  sizeof(run_globals.baryon_frac_modifier[0].mass_errl),
      sizeof(run_globals.baryon_frac_modifier[0].mass_erru),  sizeof(run_globals.baryon_frac_modifier[0].ratio),
      sizeof(run_globals.baryon_frac_modifier[0].ratio_errl), sizeof(run_globals.baryon_frac_modifier[0].ratio_erru)
    };

    const size_t dst_offset[NFIELDS] = { HOFFSET(Modifier, logMmin),    HOFFSET(Modifier, logMmax),
                                         HOFFSET(Modifier, mass_mean),  HOFFSET(Modifier, mass_errl),
                                         HOFFSET(Modifier, mass_erru),  HOFFSET(Modifier, ratio),
                                         HOFFSET(Modifier, ratio_errl), HOFFSET(Modifier, ratio_erru) };

    if (run_globals.mpi_rank == 0) {
      hid_t fd;
      char fname[STRLEN];
      char tablename[STRLEN];

      sprintf(fname, "%s", run_globals.params.BaryonFracModifier);
      fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      sprintf(tablename, "%03d", snapshot);

      H5TBread_fields_name(fd,
                           tablename,
                           "log10(mass_lower),log10(mass_upper),mass_mean,mass_errl,mass_erru,fb_mean,fb_errl,fb_erru",
                           N_START,
                           N_LOGMS,
                           dst_size,
                           dst_offset,
                           dst_sizes,
                           run_globals.baryon_frac_modifier);
      H5Fclose(fd);
    }
    MPI_Bcast(
      run_globals.baryon_frac_modifier, sizeof(run_globals.baryon_frac_modifier), MPI_BYTE, 0, run_globals.mpi_comm);
  }
}

double interpolate_modifier(Modifier* modifier_data, double logM)
{
  if (logM <= modifier_data[0].logMmin + M_OFFSET)
    return (double)modifier_data[0].ratio;

  if (logM >= modifier_data[N_LOGMS - 1].logMmin + M_OFFSET)
    return (double)modifier_data[N_LOGMS - 1].ratio;

  double logM_below, ratio_below, ratio_above, ratio;
  int i;

  i = 0;
  while (logM > modifier_data[i].logMmin + M_OFFSET)
    i++;

  logM_below = (double)modifier_data[i - 1].logMmin + M_OFFSET;
  ratio_below = (double)modifier_data[i - 1].ratio;
  ratio_above = (double)modifier_data[i].ratio;

  ratio = ratio_below + (ratio_above - ratio_below) / DELTA_M * (logM - logM_below);

  if (ratio <= 0) {
    mlog_error("Something wrong about the modifier: logM=%f, ratio=%f", logM, ratio);
    ABORT(EXIT_FAILURE);
  }

  return ratio;
}