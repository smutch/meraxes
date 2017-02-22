#include "meraxes.h"
#include "parse_paramfile.h"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <hdf5.h>

void cleanup()
{
  SID_log("Running cleanup...", SID_LOG_OPEN);

  free_grids_cache();

  if (run_globals.RequestedMassRatioModifier != -1)
    SID_free(SID_FARG run_globals.mass_ratio_modifier);

  free_halo_storage();

  cleanup_mags();

  if (run_globals.RequestedForestId)
    SID_free(SID_FARG run_globals.RequestedForestId);

  if (run_globals.params.Flag_PatchyReion)
  {
    free_reionization_grids();
    fftwf_mpi_cleanup();
  }

  if (!run_globals.params.FlagMCMC)
  {
    SID_log("Freeing hdf5 related stuff...", SID_LOG_OPEN);
    if (SID.My_rank == 0)
    {
      SID_free(SID_FARG run_globals.hdf5props.params_addr);
      SID_free(SID_FARG run_globals.hdf5props.params_type);
      for (int ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
        SID_free(SID_FARG run_globals.hdf5props.params_tag[ii]);
      SID_free(SID_FARG run_globals.hdf5props.params_tag);
    }
    SID_free(SID_FARG run_globals.hdf5props.field_h_conv);
    SID_free(SID_FARG run_globals.hdf5props.field_units);
    SID_free(SID_FARG run_globals.hdf5props.field_types);
    SID_free(SID_FARG run_globals.hdf5props.field_names);
    SID_free(SID_FARG run_globals.hdf5props.dst_field_sizes);
    SID_free(SID_FARG run_globals.hdf5props.dst_offsets);
    H5Tclose(run_globals.hdf5props.array3f_tid);
    SID_log(" ...done", SID_LOG_CLOSE);
  }

  gsl_rng_free(run_globals.random_generator);

  SID_free(SID_FARG run_globals.ListOutputSnaps);
  SID_free(SID_FARG run_globals.LTTime);
  SID_free(SID_FARG run_globals.ZZ);
  SID_free(SID_FARG run_globals.AA);

  SID_log(" ...done", SID_LOG_CLOSE);
}