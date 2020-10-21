#include <fftw3-mpi.h>

#include "magnitudes.h"
#include "meraxes.h"
#include "parse_paramfile.h"
#include "read_grids.h"
#include "read_halos.h"
#include "recombinations.h"
#include "reionization.h"

void cleanup()
{
  mlog("Running cleanup...", MLOG_OPEN);

  free_grids_cache();

  if (run_globals.RequestedMassRatioModifier != -1)
    free(run_globals.mass_ratio_modifier);
  if (run_globals.RequestedBaryonFracModifier != -1)
    free(run_globals.baryon_frac_modifier);

  free_halo_storage();

#ifdef CALC_MAGS
  cleanup_mags();
#endif

  if (run_globals.RequestedForestId)
    free(run_globals.RequestedForestId);

  if (run_globals.params.Flag_PatchyReion) {
    free_reionization_grids();
    fftwf_mpi_cleanup();
  }

  if (run_globals.params.Flag_IncludeRecombinations) {
    free_MHR();
  }

  if (!run_globals.params.FlagMCMC) {
    mlog("Freeing hdf5 related stuff...", MLOG_OPEN);
    if (run_globals.mpi_rank == 0) {
      free(run_globals.hdf5props.params_addr);
      free(run_globals.hdf5props.params_type);
      for (int ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
        free(run_globals.hdf5props.params_tag[ii]);
      free(run_globals.hdf5props.params_tag);
    }
    free(run_globals.hdf5props.field_h_conv);
    free(run_globals.hdf5props.field_units);
    free(run_globals.hdf5props.field_types);
    free(run_globals.hdf5props.field_names);
    free(run_globals.hdf5props.dst_field_sizes);
    free(run_globals.hdf5props.dst_offsets);
    H5Tclose(run_globals.hdf5props.array3f_tid);
    H5Tclose(run_globals.hdf5props.array_nhist_f_tid);
    mlog(" ...done", MLOG_CLOSE);
  }

  gsl_rng_free(run_globals.random_generator);

  free(run_globals.ListOutputSnaps);
  free(run_globals.LTTime);
  free(run_globals.ZZ);
  free(run_globals.AA);

  if (run_globals.gpu != NULL)
    free(run_globals.gpu);

  mlog(" ...done\n", MLOG_CLOSE | MLOG_FLUSH);
}
