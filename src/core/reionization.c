#include <assert.h>
#include <complex.h>
#include <fenv.h>
#include <fftw3-mpi.h>
#include <hdf5_hl.h>
#include <math.h>
#include <sys/stat.h>

#include "ComputeTs.h"
#include "find_HII_bubbles.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "read_grids.h"
#include "reionization.h"
#include "virial_properties.h"

void update_galaxy_fesc_vals(galaxy_t* gal, double new_stars, int snapshot)
{
  physics_params_t* params = &(run_globals.params.physics);

  float fesc_bh = (float)(params->EscapeFracBHNorm *
                          (powf((float)((1.0 + run_globals.ZZ[snapshot]) / 6.0), (float)params->EscapeFracBHScaling)));

  double fesc = params->EscapeFracNorm;

  // redshift
  if ((params->EscapeFracDependency > 0) && (params->EscapeFracDependency <= 6))
    if (params->EscapeFracRedshiftScaling != 0.0)
      fesc *=
        2.0 /
        (1.0 + exp(params->EscapeFracRedshiftScaling * (run_globals.ZZ[snapshot] - params->EscapeFracRedshiftOffset)));

  // galaxy properties
  switch (params->EscapeFracDependency) {
    case 0:
    case 1:
      break;
    case 2: // stellar mass (Msun)
      if (gal->StellarMass > 0.0)
        fesc *= pow((gal->StellarMass / run_globals.params.Hubble_h), params->EscapeFracPropScaling);
      else
        fesc = 1.0;
      break;
    case 3: // star formation rate (Msun / yr)
      if (gal->Sfr > 0.0)
        fesc *=
          pow(gal->Sfr * run_globals.units.UnitMass_in_g / run_globals.units.UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS,
              params->EscapeFracPropScaling);
      else
        fesc = 0.0;
      break;
    case 4: // cold gas density (Msun / pc^2)
      if ((gal->ColdGas > 0.0) && (gal->DiskScaleLength > 0.0))
        fesc *=
          pow((gal->ColdGas / gal->DiskScaleLength / gal->DiskScaleLength * 0.01 * run_globals.params.Hubble_h) / 10.,
              params->EscapeFracPropScaling);
      else
        fesc = 1.0;
      break;
    case 5: // halo mass (1e9 Msun)
      if (gal->Mvir > 0.0)
        fesc *= pow(gal->Mvir * 10. / run_globals.params.Hubble_h, params->EscapeFracPropScaling);
      else
        fesc = 1.0;
      break;
    case 6: // specific star formation rate (1 / Myr)
      if ((gal->Sfr > 0.0) && (gal->StellarMass > 0.0))
        fesc *= pow(gal->Sfr / gal->StellarMass / run_globals.units.UnitTime_in_s * SEC_PER_MEGAYEAR,
                    params->EscapeFracPropScaling);
      else
        fesc = 0.0;
      break;
    default:
      mlog_error("Unrecognised EscapeFracDependency parameter value.");
  }

  if (fesc > 1.0)
    fesc = 1.0;
  else if (fesc < 0.0)
    fesc = 0.0;

  if (fesc_bh > 1.0)
    fesc_bh = 1.0;
  else if (fesc_bh < 0.0)
    fesc_bh = 0.0;

  gal->Fesc = fesc;
  gal->FescWeightedGSM += new_stars * fesc;

  // Here we just set the black hole escape fraction for use in the next
  // timestep when calculating previous_merger_driven_BH_growth().  It is in
  // this function that the EffectiveBHM will be updated.  For the galaxy
  // that we are working on now, the EffectiveBHM has already been calculated
  // using the FescBH value from the previous snapshot.
  // The upshot is that we don't need to do anything to the EffectiveBHM
  // here.  It's confusing I know.  I intend to re-write this to make things
  // more obvious at some point in the future.
  // TODO(smutch): Check this all out and ensure that it is valid for reidentified ghosts
  gal->FescBH = fesc_bh;
}

void set_quasar_fobs()
{
  physics_params_t* params = &(run_globals.params.physics);

  params->quasar_fobs = 1. - cos(params->quasar_open_angle / 180. * M_PI / 2.);
  mlog("Quasar radiation open angle is set to be %g, corresponding to an obscure fraction of %g",
       MLOG_MESG | MLOG_FLUSH,
       params->quasar_open_angle,
       params->quasar_fobs);
}

void set_ReionEfficiency()
{
  // Use the params passed to Meraxes via the input file to set the HII ionising efficiency factor
  physics_params_t* params = &(run_globals.params.physics);

  // The following is based on Sobacchi & Messinger (2013) eqn 7
  // with f_* removed and f_b added since we define f_coll as M_*/M_tot rather than M_vir/M_tot,
  // and also with the inclusion of the effects of the Helium fraction.
  params->ReionEfficiency =
    1.0 / run_globals.params.BaryonFrac * params->ReionNionPhotPerBary / (1.0 - 0.75 * params->Y_He);

  // Account for instantaneous recycling factor so that stellar mass is cumulative
  if (params->Flag_IRA)
    params->ReionEfficiency /= params->SfRecycleFraction;

  mlog("Set value of run_globals.params.ReionEfficiency = %g", MLOG_MESG, params->ReionEfficiency);
}

void assign_slabs()
{
  mlog("Assigning slabs to MPI cores...", MLOG_OPEN);

  // Assign the slab size
  int n_rank = run_globals.mpi_size;
  int dim = run_globals.params.ReionGridDim;

  // Use fftw to find out what slab each rank should get
  ptrdiff_t local_nix, local_ix_start;
  ptrdiff_t local_n_complex =
    fftwf_mpi_local_size_3d(dim, dim, dim / 2 + 1, run_globals.mpi_comm, &local_nix, &local_ix_start);

  // let every core know...
  ptrdiff_t** slab_nix = &run_globals.reion_grids.slab_nix;
  *slab_nix = malloc(sizeof(ptrdiff_t) * n_rank); ///< array of number of x cells of every rank
  MPI_Allgather(&local_nix, sizeof(ptrdiff_t), MPI_BYTE, *slab_nix, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm);

  ptrdiff_t** slab_ix_start = &run_globals.reion_grids.slab_ix_start;
  *slab_ix_start = malloc(sizeof(ptrdiff_t) * n_rank); ///< array first x cell of every rank
  (*slab_ix_start)[0] = 0;
  for (int ii = 1; ii < n_rank; ii++)
    (*slab_ix_start)[ii] = (*slab_ix_start)[ii - 1] + (*slab_nix)[ii - 1];

  ptrdiff_t** slab_n_complex = &run_globals.reion_grids.slab_n_complex; ///< array of allocation counts for every rank
  *slab_n_complex = malloc(sizeof(ptrdiff_t) * n_rank);                 ///< array of allocation counts for every rank
  MPI_Allgather(
    &local_n_complex, sizeof(ptrdiff_t), MPI_BYTE, *slab_n_complex, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm);

  mlog("...done", MLOG_CLOSE);
}

void call_find_HII_bubbles(int snapshot, int nout_gals, timer_info* timer)
{
  // Thin wrapper round find_HII_bubbles

  int total_n_out_gals = 0;

  reion_grids_t* grids = &(run_globals.reion_grids);

  mlog("Getting ready to call find_HII_bubbles...", MLOG_OPEN);

  // Check to see if there are actually any galaxies at this snapshot
  MPI_Allreduce(&nout_gals, &total_n_out_gals, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
  //    if (total_n_out_gals == 0) {
  //        mlog("No galaxies in the simulation - skipping...", MLOG_CLOSE);
  //        return;
  //    }

  // Logic statement to avoid gridding the density field twice
  if (!run_globals.params.Flag_IncludeSpinTemp) {

    // Construct the baryon grids
    construct_baryon_grids(snapshot, nout_gals);

    // Read in the dark matter density grid
    read_grid(DENSITY, snapshot, grids->deltax);

    // save the grids prior to doing FFTs to avoid precision loss and aliasing etc.
    for (int i_out = 0; i_out < run_globals.NOutputSnaps; i_out++)
      if (snapshot == run_globals.ListOutputSnaps[i_out] && run_globals.params.Flag_OutputGrids)
        save_reion_input_grids(snapshot);
  }

  mlog("...done", MLOG_CLOSE);

  // Call find_HII_bubbles
  mlog("Calling find_HII_bubbles", MLOG_OPEN | MLOG_TIMERSTART);

  // Call find_HII_bubbles
  find_HII_bubbles(snapshot, timer);

  mlog("grids->volume_weighted_global_xH = %g", MLOG_MESG, grids->volume_weighted_global_xH);
  mlog("grids->volume_weighted_global_J_21 = %g", MLOG_MESG, grids->volume_weighted_global_J_21);
  mlog("global mass weighted xHII = %g at z = %g",
       MLOG_MESG,
       1.0 - grids->mass_weighted_global_xH,
       run_globals.ZZ[snapshot]);
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

void call_ComputeTs(int snapshot, int nout_gals, timer_info* timer)
{
  // Thin wrapper round ComputeTs

  int total_n_out_gals = 0;

  reion_grids_t* grids = &(run_globals.reion_grids);

  mlog("Getting ready to call ComputeTs...", MLOG_OPEN);

  // Check to see if there are actually any galaxies at this snapshot
  MPI_Allreduce(&nout_gals, &total_n_out_gals, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);

  // Construct the baryon grids
  construct_baryon_grids(snapshot, nout_gals);

  // Read in the dark matter density grid
  read_grid(DENSITY, snapshot, grids->deltax);

  // read in the velocity grids (only works for GBPTREES_TREES at the moment)
  if (run_globals.params.Flag_IncludePecVelsFor21cm > 0) {
    read_grid(run_globals.params.TsVelocityComponent, snapshot, grids->vel);
  }

  // save the grids prior to doing FFTs to avoid precision loss and aliasing etc.
  for (int i_out = 0; i_out < run_globals.NOutputSnaps; i_out++)
    if (snapshot == run_globals.ListOutputSnaps[i_out] && run_globals.params.Flag_OutputGrids &&
        !run_globals.params.FlagMCMC)
      save_reion_input_grids(snapshot);

  mlog("...done", MLOG_CLOSE);

  // Call Compute Ts
  mlog("Calling ComputeTs", MLOG_OPEN | MLOG_TIMERSTART);

  ComputeTs(snapshot, timer);
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

void init_reion_grids()
{

  reion_grids_t* grids = &(run_globals.reion_grids);
  int ReionGridDim = run_globals.params.ReionGridDim;
  ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
  ptrdiff_t slab_n_real = slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim;
  ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

  ptrdiff_t slab_n_real_smoothedSFR;
  if (run_globals.params.Flag_IncludeSpinTemp) {
    slab_n_real_smoothedSFR =
      slab_nix[run_globals.mpi_rank] * run_globals.params.TsNumFilterSteps * ReionGridDim * ReionGridDim;
  }

  ptrdiff_t slab_n_real_LC;
  if (run_globals.params.Flag_ConstructLightcone) {
    slab_n_real_LC = slab_nix[run_globals.mpi_rank] * ReionGridDim * run_globals.params.LightconeLength;
  }

  mlog("Initialising grids...", MLOG_MESG);

  grids->volume_weighted_global_xH = 1.0;
  grids->volume_weighted_global_J_21 = 0.0;
  grids->mass_weighted_global_xH = 1.0;
  grids->started = 0;
  grids->finished = 0;

  grids->volume_ave_J_alpha = 0.0;
  grids->volume_ave_xalpha = 0.0;
  grids->volume_ave_Xheat = 0.0;
  grids->volume_ave_Xion = 0.0;
  grids->volume_ave_TS = 0.0;
  grids->volume_ave_TK = 0.0;
  grids->volume_ave_xe = 0.0;
  grids->volume_ave_Tb = 0.0;

  for (int ii = 0; ii < slab_n_real; ii++) {
    grids->xH[ii] = 1.0;
    grids->z_at_ionization[ii] = -1;
    grids->r_bubble[ii] = 0.0;
    if (run_globals.params.Flag_IncludeSpinTemp) {
      grids->Tk_box[ii] = 0.0;
      grids->TS_box[ii] = 0.0;
    }
    if (run_globals.params.Flag_IncludeRecombinations) {
      grids->z_re[ii] = 0.0;
      grids->Gamma12[ii] = 0.0;
    }
    if (run_globals.params.Flag_Compute21cmBrightTemp) {
      grids->delta_T[ii] = 0.0;
      if (run_globals.params.Flag_ConstructLightcone) {
        grids->delta_T_prev[ii] = 0.0;
      }
    }
  }

  if (run_globals.params.Flag_IncludeSpinTemp) {

    for (int ii = 0; ii < slab_n_real_smoothedSFR; ii++) {
      grids->SMOOTHED_SFR_GAL[ii] = 0.0;

      if (run_globals.params.Flag_SeparateQSOXrays) {
        grids->SMOOTHED_SFR_QSO[ii] = 0.0;
      }
    }
  }

  if (run_globals.params.Flag_ConstructLightcone) {
    for (int ii = 0; ii < slab_n_real_LC; ii++) {
      grids->LightconeBox[ii] = 0.0;
    }

    for (int ii = 0; ii < run_globals.params.LightconeLength; ii++) {
      grids->Lightcone_redshifts[ii] = 0.0;
    }
  }

  for (int ii = 0; ii < slab_n_real; ii++)
    if (run_globals.params.ReionUVBFlag) {
      grids->J_21_at_ionization[ii] = (float)0.;
      grids->J_21[ii] = (float)0.;
      grids->Mvir_crit[ii] = (float)0.;
    }

  for (int ii = 0; ii < slab_n_complex; ii++) {
    grids->stars_filtered[ii] = 0 + 0I;
    grids->stars_unfiltered[ii] = 0 + 0I;
    grids->deltax_filtered[ii] = 0 + 0I;
    grids->deltax_unfiltered[ii] = 0 + 0I;
    grids->sfr_filtered[ii] = 0 + 0I;
    grids->sfr_unfiltered[ii] = 0 + 0I;
    if (run_globals.params.Flag_IncludeSpinTemp) {
      grids->x_e_filtered[ii] = 0 + 0I;
      grids->x_e_unfiltered[ii] = 0 + 0I;
    }
    if (run_globals.params.Flag_IncludeRecombinations) {
      grids->N_rec_filtered[ii] = 0 + 0I;
      grids->N_rec_unfiltered[ii] = 0 + 0I;
    }
    if (run_globals.params.Flag_Compute21cmBrightTemp && (run_globals.params.Flag_IncludePecVelsFor21cm > 0)) {
      grids->vel_gradient[ii] = 0 + 0I;
    }
  }

  for (int ii = 0; ii < slab_n_complex * 2; ii++) {
    grids->deltax[ii] = 0;
    grids->stars[ii] = 0;
    grids->sfr[ii] = 0;

    if (run_globals.params.Flag_IncludeSpinTemp) {
      grids->x_e_box_prev[ii] = 0;
      grids->x_e_box[ii] = 0;
    }
    if (run_globals.params.Flag_IncludeRecombinations) {
      grids->N_rec[ii] = 0;
    }
    if (run_globals.params.Flag_Compute21cmBrightTemp && (run_globals.params.Flag_IncludePecVelsFor21cm > 0)) {
      grids->vel[ii] = 0;
    }
  }

  if (run_globals.params.Flag_ComputePS) {
    for (int ii = 0; ii < run_globals.params.PS_Length; ii++) {
      grids->PS_k[ii] = (float)0.;
      grids->PS_data[ii] = (float)0.;
      grids->PS_error[ii] = (float)0.;
    }
  }
}

void malloc_reionization_grids()
{
  reion_grids_t* grids = &(run_globals.reion_grids);

  fftwf_mpi_init();

  // Load wisdom if requested
  run_globals.reion_grids.flag_wisdom = strlen(run_globals.params.FFTW3WisdomDir) > 0;
  bool save_wisdom = false;
  char wisdom_fname[STRLEN + 32];
  const int ReionGridDim = run_globals.params.ReionGridDim;
  unsigned plan_flags = FFTW_ESTIMATE;

  if (run_globals.reion_grids.flag_wisdom) {
    plan_flags = FFTW_PATIENT;
    sprintf(wisdom_fname,
            "%s/fftw3f-meraxes-N_%d-ranks_%d.wisdom",
            run_globals.params.FFTW3WisdomDir,
            ReionGridDim,
            run_globals.mpi_size);
    if (run_globals.mpi_rank == 0) {
      if (fftwf_import_wisdom_from_filename(wisdom_fname)) {
        mlog("Successfully loaded FFTW3 wisdom from %s", MLOG_MESG, wisdom_fname);
      } else {
        mlog("FFTW3 wisdom directory provided, but no suitable wisdom exists. New wisdom will be created (this make "
             "take a while).",
             MLOG_MESG);
        save_wisdom = true;
        // Check to see if the wisdom directory exists and if not, create it
        struct stat filestatus;
        if (stat(run_globals.params.FFTW3WisdomDir, &filestatus) != 0)
          mkdir(run_globals.params.FFTW3WisdomDir, 02755);
      }
    }
    MPI_Bcast(&save_wisdom, 1, MPI_C_BOOL, 0, run_globals.mpi_comm);
    if (!save_wisdom) {
      fftwf_mpi_broadcast_wisdom(run_globals.mpi_comm);
    }
  }

  // run_globals.NStoreSnapshots is set in `initialize_halo_storage`
  run_globals.SnapshotDeltax = (float**)calloc((size_t)run_globals.NStoreSnapshots, sizeof(float*));
  run_globals.SnapshotVel = (float**)calloc((size_t)run_globals.NStoreSnapshots, sizeof(float*));

  grids->galaxy_to_slab_map = NULL;

  grids->xH = NULL;
  grids->stars = NULL;
  grids->stars_unfiltered = NULL;
  grids->stars_filtered = NULL;
  grids->deltax = NULL;
  grids->deltax_unfiltered = NULL;
  grids->deltax_filtered = NULL;
  grids->sfr = NULL;
  grids->sfr_unfiltered = NULL;
  grids->sfr_filtered = NULL;
  grids->z_at_ionization = NULL;
  grids->J_21_at_ionization = NULL;
  grids->J_21 = NULL;

  // Grids required for the spin temperature calculation
  grids->x_e_box = NULL;
  grids->x_e_box_prev = NULL;
  grids->Tk_box = NULL;
  grids->TS_box = NULL;
  grids->x_e_unfiltered = NULL;
  grids->x_e_filtered = NULL;

  grids->SMOOTHED_SFR_GAL = NULL;
  grids->SMOOTHED_SFR_QSO = NULL;

  // Grids required for inhomogeneous recombinations
  grids->N_rec_unfiltered = NULL;
  grids->N_rec_filtered = NULL;
  grids->z_re = NULL;
  grids->Gamma12 = NULL;
  grids->N_rec = NULL;
  grids->N_rec_filtered = NULL;
  grids->N_rec_unfiltered = NULL;

  // Grids required for 21cm brightness temperature
  grids->delta_T = NULL;
  grids->delta_T_prev = NULL;

  // A grid for the lightcone (cuboid) box
  grids->LightconeBox = NULL;
  grids->Lightcone_redshifts = NULL;

  // Grids required for addining in peculiar velocity effects
  grids->vel = NULL;
  grids->vel_gradient = NULL;

  grids->PS_k = NULL;
  grids->PS_data = NULL;
  grids->PS_error = NULL;

  if (run_globals.params.Flag_PatchyReion) {
    assign_slabs();

    int ReionGridDim = run_globals.params.ReionGridDim;
    ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
    ptrdiff_t slab_n_real = slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim;
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

    ptrdiff_t slab_n_real_smoothedSFR;
    if (run_globals.params.Flag_IncludeSpinTemp) {
      slab_n_real_smoothedSFR =
        slab_nix[run_globals.mpi_rank] * run_globals.params.TsNumFilterSteps * ReionGridDim * ReionGridDim;
    }

    ptrdiff_t slab_n_real_LC;
    if (run_globals.params.Flag_ConstructLightcone) {
      slab_n_real_LC = slab_nix[run_globals.mpi_rank] * ReionGridDim * run_globals.params.LightconeLength;
    }

    // create a buffer on each rank which is as large as the largest LOGICAL allocation on any single rank
    int max_cells = 0;

    for (int ii = 0; ii < run_globals.mpi_size; ii++)
      if (slab_nix[ii] > max_cells)
        max_cells = (int)slab_nix[ii];

    max_cells *= ReionGridDim * ReionGridDim;
    grids->buffer_size = max_cells;

    grids->buffer = fftwf_alloc_real((size_t)max_cells);

    grids->stars = fftwf_alloc_real((size_t)slab_n_complex * 2);
    grids->stars_unfiltered = fftwf_alloc_complex((size_t)slab_n_complex);
    grids->stars_filtered = fftwf_alloc_complex((size_t)slab_n_complex);

    grids->stars_forward_plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim,
                                                          ReionGridDim,
                                                          ReionGridDim,
                                                          grids->stars,
                                                          grids->stars_unfiltered,
                                                          run_globals.mpi_comm,
                                                          plan_flags);
    grids->stars_filtered_reverse_plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                                   ReionGridDim,
                                                                   ReionGridDim,
                                                                   grids->stars_filtered,
                                                                   (float*)grids->stars_filtered,
                                                                   run_globals.mpi_comm,
                                                                   plan_flags);

    grids->deltax = fftwf_alloc_real((size_t)slab_n_complex * 2);
    grids->deltax_unfiltered = fftwf_alloc_complex((size_t)slab_n_complex);
    grids->deltax_filtered = fftwf_alloc_complex((size_t)slab_n_complex);

    grids->deltax_forward_plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim,
                                                           ReionGridDim,
                                                           ReionGridDim,
                                                           grids->deltax,
                                                           grids->deltax_unfiltered,
                                                           run_globals.mpi_comm,
                                                           plan_flags);
    grids->deltax_filtered_reverse_plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                                    ReionGridDim,
                                                                    ReionGridDim,
                                                                    grids->deltax_filtered,
                                                                    (float*)grids->deltax_filtered,
                                                                    run_globals.mpi_comm,
                                                                    plan_flags);

    grids->sfr = fftwf_alloc_real((size_t)slab_n_complex * 2);
    grids->sfr_unfiltered = fftwf_alloc_complex((size_t)slab_n_complex);
    grids->sfr_filtered = fftwf_alloc_complex((size_t)slab_n_complex);

    grids->sfr_forward_plan = fftwf_mpi_plan_dft_r2c_3d(
      ReionGridDim, ReionGridDim, ReionGridDim, grids->sfr, grids->sfr_unfiltered, run_globals.mpi_comm, plan_flags);
    grids->sfr_filtered_reverse_plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                                 ReionGridDim,
                                                                 ReionGridDim,
                                                                 grids->sfr_filtered,
                                                                 (float*)grids->sfr_filtered,
                                                                 run_globals.mpi_comm,
                                                                 plan_flags);

    grids->xH = fftwf_alloc_real((size_t)slab_n_real);
    grids->z_at_ionization = fftwf_alloc_real((size_t)slab_n_real);
    grids->r_bubble = fftwf_alloc_real((size_t)slab_n_real);

    if (run_globals.params.Flag_IncludeSpinTemp) {

      grids->x_e_box = fftwf_alloc_real((size_t)slab_n_complex * 2);
      grids->x_e_unfiltered = fftwf_alloc_complex((size_t)slab_n_complex);
      grids->x_e_filtered = fftwf_alloc_complex((size_t)slab_n_complex);
      grids->x_e_box_prev = fftwf_alloc_real((size_t)slab_n_complex * 2);

      grids->x_e_box_forward_plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim,
                                                              ReionGridDim,
                                                              ReionGridDim,
                                                              grids->x_e_box,
                                                              grids->x_e_unfiltered,
                                                              run_globals.mpi_comm,
                                                              plan_flags);
      grids->x_e_filtered_reverse_plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                                   ReionGridDim,
                                                                   ReionGridDim,
                                                                   grids->x_e_filtered,
                                                                   (float*)grids->x_e_filtered,
                                                                   run_globals.mpi_comm,
                                                                   plan_flags);

      grids->Tk_box = fftwf_alloc_real((size_t)slab_n_real);
      grids->TS_box = fftwf_alloc_real((size_t)slab_n_real);

      grids->SMOOTHED_SFR_GAL = calloc((size_t)slab_n_real_smoothedSFR, sizeof(double));
      if (run_globals.params.Flag_SeparateQSOXrays) {
        grids->SMOOTHED_SFR_QSO = calloc((size_t)slab_n_real_smoothedSFR, sizeof(double));
      }
    }

    if (run_globals.params.Flag_IncludeRecombinations) {
      grids->N_rec = fftwf_alloc_real((size_t)slab_n_complex * 2);
      grids->N_rec_unfiltered = fftwf_alloc_complex((size_t)slab_n_complex);
      grids->N_rec_filtered = fftwf_alloc_complex((size_t)slab_n_complex);

      grids->N_rec_forward_plan = fftwf_mpi_plan_dft_r2c_3d(ReionGridDim,
                                                            ReionGridDim,
                                                            ReionGridDim,
                                                            grids->N_rec,
                                                            grids->N_rec_unfiltered,
                                                            run_globals.mpi_comm,
                                                            plan_flags);
      grids->N_rec_filtered_reverse_plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                                     ReionGridDim,
                                                                     ReionGridDim,
                                                                     grids->N_rec_filtered,
                                                                     (float*)grids->N_rec_filtered,
                                                                     run_globals.mpi_comm,
                                                                     plan_flags);

      grids->z_re = fftwf_alloc_real((size_t)slab_n_real);
      grids->Gamma12 = fftwf_alloc_real((size_t)slab_n_real);
    }

    if (run_globals.params.Flag_Compute21cmBrightTemp) {
      grids->delta_T = fftwf_alloc_real((size_t)slab_n_real);

      if (run_globals.params.Flag_IncludePecVelsFor21cm > 0) {
        grids->vel = fftwf_alloc_real((size_t)slab_n_complex * 2);
        grids->vel_gradient = fftwf_alloc_complex((size_t)slab_n_complex);

        grids->vel_forward_plan = fftwf_mpi_plan_dft_r2c_3d(
          ReionGridDim, ReionGridDim, ReionGridDim, grids->vel, grids->vel_gradient, run_globals.mpi_comm, plan_flags);
        grids->vel_gradient_reverse_plan = fftwf_mpi_plan_dft_c2r_3d(ReionGridDim,
                                                                     ReionGridDim,
                                                                     ReionGridDim,
                                                                     grids->vel_gradient,
                                                                     (float*)grids->vel_gradient,
                                                                     run_globals.mpi_comm,
                                                                     plan_flags);
      }

      if (run_globals.params.Flag_ConstructLightcone) {
        grids->delta_T_prev = fftwf_alloc_real((size_t)slab_n_real);
      }
    }

    if (run_globals.params.ReionUVBFlag) {
      grids->J_21_at_ionization = fftwf_alloc_real((size_t)slab_n_real);
      grids->J_21 = fftwf_alloc_real((size_t)slab_n_real);
      grids->Mvir_crit = fftwf_alloc_real((size_t)slab_n_real);
    }

    if (run_globals.params.Flag_ConstructLightcone) {
      grids->LightconeBox = fftwf_alloc_real((size_t)slab_n_real_LC);
      grids->Lightcone_redshifts = fftwf_alloc_real((size_t)run_globals.params.LightconeLength);
    }

    if (run_globals.params.Flag_ComputePS) {
      grids->PS_k = fftwf_alloc_real((size_t)run_globals.params.PS_Length);
      grids->PS_data = fftwf_alloc_real((size_t)run_globals.params.PS_Length);
      grids->PS_error = fftwf_alloc_real((size_t)run_globals.params.PS_Length);
    }

    init_reion_grids();

    if (run_globals.reion_grids.flag_wisdom && save_wisdom) {
      fftwf_mpi_gather_wisdom(run_globals.mpi_comm);
      if (run_globals.mpi_rank == 0) {
        if (fftwf_export_wisdom_to_filename(wisdom_fname)) {
          mlog("Successfully saved FFTW3 wisdom to %s", MLOG_MESG, wisdom_fname);
        }
      }
    }

  } // if (run_globals.params.Flag_PatchyReion)
}

void free_reionization_grids()
{
  mlog("Freeing reionization grids...", MLOG_OPEN);

  reion_grids_t* grids = &(run_globals.reion_grids);

  free(run_globals.reion_grids.slab_n_complex);
  free(run_globals.reion_grids.slab_ix_start);
  free(run_globals.reion_grids.slab_nix);

  if (run_globals.params.Flag_ComputePS) {
    fftwf_free(grids->PS_error);
    fftwf_free(grids->PS_data);
    fftwf_free(grids->PS_k);
  }

  if (run_globals.params.Flag_ConstructLightcone) {
    fftwf_free(grids->LightconeBox);
  }

  if (run_globals.params.ReionUVBFlag) {
    fftwf_free(grids->Mvir_crit);
    fftwf_free(grids->J_21);
    fftwf_free(grids->J_21_at_ionization);
  }

  if (run_globals.params.Flag_Compute21cmBrightTemp) {

    if (run_globals.params.Flag_ConstructLightcone) {
      fftwf_free(grids->delta_T_prev);
      fftwf_free(grids->Lightcone_redshifts);
    }

    if (run_globals.params.Flag_IncludePecVelsFor21cm > 0) {
      fftwf_destroy_plan(grids->vel_gradient_reverse_plan);
      fftwf_destroy_plan(grids->vel_forward_plan);
      fftwf_free(grids->vel_gradient);
      fftwf_free(grids->vel);
    }

    fftwf_free(grids->delta_T);
  }

  if (run_globals.params.Flag_IncludeRecombinations) {

    fftwf_free(grids->Gamma12);
    fftwf_free(grids->z_re);

    fftwf_destroy_plan(grids->N_rec_filtered_reverse_plan);
    fftwf_destroy_plan(grids->N_rec_forward_plan);
    fftwf_free(grids->N_rec_filtered);
    fftwf_free(grids->N_rec_unfiltered);
    fftwf_free(grids->N_rec);
  }

  if (run_globals.params.Flag_IncludeSpinTemp) {
    free(grids->SMOOTHED_SFR_QSO);
    free(grids->SMOOTHED_SFR_GAL);

    fftwf_free(grids->Tk_box);
    fftwf_free(grids->TS_box);

    fftwf_destroy_plan(grids->x_e_filtered_reverse_plan);
    fftwf_destroy_plan(grids->x_e_box_forward_plan);
    fftwf_free(grids->x_e_filtered);
    fftwf_free(grids->x_e_unfiltered);
    fftwf_free(grids->x_e_box_prev);
    fftwf_free(grids->x_e_box);
  }

  fftwf_free(grids->r_bubble);
  fftwf_free(grids->z_at_ionization);
  fftwf_free(grids->xH);

  fftwf_destroy_plan(grids->sfr_filtered_reverse_plan);
  fftwf_destroy_plan(grids->sfr_forward_plan);
  fftwf_free(grids->sfr_filtered);
  fftwf_free(grids->sfr_unfiltered);
  fftwf_free(grids->sfr);

  fftwf_destroy_plan(grids->deltax_filtered_reverse_plan);
  fftwf_destroy_plan(grids->deltax_forward_plan);
  fftwf_free(grids->deltax_filtered);
  fftwf_free(grids->deltax_unfiltered);
  fftwf_free(grids->deltax);

  fftwf_destroy_plan(grids->stars_filtered_reverse_plan);
  fftwf_destroy_plan(grids->stars_forward_plan);
  fftwf_free(grids->stars_filtered);
  fftwf_free(grids->stars_unfiltered);
  fftwf_free(grids->stars);

  fftwf_free(grids->buffer);

  mlog(" ...done", MLOG_CLOSE);
}

int map_galaxies_to_slabs(int ngals)
{
  double box_size = run_globals.params.BoxSize;
  int ReionGridDim = run_globals.params.ReionGridDim;

  mlog("Mapping galaxies to slabs...", MLOG_OPEN);

  // Loop through each valid galaxy and find what slab it sits in
  if (ngals > 0)
    run_globals.reion_grids.galaxy_to_slab_map = malloc(sizeof(gal_to_slab_t) * ngals);
  else
    run_globals.reion_grids.galaxy_to_slab_map = NULL;

  gal_to_slab_t* galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
  ptrdiff_t* slab_ix_start = run_globals.reion_grids.slab_ix_start;

  galaxy_t* gal = run_globals.FirstGal;
  int gal_counter = 0;
  while (gal != NULL) {
    // TODO: Note that I am including ghosts here.  We will need to check the
    // validity of this.  By definition, if they are ghosts then their host
    // halo hasn't been identified at this time step and hence they haven't
    // been evolved.  Their properties (Sfr, StellarMass, etc.) will all have
    // been set when they were last identified.
    if (gal->Type < 3) {
      // TODO: for type 2 galaxies these positions will be set from the last
      // time they were identified.  If their host halo has moved significantly
      // since then, these positions won't reflect that and the satellites will
      // be spatially disconnected from their hosts.  We will need to fix this
      // at some point.

      ptrdiff_t ix = pos_to_ngp(gal->Pos[0], box_size, ReionGridDim);

      assert((ix >= 0) && (ix < ReionGridDim));

      galaxy_to_slab_map[gal_counter].index = gal_counter;
      galaxy_to_slab_map[gal_counter].slab_ind =
        searchsorted(&ix, slab_ix_start, run_globals.mpi_size, sizeof(ptrdiff_t), compare_ptrdiff, -1, -1);
      galaxy_to_slab_map[gal_counter++].galaxy = gal;
    }

    gal = gal->Next;
  }

  // sort the slab indices IN PLACE (n.b. compare_slab_assign is a stable comparison)
  if (galaxy_to_slab_map != NULL)
    qsort(galaxy_to_slab_map, (size_t)gal_counter, sizeof(gal_to_slab_t), compare_slab_assign);

  assert(gal_counter == ngals);

  mlog("...done.", MLOG_CLOSE);

  return gal_counter;
}

void assign_Mvir_crit_to_galaxies(int ngals_in_slabs)
{
  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  gal_to_slab_t* galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
  float* Mvir_crit = run_globals.reion_grids.Mvir_crit;
  float* buffer = run_globals.reion_grids.buffer;
  ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
  ptrdiff_t* slab_ix_start = run_globals.reion_grids.slab_ix_start;
  int ReionGridDim = run_globals.params.ReionGridDim;
  double box_size = run_globals.params.BoxSize;
  int total_assigned = 0;

  mlog("Assigning Mvir_crit to galaxies...", MLOG_OPEN);

  // Work out the index of the galaxy_to_slab_map where each slab begins.
  // TODO: This needs checked...
  int slab_map_offsets[run_globals.mpi_size];
  for (int ii = 0, i_gal = 0; ii < run_globals.mpi_size; ii++) {
    if (galaxy_to_slab_map != NULL) {
      while ((i_gal < (ngals_in_slabs - 1)) && (galaxy_to_slab_map[i_gal].slab_ind < ii))
        i_gal++;

      if (galaxy_to_slab_map[i_gal].slab_ind == ii)
        slab_map_offsets[ii] = i_gal;
      else
        slab_map_offsets[ii] = -1;
    } else
      // if this core has no galaxies then the offsets are -1 everywhere
      slab_map_offsets[ii] = -1;
  }

  // DEBUG
  // for (int ii = 0; ii < run_globals.mpi_size; ii++) {
  //     if (run_globals.mpi_rank == ii) {
  //         mlog("slab_map_offsets[%d] = [ ", MLOG_MESG|MLOG_ALLRANKS, ii);
  //         for (int jj = 0; jj < run_globals.mpi_size; ++jj)
  //             printf("%d ", slab_map_offsets[jj]);
  //         printf("]\n");
  //     }
  //     MPI_Barrier(run_globals.mpi_comm);
  // }

  // do a ring exchange of slabs between all cores
  for (int i_skip = 0; i_skip < run_globals.mpi_size; i_skip++) {
    int recv_from_rank = (run_globals.mpi_rank + i_skip) % run_globals.mpi_size;
    int send_to_rank = (run_globals.mpi_rank - i_skip + run_globals.mpi_size) % run_globals.mpi_size;

    bool send_flag = false;
    bool recv_flag = (slab_map_offsets[recv_from_rank] > -1);

    if (i_skip > 0) {
      MPI_Sendrecv(&recv_flag,
                   sizeof(bool),
                   MPI_BYTE,
                   recv_from_rank,
                   6393762,
                   &send_flag,
                   sizeof(bool),
                   MPI_BYTE,
                   send_to_rank,
                   6393762,
                   run_globals.mpi_comm,
                   MPI_STATUS_IGNORE);

      // need to ensure sends and receives do not clash!
      if (send_to_rank > run_globals.mpi_rank) {
        if (send_flag) {
          int n_cells = (int)(slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim);
          MPI_Send(Mvir_crit, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
        }
        if (recv_flag) {
          int n_cells = (int)(slab_nix[recv_from_rank] * ReionGridDim * ReionGridDim);
          MPI_Recv(buffer, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
        }
      } else {
        if (recv_flag) {
          int n_cells = (int)(slab_nix[recv_from_rank] * ReionGridDim * ReionGridDim);
          MPI_Recv(buffer, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
        }
        if (send_flag) {
          int n_cells = (int)(slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim);
          MPI_Send(Mvir_crit, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
        }
      }
    } else {
      int n_cells = (int)(slab_nix[recv_from_rank] * ReionGridDim * ReionGridDim);
      memcpy(buffer, Mvir_crit, sizeof(float) * n_cells);
    }

    // if this core has received a slab of Mvir_crit then assign values to the
    // galaxies which belong to this slab
    if (recv_flag) {
      int i_gal = slab_map_offsets[recv_from_rank];
      int ix_start = (int)slab_ix_start[recv_from_rank];
      while ((i_gal < ngals_in_slabs) && (galaxy_to_slab_map[i_gal].slab_ind == recv_from_rank)) {
        // TODO: We should use the position of the FOF group here...
        galaxy_t* gal = galaxy_to_slab_map[i_gal].galaxy;
        int ix = pos_to_ngp(gal->Pos[0], box_size, ReionGridDim) - ix_start;
        int iy = pos_to_ngp(gal->Pos[1], box_size, ReionGridDim);
        int iz = pos_to_ngp(gal->Pos[2], box_size, ReionGridDim);

        assert(ix >= 0);
        assert(ix < slab_nix[recv_from_rank]);

        // Record the Mvir_crit (filtering mass) value
        gal->MvirCrit = (double)buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];

        // increment counters
        i_gal++;
        total_assigned++;
      }
    }
  }

  if (total_assigned != ngals_in_slabs)
    ABORT(EXIT_FAILURE);

  mlog("...done.", MLOG_CLOSE);
}

void construct_baryon_grids(int snapshot, int local_ngals)
{
  double box_size = run_globals.params.BoxSize;
  float* stellar_grid = run_globals.reion_grids.stars;
  float* sfr_grid = run_globals.reion_grids.sfr;
  int ReionGridDim = run_globals.params.ReionGridDim;
  double sfr_timescale = run_globals.params.ReionSfrTimescale * hubble_time(snapshot);

  gal_to_slab_t* galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
  ptrdiff_t* slab_ix_start = run_globals.reion_grids.slab_ix_start;
  int local_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);

  mlog("Constructing stellar mass and sfr grids...", MLOG_OPEN | MLOG_TIMERSTART);

  // init the grid
  for (int ii = 0; ii < local_n_complex * 2; ii++) {
    stellar_grid[ii] = 0.0;
    sfr_grid[ii] = 0.0;
  }

  // loop through each slab
  //
  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
  ptrdiff_t buffer_size = run_globals.reion_grids.buffer_size;
  float* buffer = run_globals.reion_grids.buffer;

  enum property
  {
    prop_stellar,
    prop_sfr
  };
  for (int prop = prop_stellar; prop <= prop_sfr; prop++) {
    int i_gal = 0;
    int skipped_gals = 0;
    long N_BlackHoleMassLimitReion = 0;

    for (int i_r = 0; i_r < run_globals.mpi_size; i_r++) {
      // init the buffer
      for (int ii = 0; ii < buffer_size; ii++)
        buffer[ii] = (float)0.;

      // if this core holds no galaxies then we don't need to fill the buffer
      if (local_ngals != 0)
        // fill the local buffer for this slab
        while (((i_gal - skipped_gals) < local_ngals) && (galaxy_to_slab_map[i_gal].slab_ind == i_r)) {
          galaxy_t* gal = galaxy_to_slab_map[i_gal].galaxy;

          // Dead galaxies should not be included here and are not in the
          // local_ngals count.  They will, however, have been assigned to a
          // slab so we will need to ignore them here...
          if (gal->Type > 2) {
            i_gal++;
            skipped_gals++;
            continue;
          }

          assert(galaxy_to_slab_map[i_gal].index >= 0);
          assert((galaxy_to_slab_map[i_gal].slab_ind >= 0) &&
                 (galaxy_to_slab_map[i_gal].slab_ind < run_globals.mpi_size));

          int ix = (int)(pos_to_ngp(gal->Pos[0], box_size, ReionGridDim) - slab_ix_start[i_r]);
          int iy = pos_to_ngp(gal->Pos[1], box_size, ReionGridDim);
          int iz = pos_to_ngp(gal->Pos[2], box_size, ReionGridDim);

          assert((ix < slab_nix[i_r]) && (ix >= 0));
          assert((iy < ReionGridDim) && (iy >= 0));
          assert((iz < ReionGridDim) && (iz >= 0));

          int ind = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

          assert((ind >= 0) && (ind < slab_nix[i_r] * ReionGridDim * ReionGridDim));

          // They are the same just now, but may be different in the future once the model is improved.
          switch (prop) {
            case prop_stellar:

              buffer[ind] += gal->FescWeightedGSM;
              // a trick to include quasar radiation using current 21cmFAST code
              if (run_globals.params.physics.Flag_BHFeedback) {
                if (gal->BlackHoleMass >= run_globals.params.physics.BlackHoleMassLimitReion)
                  buffer[ind] += gal->EffectiveBHM;
                else
                  N_BlackHoleMassLimitReion += 1;
              }
              break;

            case prop_sfr:
              buffer[ind] += gal->FescWeightedGSM;
              // for ionizing_source_formation_rate_grid, need further convertion due to different UV spectral index of
              // quasar and stellar component
              if (run_globals.params.physics.Flag_BHFeedback)
                if (gal->BlackHoleMass >= run_globals.params.physics.BlackHoleMassLimitReion)
                  buffer[ind] += gal->EffectiveBHM * run_globals.params.physics.ReionAlphaUVBH /
                                 run_globals.params.physics.ReionAlphaUV;
              break;

            default:
              mlog_error("Unrecognised property in slab creation.");
              ABORT(EXIT_FAILURE);
              break;
          }

          i_gal++;
        }

      // reduce on to the correct rank
      if (run_globals.mpi_rank == i_r)
        MPI_Reduce(MPI_IN_PLACE, buffer, (int)buffer_size, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);
      else
        MPI_Reduce(buffer, buffer, (int)buffer_size, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);

      if (run_globals.mpi_rank == i_r)

        // Do one final pass and divide the sfr_grid by the sfr timescale
        // in order to convert the stellar masses recorded into SFRs before
        // finally copying the values into the appropriate slab.
        switch (prop) {
          case prop_sfr:
            for (int ix = 0; ix < slab_nix[i_r]; ix++)
              for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                  double val = (double)buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];
                  val = (val > 0) ? val / sfr_timescale : 0;
                  sfr_grid[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = (float)val;
                }
            break;

          case prop_stellar:
            for (int ix = 0; ix < slab_nix[i_r]; ix++)
              for (int iy = 0; iy < ReionGridDim; iy++)
                for (int iz = 0; iz < ReionGridDim; iz++) {
                  float val = buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];
                  if (val < 0)
                    val = 0;
                  stellar_grid[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = val;
                }
            break;

          default:
            mlog_error("Eh!?!");
            ABORT(EXIT_FAILURE);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &N_BlackHoleMassLimitReion, 1, MPI_LONG, MPI_SUM, run_globals.mpi_comm);
    mlog("%d quasars are smaller than %g",
         MLOG_MESG,
         N_BlackHoleMassLimitReion,
         run_globals.params.physics.BlackHoleMassLimitReion);
  }

  mlog("done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

static void write_grid_float(const char* name,
                             float* data,
                             hid_t file_id,
                             hid_t fspace_id,
                             hid_t memspace_id,
                             hid_t dcpl_id)
{
  // create the dataset
  hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, fspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

  // create the property list
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // write the dataset
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id, fspace_id, plist_id, data);

  // cleanup
  H5Pclose(plist_id);
  H5Dclose(dset_id);
}

void gen_grids_fname(const int snapshot, char* name, const bool relative)
{
  if (!relative)
    sprintf(name, "%s/%s_grids_%d.hdf5", run_globals.params.OutputDir, run_globals.params.FileNameGalaxies, snapshot);
  else
    sprintf(name, "%s_grids_%d.hdf5", run_globals.params.FileNameGalaxies, snapshot);
}

void save_reion_input_grids(int snapshot)
{
  reion_grids_t* grids = &(run_globals.reion_grids);
  int ReionGridDim = run_globals.params.ReionGridDim;
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
  double UnitTime_in_s = run_globals.units.UnitTime_in_s;
  double UnitMass_in_g = run_globals.units.UnitMass_in_g;

  mlog("Saving tocf input grids...", MLOG_OPEN);

  char name[STRLEN];
  gen_grids_fname(snapshot, name, false);

  // create the file (in parallel)
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // create the filespace
  hsize_t dims[3] = { (hsize_t)ReionGridDim, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim };
  hid_t fspace_id = H5Screate_simple(3, dims, NULL);

  // create the memspace
  hsize_t mem_dims[3] = { (hsize_t)local_nix, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim };
  hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

  // select a hyperslab in the filespace
  hsize_t start[3] = { (hsize_t)run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank], 0, 0 };
  hsize_t count[3] = { (hsize_t)local_nix, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim };
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  // set the dataset creation property list to use chunking along x-axis
  hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl_id, 3, (hsize_t[3]){ 1, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim });

  // fftw padded grids
  float* grid = (float*)calloc((size_t)local_nix * (size_t)ReionGridDim * (size_t)ReionGridDim, sizeof(float));

  for (int ii = 0; ii < local_nix; ii++)
    for (int jj = 0; jj < ReionGridDim; jj++)
      for (int kk = 0; kk < ReionGridDim; kk++)
        grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] =
          (grids->deltax)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];
  write_grid_float("deltax", grid, file_id, fspace_id, memspace_id, dcpl_id);

  for (int ii = 0; ii < local_nix; ii++)
    for (int jj = 0; jj < ReionGridDim; jj++)
      for (int kk = 0; kk < ReionGridDim; kk++)
        grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] =
          (grids->stars)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];
  write_grid_float("stars", grid, file_id, fspace_id, memspace_id, dcpl_id);

  for (int ii = 0; ii < local_nix; ii++)
    for (int jj = 0; jj < ReionGridDim; jj++)
      for (int kk = 0; kk < ReionGridDim; kk++)
        grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] =
          (float)((grids->sfr)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)] * UnitMass_in_g / UnitTime_in_s *
                  SEC_PER_YEAR / SOLAR_MASS);
  write_grid_float("sfr", grid, file_id, fspace_id, memspace_id, dcpl_id);

  // tidy up
  free(grid);
  H5Pclose(dcpl_id);
  H5Sclose(memspace_id);
  H5Sclose(fspace_id);
  H5Fclose(file_id);

  mlog("...done", MLOG_CLOSE);
}

void save_reion_output_grids(int snapshot)
{

  reion_grids_t* grids = &(run_globals.reion_grids);
  int ReionGridDim = run_globals.params.ReionGridDim;
  int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);

  // float *ps;
  // int   ps_nbins;
  // float average_deltaT;
  // double Hubble_h = run_globals.params.Hubble_h;

  // Save tocf grids
  // ----------------------------------------------------------------------------------------------------

  mlog("Saving tocf output grids...", MLOG_OPEN);

  char name[STRLEN];
  gen_grids_fname(snapshot, name, false);

  // open the file (in parallel)
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(name, H5F_ACC_RDWR, plist_id);
  H5Pclose(plist_id);

  // create the filespace
  hsize_t dims[3] = { (hsize_t)ReionGridDim, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim };
  hid_t fspace_id = H5Screate_simple(3, dims, NULL);

  // create the memspace
  hsize_t mem_dims[3] = { (hsize_t)local_nix, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim };
  hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

  // select a hyperslab in the filespace
  hsize_t start[3] = { (hsize_t)run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank], 0, 0 };
  hsize_t count[3] = { (hsize_t)local_nix, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim };
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  // set the dataset creation property list to use chunking along x-axis
  hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl_id, 3, (hsize_t[3]){ 1, (hsize_t)ReionGridDim, (hsize_t)ReionGridDim });

  // create and write the datasets
  write_grid_float("xH", grids->xH, file_id, fspace_id, memspace_id, dcpl_id);
  write_grid_float("z_at_ionization", grids->z_at_ionization, file_id, fspace_id, memspace_id, dcpl_id);
  write_grid_float("r_bubble", grids->r_bubble, file_id, fspace_id, memspace_id, dcpl_id);

  if (run_globals.params.ReionUVBFlag) {
    write_grid_float("J_21", grids->J_21, file_id, fspace_id, memspace_id, dcpl_id);
    H5LTset_attribute_double(file_id, "J_21", "volume_weighted_global_J_21", &(grids->volume_weighted_global_J_21), 1);
    write_grid_float("J_21_at_ionization", grids->J_21_at_ionization, file_id, fspace_id, memspace_id, dcpl_id);
    write_grid_float("Mvir_crit", grids->Mvir_crit, file_id, fspace_id, memspace_id, dcpl_id);
  }

  // fftw padded grids
  float* grid = (float*)calloc((size_t)(local_nix * ReionGridDim * ReionGridDim), sizeof(float));

  if (run_globals.params.Flag_IncludeSpinTemp) {
    write_grid_float("TS_box", grids->TS_box, file_id, fspace_id, memspace_id, dcpl_id);
    write_grid_float("Tk_box", grids->Tk_box, file_id, fspace_id, memspace_id, dcpl_id);

    for (int ii = 0; ii < local_nix; ii++)
      for (int jj = 0; jj < ReionGridDim; jj++)
        for (int kk = 0; kk < ReionGridDim; kk++)
          grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] =
            (grids->x_e_box_prev)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];

    write_grid_float("x_e_box", grid, file_id, fspace_id, memspace_id, dcpl_id);
  }

  if (run_globals.params.Flag_Compute21cmBrightTemp) {
    write_grid_float("delta_T", grids->delta_T, file_id, fspace_id, memspace_id, dcpl_id);
  }

  if (run_globals.params.Flag_ConstructLightcone && run_globals.params.EndSnapshotLightcone == snapshot &&
      snapshot != 0) {

    // create the filespace
    hsize_t dims_LC[3] = { (hsize_t)ReionGridDim, (hsize_t)ReionGridDim, (hsize_t)run_globals.params.LightconeLength };
    hid_t fspace_id_LC = H5Screate_simple(3, dims_LC, NULL);

    // create the memspace
    hsize_t mem_dims_LC[3] = { (hsize_t)local_nix, (hsize_t)ReionGridDim, (hsize_t)run_globals.params.LightconeLength };
    hid_t memspace_id_LC = H5Screate_simple(3, mem_dims_LC, NULL);

    // select a hyperslab in the filespace
    hsize_t start_LC[3] = { (hsize_t)run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank], 0, 0 };
    hsize_t count_LC[3] = { (hsize_t)local_nix, (hsize_t)ReionGridDim, (hsize_t)run_globals.params.LightconeLength };
    H5Sselect_hyperslab(fspace_id_LC, H5S_SELECT_SET, start_LC, NULL, count_LC, NULL);

    // set the dataset creation property list to use chunking along x-axis
    hid_t dcpl_id_LC = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id_LC, 3, (hsize_t[3]){ 1, (hsize_t)ReionGridDim, (hsize_t)run_globals.params.LightconeLength });

    mlog("Outputting light-cone", MLOG_MESG);
    write_grid_float("LightconeBox", grids->LightconeBox, file_id, fspace_id_LC, memspace_id_LC, dcpl_id_LC);

    // create the filespace
    hsize_t dims_LCz[1] = { (hsize_t)run_globals.params.LightconeLength };
    hid_t fspace_id_LCz = H5Screate_simple(1, dims_LCz, NULL);

    // create the memspace
    hsize_t mem_dims_LCz[1] = { (hsize_t)run_globals.params.LightconeLength };
    hid_t memspace_id_LCz = H5Screate_simple(1, mem_dims_LCz, NULL);

    hid_t dcpl_id_LCz = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id =
      H5Dcreate(file_id, "lightcone-z", H5T_NATIVE_FLOAT, fspace_id_LCz, H5P_DEFAULT, dcpl_id_LCz, H5P_DEFAULT);

    plist_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // write the dataset
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id_LCz, fspace_id_LCz, plist_id, grids->Lightcone_redshifts);

    // cleanup
    H5Pclose(plist_id);
    H5Dclose(dset_id);
  }

  H5LTset_attribute_double(file_id, "xH", "volume_weighted_global_xH", &(grids->volume_weighted_global_xH), 1);
  H5LTset_attribute_double(file_id, "xH", "mass_weighted_global_xH", &(grids->mass_weighted_global_xH), 1);

  if (run_globals.params.Flag_IncludeSpinTemp) {
    H5LTset_attribute_double(file_id, "TS_box", "volume_ave_TS", &(grids->volume_ave_TS), 1);
    H5LTset_attribute_double(file_id, "Tk_box", "volume_ave_TK", &(grids->volume_ave_TK), 1);
    H5LTset_attribute_double(file_id, "x_e_box", "volume_ave_xe", &(grids->volume_ave_xe), 1);

    H5LTset_attribute_double(file_id, "TS_box", "volume_ave_J_alpha", &(grids->volume_ave_J_alpha), 1);
    H5LTset_attribute_double(file_id, "TS_box", "volume_ave_xalpha", &(grids->volume_ave_xalpha), 1);
    H5LTset_attribute_double(file_id, "TS_box", "volume_ave_Xheat", &(grids->volume_ave_Xheat), 1);
    H5LTset_attribute_double(file_id, "TS_box", "volume_ave_Xion", &(grids->volume_ave_Xion), 1);
  }

  if (run_globals.params.Flag_Compute21cmBrightTemp) {
    H5LTset_attribute_double(file_id, "delta_T", "volume_ave_Tb", &(grids->volume_ave_Tb), 1);
  }

  if (run_globals.params.Flag_ComputePS) {

    // create the filespace
    hsize_t dims_PS[1] = { (hsize_t)run_globals.params.PS_Length };
    hid_t fspace_id_PS = H5Screate_simple(1, dims_PS, NULL);

    // create the memspace
    hsize_t mem_dims_PS[1] = { (hsize_t)run_globals.params.PS_Length };
    hid_t memspace_id_PS = H5Screate_simple(1, mem_dims_PS, NULL);

    hid_t dcpl_id_PS = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate(file_id, "k_bins", H5T_NATIVE_FLOAT, fspace_id_PS, H5P_DEFAULT, dcpl_id_PS, H5P_DEFAULT);

    plist_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // write the dataset
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id_PS, fspace_id_PS, plist_id, grids->PS_k);

    // cleanup
    H5Pclose(plist_id);
    H5Dclose(dset_id);

    dset_id = H5Dcreate(file_id, "PS_data", H5T_NATIVE_FLOAT, fspace_id_PS, H5P_DEFAULT, dcpl_id_PS, H5P_DEFAULT);

    plist_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id_PS, fspace_id_PS, plist_id, grids->PS_data);

    // cleanup
    H5Pclose(plist_id);
    H5Dclose(dset_id);

    dset_id = H5Dcreate(file_id, "PS_error", H5T_NATIVE_FLOAT, fspace_id_PS, H5P_DEFAULT, dcpl_id_PS, H5P_DEFAULT);

    plist_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id_PS, fspace_id_PS, plist_id, grids->PS_error);

    // cleanup
    H5Pclose(plist_id);
    H5Dclose(dset_id);
  }

  // tidy up
  free(grid);
  H5Pclose(dcpl_id);
  H5Sclose(memspace_id);
  H5Sclose(fspace_id);
  H5Fclose(file_id);

  mlog("...done", MLOG_CLOSE); // Saving tocf grids
}

bool check_if_reionization_ongoing(int snapshot)
{
  int started = run_globals.reion_grids.started;
  int finished = run_globals.reion_grids.finished;

  // First check if we've already finished on all cores.
  if (finished)
    return false;

  // Ok, so we haven't finished.  Have we started then?
  if (started) {
    // whether we want to continue even when reionization is finished
    // In order to keep outputting meraxes_grids_%d.hdf5
    if (run_globals.params.Flag_OutputGridsPostReion)
      return true;

    if (run_globals.params.Flag_ConstructLightcone && snapshot <= run_globals.params.EndSnapshotLightcone) {
      return true;
    }

    // So we have started, but have not previously found to be finished.  Have
    // we now finished though?
    float* xH = run_globals.reion_grids.xH;
    int ReionGridDim = run_globals.params.ReionGridDim;
    int slab_n_real = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]) * ReionGridDim * ReionGridDim;

    // If not all cells are ionised then reionization is still progressing...
    finished = 1;
    for (int ii = 0; ii < slab_n_real; ii++)
      if (xH[ii] != 0.0) {
        finished = 0;
        break;
      }
  } else {

    // Here we haven't finished or previously started.  Should we start then?
    if (run_globals.params.Flag_IncludeSpinTemp || run_globals.params.Flag_ConstructLightcone) {
      started = 1;
    } else {
      if (run_globals.FirstGal != NULL) {
        started = 1;
      }
    }
  }

  // At this stage, `started` and `finished` should be set accordingly for each
  // individual core.  Now we need to combine them on all cores.
  MPI_Allreduce(MPI_IN_PLACE, &started, 1, MPI_INT, MPI_LOR, run_globals.mpi_comm);
  run_globals.reion_grids.started = started;
  MPI_Allreduce(MPI_IN_PLACE, &finished, 1, MPI_INT, MPI_LAND, run_globals.mpi_comm);
  run_globals.reion_grids.finished = finished;

  if (started && (!finished))
    return true;
  else
    return false;
}

void filter(fftwf_complex* box, int local_ix_start, int slab_nx, int grid_dim, float R, int filter_type)
{
  int middle = grid_dim / 2;
  float box_size = (float)run_globals.params.BoxSize;
  float delta_k = (float)(2.0 * M_PI / box_size);

  // Loop through k-box
  for (int n_x = 0; n_x < slab_nx; n_x++) {
    float k_x;
    int n_x_global = n_x + local_ix_start;

    if (n_x_global > middle)
      k_x = (n_x_global - grid_dim) * delta_k;
    else
      k_x = n_x_global * delta_k;

    for (int n_y = 0; n_y < grid_dim; n_y++) {
      float k_y;

      if (n_y > middle)
        k_y = (n_y - grid_dim) * delta_k;
      else
        k_y = n_y * delta_k;

      for (int n_z = 0; n_z <= middle; n_z++) {
        float k_z = n_z * delta_k;

        float k_mag = sqrtf(k_x * k_x + k_y * k_y + k_z * k_z);

        float kR = k_mag * R; // Real space top-hat

        switch (filter_type) {
          case 0: // Real space top-hat
            if (kR > 1e-4)
              box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] *=
                (fftwf_complex)(3.0 * (sinf(kR) / powf(kR, 3) - cosf(kR) / powf(kR, 2)));
            break;

          case 1:              // k-space top hat
            kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
            if (kR > 1)
              box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] = (fftwf_complex)0.0;
            break;

          case 2:        // Gaussian
            kR *= 0.643; // Equates integrated volume to the real space top-hat
            box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] *=
              (fftwf_complex)(powf((float)M_E, (float)(-kR * kR / 2.0)));
            break;

          default:
            if ((n_x == 0) && (n_y == 0) && (n_z == 0)) {
              mlog_error("ReionFilterType.c: Warning, ReionFilterType type %d is undefined!", filter_type);
              ABORT(EXIT_FAILURE);
            }
            break;
        }
      }
    }
  } // End looping through k box
}

void velocity_gradient(fftwf_complex* box, int local_ix_start, int slab_nx, int grid_dim)
{
  int middle = grid_dim / 2;
  float box_size = (float)run_globals.params.BoxSize;
  float delta_k = (float)(2.0 * M_PI / box_size);

  // Loop through k-box
  for (int n_x = 0; n_x < slab_nx; n_x++) {
    for (int n_y = 0; n_y < grid_dim; n_y++) {
      for (int n_z = 0; n_z <= middle; n_z++) {
        float k_z = n_z * delta_k;
        box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] *= (fftwf_complex)(k_z * 1I);
      }
    }
  } // End looping through k box
}
