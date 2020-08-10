#include "mlog.h"
#ifdef CALC_MAGS

#include "debug.h"
#include "magnitudes.h"
#include "meraxes.h"
#include "misc_tools.h"

void init_luminosities(galaxy_t* gal)
{
  // Initialise all elements of flux arrays to TOL.
  double* inBCFlux = gal->inBCFlux;
  double* outBCFlux = gal->outBCFlux;

  for (int iSF = 0; iSF < MAGS_N; ++iSF) {
    inBCFlux[iSF] = TOL;
    outBCFlux[iSF] = TOL;
  }
}

void add_luminosities(mag_params_t* miniSpectra, galaxy_t* gal, int snapshot, double metals, double sfr)
{
  // Add luminosities when there is a burst. SFRs in principal should be in a
  // unit of M_solar/yr. However, one can convert the unit on final results
  // rather than here in order to achieve better performance.

  // Compute integer metallicity
  int Z = (int)(metals * 1000 - .5);
  if (Z < miniSpectra->minZ)
    Z = miniSpectra->minZ;
  else if (Z > miniSpectra->maxZ)
    Z = miniSpectra->maxZ;

  // Add luminosities
  int iA, iF, iS, iAgeBC;
  int offset;
  int nAgeStep;
  int nZF = miniSpectra->nMaxZ * MAGS_N_BANDS;
  double* pWorking = miniSpectra->working;
  double* pInBC = miniSpectra->inBC;
  double* pOutBC = miniSpectra->outBC;
  double* pInBCFlux = gal->inBCFlux;
  double* pOutBCFlux = gal->outBCFlux;

  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nAgeStep = miniSpectra->targetSnap[iS];
    iA = nAgeStep - snapshot;
    if (iA >= 0) {
      iAgeBC = miniSpectra->iAgeBC[iS];
      if (iA > iAgeBC) {
        offset = (Z * nAgeStep + iA) * MAGS_N_BANDS;
        for (iF = 0; iF < MAGS_N_BANDS; ++iF)
          pOutBCFlux[iF] += sfr * pWorking[offset + iF];
      } else if (iA == iAgeBC) {
        offset = Z * MAGS_N_BANDS;
        for (iF = 0; iF < MAGS_N_BANDS; ++iF) {
          pInBCFlux[iF] += sfr * pInBC[offset + iF];
          pOutBCFlux[iF] += sfr * pOutBC[offset + iF];
        }
      } else {
        offset = (Z * nAgeStep + iA) * MAGS_N_BANDS;
        for (iF = 0; iF < MAGS_N_BANDS; ++iF)
          pInBCFlux[iF] += sfr * pWorking[offset + iF];
      }
    }
    pWorking += nAgeStep * nZF;
    pInBC += nZF;
    pOutBC += nZF;
    pInBCFlux += MAGS_N_BANDS;
    pOutBCFlux += MAGS_N_BANDS;
  }
}

void merge_luminosities(galaxy_t* target, galaxy_t* gal)
{
  // Sum fluexs together when a merge happens.

  double* inBCFluxTgt = target->inBCFlux;
  double* outBCFluxTgt = target->outBCFlux;
  double* inBCFlux = gal->inBCFlux;
  double* outBCFlux = gal->outBCFlux;

  for (int iSF = 0; iSF < MAGS_N; ++iSF) {
    inBCFluxTgt[iSF] += inBCFlux[iSF];
    outBCFluxTgt[iSF] += outBCFlux[iSF];
  }
}

void init_templates_mini(mag_params_t* miniSpectra,
                         char* fName,
                         double* LTTime,
                         int* targetSnap,
                         double* redshifts,
                         double* betaBands,
                         int nBeta,
                         double* restBands,
                         int nRest,
                         double tBC)
{
  // This function first initialises all the full SED templates defined by
  // ``sed_params_t`` at given snapshots, and then transfer them to
  // ``mag_params_t``, which only contains necessary data for on-the-fly
  // luminosity calculations. It stores all arrays in a contiguous memory
  // block.

  // Initialise full templates
  int iS;
  struct sed_params_t spectra[MAGS_N_SNAPS];
  int nAgeStep;
  double* ageStep;

  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nAgeStep = targetSnap[iS];
    // Initialise raw templates
    init_templates_raw(spectra + iS, fName);
    // Initialise filters
    init_filters(spectra + iS, betaBands, nBeta, restBands, nRest, NULL, NULL, NULL, 0, 1. + redshifts[iS]);
    if (spectra[iS].nFlux != MAGS_N_BANDS) {
      mlog_error("MAGS_N_BANDS does not match!\n");
      exit(EXIT_FAILURE);
    }
    // Initialise time step
    spectra[iS].nAgeStep = nAgeStep;
    ageStep = (double*)malloc(nAgeStep * sizeof(double));
    //   -Should be in a unit of yr
    for (int iA = 0; iA < nAgeStep; ++iA)
      ageStep[iA] = LTTime[nAgeStep - iA - 1] - LTTime[nAgeStep];
    spectra[iS].ageStep = ageStep;
    //   -This function may be omitted
    shrink_templates_raw(spectra + iS, ageStep[nAgeStep - 1]);
    //   -Disable IGM absorption
    spectra[iS].igm = 0;
    // Integrate templates over given time steps
    init_templates_integrated(spectra + iS);
    // Initialise working templates
    spectra[iS].ready = (double*)malloc(spectra[iS].nZ * nAgeStep * spectra[iS].nWaves * sizeof(double));
    spectra[iS].working = (double*)malloc(spectra[iS].nMaxZ * nAgeStep * spectra[iS].nFlux * sizeof(double));
    init_templates_working(spectra + iS, NULL, NULL, -1);
    // Initialise special templates for birth cloud
    init_templates_special(spectra + iS, tBC, 1);
  }

  // Initialise mini templates
  int nSize = 0;
  int nMaxZ = spectra->nMaxZ;
  double* working;
  size_t totalSize = 0;
  int offsetWorking = 0;
  int offsetInBC = 0;
  int offsetOutBC = 0;
  int offsetWaves = 0;

  // Compute size of working templates
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS)
    totalSize += targetSnap[iS];
  totalSize *= nMaxZ * MAGS_N_BANDS;
  // Compute size of special templates
  totalSize += 2 * MAGS_N_SNAPS * nMaxZ * MAGS_N_BANDS;
  //  Compute size of wavelengths
  totalSize += 2 * MAGS_N_BANDS;
  totalSize *= sizeof(double);
  //
  working = (double*)malloc(totalSize);
  // Copy working templates
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nSize = targetSnap[iS] * nMaxZ * MAGS_N_BANDS;
    memcpy(working + offsetWorking, spectra[iS].working, nSize * sizeof(double));
    offsetWorking += nSize;
  }
  // Copy special templates
  offsetInBC = offsetWorking;
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nSize = nMaxZ * MAGS_N_BANDS;
    memcpy(working + offsetInBC, spectra[iS].inBC, nSize * sizeof(double));
    offsetInBC += nSize;
  }
  offsetOutBC = offsetInBC;
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nSize = nMaxZ * MAGS_N_BANDS;
    memcpy(working + offsetOutBC, spectra[iS].outBC, nSize * sizeof(double));
    offsetOutBC += nSize;
  }
  // Copy wavelengths (same at each target snapshot)
  offsetWaves = offsetOutBC;
  memcpy(working + offsetWaves, spectra->centreWaves, MAGS_N_BANDS * sizeof(double));
  offsetWaves += MAGS_N_BANDS;
  memcpy(working + offsetWaves, spectra->logWaves, MAGS_N_BANDS * sizeof(double));
  // Set attributes
  memcpy(miniSpectra->targetSnap, targetSnap, MAGS_N_SNAPS * sizeof(int));
  miniSpectra->nBeta = nBeta;
  miniSpectra->nRest = nRest;
  miniSpectra->minZ = spectra->minZ;
  miniSpectra->maxZ = spectra->maxZ;
  miniSpectra->nMaxZ = nMaxZ;
  miniSpectra->tBC = tBC;
  // Find the interval for birth cloud
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS)
    miniSpectra->iAgeBC[iS] = birth_cloud_interval(tBC, spectra[iS].ageStep, spectra[iS].nAgeStep);
  miniSpectra->totalSize = totalSize;
  miniSpectra->working = working;
  miniSpectra->inBC = working + offsetWorking;
  miniSpectra->outBC = working + offsetInBC;
  miniSpectra->centreWaves = working + offsetOutBC;
  miniSpectra->logWaves = working + offsetWaves;

  // Free full templates
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    free(spectra[iS].Z);
    free(spectra[iS].waves);
    free(spectra[iS].age);
    free(spectra[iS].raw);
    free(spectra[iS].nFilterWaves);
    free(spectra[iS].filterWaves);
    free(spectra[iS].filters);
    free(spectra[iS].integrated);
    free(spectra[iS].ready);
    free(spectra[iS].working);
    free(spectra[iS].inBC);
    free(spectra[iS].outBC);
    free(spectra[iS].centreWaves);
    free(spectra[iS].logWaves);
  }
}

void init_magnitudes(void)
{
  // Initalise the primary data strcture (``mag_params_t``) for on-the-fly
  // luminosity calcuations.

  int mpi_rank = run_globals.mpi_rank;
  mag_params_t* mag_params = &run_globals.mag_params;

  // Initalise all relevant parameters at the master core
  if (mpi_rank == MASTER) {
#ifdef DEBUG
    mlog("#***********************************************************", MLOG_MESG);
    mlog("# Compute magnitudes", MLOG_MESG);
#endif

    // Read target snapshots
    run_params_t* params = &run_globals.params;
    char str[STRLEN];
    char delim[] = ",";
    char* token;
    int target_snaps[MAGS_N_SNAPS];

    memcpy(str, params->TargetSnaps, sizeof(str));
    token = strtok(str, delim);
    for (int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap) {
      if (token != NULL) {
        target_snaps[i_snap] = atoi(token);
        token = strtok(NULL, delim);
      } else if (i_snap != MAGS_N_SNAPS - 1) {
        mlog_error("TargetSnaps does not match MAGS_N_SNAPS!");
        ABORT(EXIT_FAILURE);
      }
    }

#ifdef DEBUG
    mlog("# Target snapshots: ", MLOG_MESG);
    for (int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap)
      mlog("%d ", target_snaps[i_snap], MLOG_CONT);
#endif

    // Read beta filters
    double beta_bands[2 * MAGS_N_BANDS];
    int n_beta = 0;

    memcpy(str, params->BetaBands, sizeof(str));
    token = strtok(str, delim);
    for (int i_band = 0; i_band < 2 * MAGS_N_BANDS; ++i_band) {
      if (token != NULL) {
        beta_bands[i_band] = atof(token);
        token = strtok(NULL, delim);
        ++n_beta;
      } else
        break;
    }
    if (n_beta % 2 == 0)
      n_beta /= 2;
    else {
      mlog_error("Wrong BetaBands!");
      ABORT(EXIT_FAILURE);
    }

#ifdef DEBUG
    mlog("# Beta filters:", MLOG_MESG);
    for (int i_band = 0; i_band < n_beta; ++i_band)
      mlog("#\t%.1f AA to %.1f", MLOG_MESG, beta_bands[2 * i_band], beta_bands[2 * i_band + 1]);
#endif

    // Read rest-frame filters
    double rest_bands[2 * MAGS_N_BANDS];
    int n_rest = 0;

    memcpy(str, params->RestBands, sizeof(str));
    token = strtok(str, delim);
    for (int i_band = 0; i_band < 2 * MAGS_N_BANDS; ++i_band) {
      if (token != NULL) {
        rest_bands[i_band] = atof(token);
        token = strtok(NULL, delim);
        ++n_rest;
      } else
        break;
    }
    if (n_rest % 2 == 0)
      n_rest /= 2;
    else {
      mlog_error("Wrong RestBands!");
      ABORT(EXIT_FAILURE);
    }
    mlog("# Rest-frame filters:", MLOG_MESG);
    for (int i_band = 0; i_band < n_rest; ++i_band)
      mlog("#\t%.1f AA to %.1f", MLOG_MESG, rest_bands[2 * i_band], rest_bands[2 * i_band + 1]);
    //
    if (n_beta + n_rest != MAGS_N_BANDS) {
      mlog_error("Number of beta and rest-frame filters do not match MAGS_N_BANDS!", MLOG_MESG);
      ABORT(EXIT_FAILURE);
    }
    mlog("#***********************************************************", MLOG_MESG);

    // Initialise SED templates
    char* fname = params->PhotometricTablesDir;
    strcat(fname, "/sed_library.hdf5");
    // Convert time unit to yr
    int snaplist_len = params->SnaplistLength;
    double* LTTime = malloc(snaplist_len * sizeof(double));
    double time_unit = run_globals.units.UnitTime_in_Megayears / params->Hubble_h * 1e6;

    memcpy(LTTime, run_globals.LTTime, snaplist_len * sizeof(double));
    for (int i_time = 0; i_time < snaplist_len; ++i_time)
      LTTime[i_time] *= time_unit;
    //
    init_templates_mini(mag_params,
                        fname,
                        LTTime,
                        target_snaps,
                        run_globals.ZZ,
                        beta_bands,
                        n_beta,
                        rest_bands,
                        n_rest,
                        params->BirthCloudLifetime);
  }

  // Broadcast parameters to all cores
  MPI_Comm mpi_comm = run_globals.mpi_comm;
  double* working;
  ptrdiff_t offset_inBC;
  ptrdiff_t offset_outBC;
  ptrdiff_t offset_waves;
  ptrdiff_t offset_logWaves;

  if (mpi_rank == MASTER) {
    working = mag_params->working;
    offset_inBC = mag_params->inBC - working;
    offset_outBC = mag_params->outBC - working;
    offset_waves = mag_params->centreWaves - working;
    offset_logWaves = mag_params->logWaves - working;

    mag_params->working = NULL;
    mag_params->inBC = NULL;
    mag_params->outBC = NULL;
    mag_params->centreWaves = NULL;
    mag_params->logWaves = NULL;
  }

  MPI_Bcast(mag_params, sizeof(mag_params_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_inBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_outBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_waves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_logWaves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  if (mpi_rank != MASTER)
    working = (double*)malloc(mag_params->totalSize);
  MPI_Bcast(working, mag_params->totalSize, MPI_BYTE, MASTER, mpi_comm);

  mag_params->working = working;
  mag_params->inBC = working + offset_inBC;
  mag_params->outBC = working + offset_outBC;
  mag_params->centreWaves = working + offset_waves;
  mag_params->logWaves = working + offset_logWaves;
}

void cleanup_mags(void)
{
  if (!run_globals.params.FlagMCMC)
    H5Tclose(run_globals.hdf5props.array_nmag_f_tid);
  free(run_globals.mag_params.working);
}

void get_output_magnitudes(float* mags, float* dusty_mags, galaxy_t* gal, int snapshot)
{
  // Convert fluxes to AB magnitudes at all target snapshots.

  // Check if ``snapshot`` is a target snapshot
  int iS;
  int* targetSnap = run_globals.mag_params.targetSnap;
  double* pInBCFlux = gal->inBCFlux;
  double* pOutBCFlux = gal->outBCFlux;

  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    if (snapshot == targetSnap[iS])
      break;
    else {
      pInBCFlux += MAGS_N_BANDS;
      pOutBCFlux += MAGS_N_BANDS;
    }
  }
  // Correct the unit of SFRs and convert fluxes to magnitudes
  if (iS != MAGS_N_SNAPS) {
    double redshift = run_globals.ZZ[snapshot];
    double sfr_unit =
      -2.5 * log10(run_globals.units.UnitMass_in_g / run_globals.units.UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
    for (int i_band = 0; i_band < MAGS_N_BANDS; ++i_band) {
      mags[i_band] = (float)(-2.5 * log10(pInBCFlux[i_band] + pOutBCFlux[i_band]) + 8.9 + sfr_unit);
    }

    // Best fit dust--gas model from Qiu, Mutch, da Cunha et al. 2019, MNRAS, 489, 1357
    double factor = pow(calc_metallicity(gal->ColdGas, gal->MetalsColdGas) / 0.02, 1.2) * gal->ColdGas *
                    pow(gal->DiskScaleLength * 1e3, -2.0) * exp(-0.35 * redshift);
    dust_params_t dust_params = { .tauUV_ISM = 13.5 * factor,
                                  .nISM = -1.6,
                                  .tauUV_BC = 381.3 * factor,
                                  .nBC = -1.6,
                                  .tBC = run_globals.mag_params.tBC };

    double local_InBCFlux[MAGS_N_BANDS], local_OutBCFlux[MAGS_N_BANDS];
    memcpy(local_InBCFlux, pInBCFlux, sizeof(local_InBCFlux));
    memcpy(local_OutBCFlux, pOutBCFlux, sizeof(local_OutBCFlux));

    dust_absorption_approx(
      local_InBCFlux, local_OutBCFlux, run_globals.mag_params.centreWaves, MAGS_N_BANDS, &dust_params);

    for (int i_band = 0; i_band < MAGS_N_BANDS; ++i_band) {
      dusty_mags[i_band] = (float)(-2.5 * log10(local_InBCFlux[i_band] + local_OutBCFlux[i_band]) + 8.9 + sfr_unit);
    }

  } else {
    for (int i_band = 0; i_band < MAGS_N_BANDS; ++i_band) {
      mags[i_band] = 999.999f;
      dusty_mags[i_band] = 999.999f;
    }
  }
}
#endif
