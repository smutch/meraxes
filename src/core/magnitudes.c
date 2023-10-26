#include "mlog.h"
#ifdef CALC_MAGS

#include <assert.h>
#include "debug.h"
#include "magnitudes.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "parse_paramfile.h"

void init_luminosities(galaxy_t* gal)
{
  // Initialise all elements of flux arrays to TOL.
  double* inBCFlux = gal->inBCFlux;
  double* outBCFlux = gal->outBCFlux;
#if USE_MINI_HALOS
  double* inBCFluxIII = gal->inBCFluxIII;
  double* outBCFluxIII = gal->outBCFluxIII;
#endif

  for (int iSF = 0; iSF < MAGS_N; ++iSF) {
    inBCFlux[iSF] = TOL;
    outBCFlux[iSF] = TOL;
#if USE_MINI_HALOS
    inBCFluxIII[iSF] = TOL;
    outBCFluxIII[iSF] = TOL;
#endif
  }
}

void add_luminosities(mag_params_t* miniSpectra, galaxy_t* gal, int snapshot, double metals, double sfr, double new_stars)
{
  // Add luminosities when there is a burst. SFRs in principal should be in a
  // unit of M_solar/yr. However, one can convert the unit on final results
  // rather than here in order to achieve better performance.

  // Compute integer metallicity
  int Z = (int)(metals * 1000 - .5);
  if (Z < miniSpectra->minZ + 1)
    Z = miniSpectra->minZ + 1;
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
#if USE_MINI_HALOS
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6;
  if (gal->Galaxy_Population == 3){
      nZF = MAGS_N_BANDS;
      pWorking = miniSpectra->workingIII;
      pInBC = miniSpectra->inBCIII;
      pOutBC = miniSpectra->outBCIII;
      Z = 0;
      if ((bool)run_globals.params.physics.InstantSfIII)
         sfr = new_stars * time_unit; // a bit hacky... (we want new_stars / sfr is in units of year)
  }
  double* pInBCFluxIII = gal->inBCFluxIII;
  double* pOutBCFluxIII = gal->outBCFluxIII;
#endif
  double* pInBCFlux = gal->inBCFlux;
  double* pOutBCFlux = gal->outBCFlux;

  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nAgeStep = miniSpectra->targetSnap[iS];
    iA = nAgeStep - snapshot;
    if (iA >= 0) {
      iAgeBC = miniSpectra->iAgeBC[iS];
      if (iA > iAgeBC) {
        offset = (Z * nAgeStep + iA) * MAGS_N_BANDS;
        for (iF = 0; iF < MAGS_N_BANDS; ++iF){
          pOutBCFlux[iF] += sfr * pWorking[offset + iF];
#if USE_MINI_HALOS
            if (gal->Galaxy_Population == 3)
            pOutBCFluxIII[iF] += sfr * pWorking[offset + iF];
#endif 
        }
      } 
      else if (iA == iAgeBC) {
        offset = Z * MAGS_N_BANDS;
        for (iF = 0; iF < MAGS_N_BANDS; ++iF) {
          pInBCFlux[iF] += sfr * pInBC[offset + iF];
          pOutBCFlux[iF] += sfr * pOutBC[offset + iF];
#if USE_MINI_HALOS
            if (gal->Galaxy_Population == 3){
            pInBCFluxIII[iF] += sfr * pInBC[offset + iF];
            pOutBCFluxIII[iF] += sfr * pWorking[offset + iF];
          }
#endif 
        }
      } 
      else {
        offset = (Z * nAgeStep + iA) * MAGS_N_BANDS;
        for (iF = 0; iF < MAGS_N_BANDS; ++iF){
          pInBCFlux[iF] += sfr * pWorking[offset + iF];
#if USE_MINI_HALOS
            if (gal->Galaxy_Population == 3)
            pInBCFluxIII[iF] += sfr * pWorking[offset + iF];
#endif 
        }
      }
    }
    pWorking += nAgeStep * nZF;
    pInBC += nZF;
    pOutBC += nZF;
    pInBCFlux += MAGS_N_BANDS;
    pOutBCFlux += MAGS_N_BANDS;
#if USE_MINI_HALOS
    pInBCFluxIII += MAGS_N_BANDS;
    pOutBCFluxIII += MAGS_N_BANDS;
#endif
  }
}

void merge_luminosities(galaxy_t* target, galaxy_t* gal)
{
  // Sum fluexs together when a merge happens.

  double* inBCFluxTgt = target->inBCFlux;
  double* outBCFluxTgt = target->outBCFlux;
  double* inBCFlux = gal->inBCFlux;
  double* outBCFlux = gal->outBCFlux;

#if USE_MINI_HALOS
  double* inBCFluxTgtIII = target->inBCFluxIII;
  double* outBCFluxTgtIII = target->outBCFluxIII;
  double* inBCFluxIII = gal->inBCFluxIII;
  double* outBCFluxIII = gal->outBCFluxIII;
#endif

  for (int iSF = 0; iSF < MAGS_N; ++iSF) {
    inBCFluxTgt[iSF] += inBCFlux[iSF];
    outBCFluxTgt[iSF] += outBCFlux[iSF];
#if USE_MINI_HALOS
    inBCFluxTgtIII[iSF] += inBCFluxIII[iSF];
    outBCFluxTgtIII[iSF] += outBCFluxIII[iSF];
#endif
  }
}

void init_templates_mini(mag_params_t* miniSpectra,
                         char* fName,
                         char* fNameIII,
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
#if USE_MINI_HALOS
  struct sed_params_t spectraIII[MAGS_N_SNAPS];
#endif
  int nAgeStep;
  double* ageStep;
  
#if USE_MINI_HALOS
  double* ageStepIII;
#endif
  
  run_params_t* params = &run_globals.params;
  double deltaT = params->DeltaT * 1e6; // Input value in Myrs

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
#if USE_MINI_HALOS
    init_templates_rawIII(spectraIII + iS, fNameIII);
    init_filters(spectraIII + iS, betaBands, nBeta, restBands, nRest, NULL, NULL, NULL, 0, 1. + redshifts[iS]);
    spectraIII[iS].nAgeStep = nAgeStep;
    ageStepIII = (double*)malloc(nAgeStep * sizeof(double));
    if ((bool)run_globals.params.physics.InstantSfIII){
      for (int iA = 0; iA < nAgeStep; ++iA) {
        ageStepIII[iA] = LTTime[nAgeStep - iA] - LTTime[nAgeStep];
        ageStepIII[iA] += deltaT; // deltaT defined in input parameters and already converted into yrs. 
      }
	}
	else{
      for (int iA = 0; iA < nAgeStep; ++iA)
        ageStepIII[iA] = LTTime[nAgeStep - iA - 1] - LTTime[nAgeStep];
	}
    assert(ageStepIII[0] > 0.); 
    spectraIII[iS].ageStep = ageStepIII;
    shrink_templates_raw(spectraIII + iS, ageStepIII[nAgeStep - 1]);
    spectraIII[iS].igm = 0;
    if ((bool)run_globals.params.physics.InstantSfIII)
        init_templates_interpolate(spectraIII + iS);
    else
        init_templates_integrated(spectra + iS);
    spectraIII[iS].ready = (double*)malloc(nAgeStep * spectraIII[iS].nWaves * sizeof(double));
    spectraIII[iS].working = (double*)malloc(nAgeStep * spectraIII[iS].nFlux * sizeof(double));
    init_templates_workingIII(spectraIII + iS);
    init_templates_specialIII(spectraIII + iS, 1);
#endif
  }

  // Initialise mini templates
  int nSize = 0;
  int nMaxZ = spectra->nMaxZ;
  double* working;
  double* workingIII;
  size_t totalSize = 0;
  size_t totalSizeIII = 0;
  int offsetWorking = 0;
  int offsetWorkingIII = 0;
  int offsetInBC = 0;
  int offsetInBCIII = 0;
  int offsetOutBC = 0;
  int offsetOutBCIII = 0;
  int offsetWaves = 0;
  int offsetWavesIII = 0;

  // Compute size of working templates
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS)
    totalSize += targetSnap[iS];
#if USE_MINI_HALOS
  totalSizeIII = totalSize;
#endif

  totalSize *= nMaxZ * MAGS_N_BANDS;
  // Compute size of special templates
  totalSize += 2 * MAGS_N_SNAPS * nMaxZ * MAGS_N_BANDS;
  //  Compute size of wavelengths
  totalSize += 2 * MAGS_N_BANDS;
  totalSize *= sizeof(double);
  //
  working = (double*)malloc(totalSize);

#if USE_MINI_HALOS
  totalSizeIII *= MAGS_N_BANDS;
  totalSizeIII += 2 * MAGS_N_SNAPS * MAGS_N_BANDS;
  totalSizeIII += 2 * MAGS_N_BANDS;
  totalSizeIII *= sizeof(double);
  workingIII = (double*)malloc(totalSizeIII);
#endif

  // Copy working templates
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nSize = targetSnap[iS] * nMaxZ * MAGS_N_BANDS;
    memcpy(working + offsetWorking, spectra[iS].working, nSize * sizeof(double));
    offsetWorking += nSize;
#if USE_MINI_HALOS
    nSize = targetSnap[iS] * MAGS_N_BANDS;
    memcpy(workingIII + offsetWorkingIII, spectraIII[iS].working, nSize * sizeof(double));
    offsetWorkingIII += nSize;
#endif
  }
  // Copy special templates
  offsetInBC = offsetWorking;
  offsetInBCIII = offsetWorkingIII;
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nSize = nMaxZ * MAGS_N_BANDS;
    memcpy(working + offsetInBC, spectra[iS].inBC, nSize * sizeof(double));
    offsetInBC += nSize;
#if USE_MINI_HALOS
    nSize = MAGS_N_BANDS;
    memcpy(workingIII + offsetInBCIII, spectraIII[iS].inBC, nSize * sizeof(double));
    offsetInBCIII += nSize;
#endif
  }
  offsetOutBC = offsetInBC;
  offsetOutBCIII = offsetInBCIII;
  for (iS = 0; iS < MAGS_N_SNAPS; ++iS) {
    nSize = nMaxZ * MAGS_N_BANDS;
    memcpy(working + offsetOutBC, spectra[iS].outBC, nSize * sizeof(double));
    offsetOutBC += nSize;
#if USE_MINI_HALOS
    nSize = MAGS_N_BANDS;
    memcpy(workingIII + offsetOutBCIII, spectraIII[iS].outBC, nSize * sizeof(double));
    offsetOutBCIII += nSize;
#endif
  }

  // Copy wavelengths (same at each target snapshot)
  offsetWaves = offsetOutBC;
  memcpy(working + offsetWaves, spectra->centreWaves, MAGS_N_BANDS * sizeof(double));
  offsetWaves += MAGS_N_BANDS;
  memcpy(working + offsetWaves, spectra->logWaves, MAGS_N_BANDS * sizeof(double));

#if USE_MINI_HALOS
  offsetWavesIII = offsetOutBCIII;
  memcpy(workingIII + offsetWavesIII, spectraIII->centreWaves, MAGS_N_BANDS * sizeof(double));
  offsetWavesIII += MAGS_N_BANDS;
  memcpy(workingIII + offsetWavesIII, spectraIII->logWaves, MAGS_N_BANDS * sizeof(double));
#endif

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

#ifdef USE_MINI_HALOS
  miniSpectra->totalSizeIII = totalSizeIII;
  miniSpectra->workingIII = workingIII;
  miniSpectra->inBCIII = workingIII + offsetWorkingIII;
  miniSpectra->outBCIII = workingIII + offsetInBCIII;
#endif

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
#ifdef USE_MINI_HALOS
    free(spectraIII[iS].Z);
    free(spectraIII[iS].waves);
    free(spectraIII[iS].age);
    free(spectraIII[iS].raw);
    free(spectraIII[iS].nFilterWaves);
    free(spectraIII[iS].filterWaves);
    free(spectraIII[iS].filters);
    free(spectraIII[iS].integrated);
    free(spectraIII[iS].ready);
    free(spectraIII[iS].working);
    free(spectraIII[iS].inBC);
    free(spectraIII[iS].outBC);
    free(spectraIII[iS].centreWaves);
    free(spectraIII[iS].logWaves);
#endif
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
    int target_snaps[MAGS_N_SNAPS];
    int* indices;
    int count = parse_slices(params->TargetSnaps, MAGS_N_SNAPS, &indices);

    if (count != MAGS_N_SNAPS) {
      mlog_error("TargetSnaps (%d) does not match MAGS_N_SNAPS (%d)!", count, MAGS_N_SNAPS);
      ABORT(EXIT_FAILURE);
    }

    memcpy(&target_snaps, indices, sizeof(int) * count);
    free(indices);

#ifdef DEBUG
    mlog("# Target snapshots: [ ", MLOG_MESG);
    for (int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap)
      mlog("%d ", MLOG_CONT, target_snaps[i_snap]);
    mlog("]", MLOG_CONT);
#endif

    // Read beta filters
    double beta_bands[2 * MAGS_N_BANDS];
    int n_beta = 0;

    char str[STRLEN];
    char delim[] = ",";
    char* token;

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
    memcpy(str, params->PhotometricTablesDir, sizeof(str));
    strcat(str, "/sed_library.hdf5"); //Any file would be good
    char* fname = str;

    char* fnameIII = params->PhotometricTablesDir;
#if USE_MINI_HALOS
    int IMF_Type = run_globals.params.physics.PopIII_IMF;
    if (IMF_Type == 1)
      strcat(fnameIII, "/Sal500_001.hdf5"); // Expand here once you have more IMFs
    else if (IMF_Type == 2)
      strcat(fnameIII, "/Sal500_050.hdf5");
    else if (IMF_Type == 3)
      strcat(fnameIII, "/logA500_001.hdf5");
    else if (IMF_Type == 4)
      strcat(fnameIII, "/logE500_001.hdf5");
#endif
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
                        fnameIII,
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
  double* workingIII;
  ptrdiff_t offset_inBC;
  ptrdiff_t offset_inBCIII;
  ptrdiff_t offset_outBC;
  ptrdiff_t offset_outBCIII;
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

#if USE_MINI_HALOS
    workingIII = mag_params->workingIII;
    offset_inBCIII = mag_params->inBCIII - workingIII;
    offset_outBCIII = mag_params->outBCIII - workingIII;
    mag_params->workingIII = NULL;
    mag_params->inBCIII = NULL;
    mag_params->outBCIII = NULL;
#endif
  }

  MPI_Bcast(mag_params, sizeof(mag_params_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_inBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_outBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_waves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_logWaves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
#if USE_MINI_HALOS
  MPI_Bcast(&offset_inBCIII, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
  MPI_Bcast(&offset_outBCIII, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
#endif
  if (mpi_rank != MASTER){
    working = (double*)malloc(mag_params->totalSize);
#if USE_MINI_HALOS
    workingIII = (double*)malloc(mag_params->totalSizeIII);
#endif
  }
  MPI_Bcast(working, mag_params->totalSize, MPI_BYTE, MASTER, mpi_comm);

  mag_params->working = working;
  mag_params->inBC = working + offset_inBC;
  mag_params->outBC = working + offset_outBC;
  mag_params->centreWaves = working + offset_waves;
  mag_params->logWaves = working + offset_logWaves;

#if USE_MINI_HALOS
  MPI_Bcast(workingIII, mag_params->totalSizeIII, MPI_BYTE, MASTER, mpi_comm);
  mag_params->workingIII = workingIII;
  mag_params->inBCIII = workingIII + offset_inBCIII;
  mag_params->outBCIII = workingIII + offset_outBCIII;
#endif
}

void cleanup_mags(void)
{
  if (!run_globals.params.FlagMCMC)
    H5Tclose(run_globals.hdf5props.array_nmag_f_tid);
  free(run_globals.mag_params.working);
#if USE_MINI_HALOS
  free(run_globals.mag_params.workingIII);
#endif
}

#if USE_MINI_HALOS
void get_output_magnitudesIII(float* mags, galaxy_t* gal, int snapshot){
  // Convert fluxes to AB magnitudes at all target snapshots.

  // Check if ``snapshot`` is a target snapshot
  int iS;
  int* targetSnap = run_globals.mag_params.targetSnap;
  double* pInBCFlux = gal->inBCFluxIII;
  double* pOutBCFlux = gal->outBCFluxIII;

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

  } else {
    for (int i_band = 0; i_band < MAGS_N_BANDS; ++i_band) {
      mags[i_band] = 999.999f;
    }
  }
}
#endif

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
