#ifdef CALC_MAGS
#include "meraxes.h"

#define TOL 1e-30 // Minimum Flux


enum core {MASTER};

void init_luminosities(galaxy_t *gal) {
    double *inBCFlux = gal->inBCFlux;
    double *outBCFlux = gal->outBCFlux;

    for(int iSF = 0; iSF < MAGS_N; ++iSF) {
        inBCFlux[iSF] = TOL;
        outBCFlux[iSF] = TOL;
    }
}

void add_luminosities(mag_params_t *miniSpectra, galaxy_t *gal,
                      int snapshot, double metals, double sfr) {
    /* Add luminosities when there is a burst
     *   -SFRs should be in a unit of M_solar/yr. However, one can convert the unit on
     *    final results rather than here in order to achieve better performance 
     */
    // Compute integer metallicity
    int Z = (int)(metals*1000 - .5);
    if (Z < miniSpectra->minZ)
        Z = miniSpectra->minZ;
    else if (Z > miniSpectra->maxZ)
        Z = miniSpectra->maxZ;

    // Add luminosities
    int iA, iF, iS, iAgeBC;
    int offset;
    int nAgeStep;
    int nZF = miniSpectra->nMaxZ*MAGS_N_BANDS;
    double *pWorking = miniSpectra->working;
    double *pInBC = miniSpectra->inBC;
    double *pOutBC = miniSpectra->outBC;
    double *pInBCFlux = gal->inBCFlux;
    double *pOutBCFlux = gal->outBCFlux;

    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nAgeStep = miniSpectra->targetSnap[iS];
        iA = nAgeStep - snapshot;
        if(iA >= 0) {
            iAgeBC = miniSpectra->iAgeBC[iS];
            if (iA > iAgeBC) {
                offset = (Z*nAgeStep + iA)*MAGS_N_BANDS;
                for(iF = 0; iF < MAGS_N_BANDS; ++iF)
                    pOutBCFlux[iF] += sfr*pWorking[offset + iF];
            }
            else if (iA == iAgeBC) {
                offset = Z*MAGS_N_BANDS;
                for(iF = 0; iF < MAGS_N_BANDS; ++iF) {
                    pInBCFlux[iF] += sfr*pInBC[offset + iF];
                    pOutBCFlux[iF] += sfr*pOutBC[offset + iF];
                }
            }
            else {
                offset = (Z*nAgeStep + iA)*MAGS_N_BANDS;
                for(iF = 0; iF < MAGS_N_BANDS; ++iF)
                    pInBCFlux[iF] += sfr*pWorking[offset + iF];
            }
        }
        pWorking += nAgeStep*nZF;
        pInBC += nZF;
        pOutBC += nZF;
        pInBCFlux += MAGS_N_BANDS;
        pOutBCFlux += MAGS_N_BANDS;
    }
}

void merge_luminosities(galaxy_t *target, galaxy_t *gal) {
    /* Sum fluexs together when a merge happens
     */
    double *inBCFluxTgt = target->inBCFlux;
    double *outBCFluxTgt = target->outBCFlux;
    double *inBCFlux = gal->inBCFlux;
    double *outBCFlux = gal->outBCFlux;

    for(int iSF = 0; iSF < MAGS_N; ++iSF) {
        inBCFluxTgt[iSF] += inBCFlux[iSF];
        outBCFluxTgt[iSF] += outBCFlux[iSF];
    }
}

void init_templates_mini(mag_params_t *miniSpectra, char *fName,
                         double *LTTime, int *targetSnap, double *redshifts,
                         double *betaBands, int nBeta, double *restBands, int nRest,
                         double tBC) {
    // Initialise full templates
    int iA, iS;
    struct sed_params_t spectra[MAGS_N_SNAPS];
    int nAgeStep;
    double *ageStep;

    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nAgeStep = targetSnap[iS];
        //// Initialise raw templates
        init_templates_raw(spectra + iS, fName);
        //// Initialise filters
        init_filters(spectra + iS, betaBands, nBeta, restBands, nRest,
                     NULL, NULL, NULL, 0, 1. + redshifts[iS]);
        if (spectra[iS].nFlux != MAGS_N_BANDS) {
            printf("MAGS_N_BANDS does not match!\n");
            exit(EXIT_FAILURE);
        }
        //// Initialise time step
        spectra[iS].nAgeStep = nAgeStep;
        ageStep = (double*)malloc(nAgeStep*sizeof(double));
        ////   -Should be in a unit of yr
        for(int iA = 0; iA < nAgeStep; ++iA)
            ageStep[iA] = LTTime[nAgeStep - iA - 1] - LTTime[nAgeStep];
        spectra[iS].ageStep = ageStep;
        ////   -This function may be omitted
        shrink_templates_raw(spectra + iS, ageStep[nAgeStep - 1]);
        ////   -Disable IGM absorption
        spectra[iS].igm = 0;
        //// Integrate templates over given time steps
        init_templates_integrated(spectra + iS);
        //// Initialise working templates
        spectra[iS].ready = \
        (double*)malloc(spectra[iS].nZ*nAgeStep*spectra[iS].nWaves*sizeof(double));
        spectra[iS].working = \
        (double*)malloc(spectra[iS].nMaxZ*nAgeStep*spectra[iS].nFlux*sizeof(double));
        init_templates_working(spectra + iS, NULL, NULL, -1);
        // Initialise special templates for birth cloud
        init_templates_special(spectra + iS, tBC, 1);
    }

    // Initialise mini templates
    int nSize = 0;
    int nMaxZ = spectra->nMaxZ;
    double *working;
    size_t totalSize = 0;
    int offsetWorking = 0;
    int offsetInBC = 0;
    int offsetOutBC = 0;
    int offsetWaves = 0;

    //// Compute size of working templates
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS)
        totalSize += targetSnap[iS];
    totalSize *= nMaxZ*MAGS_N_BANDS;
    //// Compute size of special templates
    totalSize += 2*MAGS_N_SNAPS*nMaxZ*MAGS_N_BANDS;
    ///  Compute size of wavelengths
    totalSize += 2*MAGS_N_BANDS;
    totalSize *= sizeof(double);
    ////
    working = (double*)malloc(totalSize);
    //// Copy working templates
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nSize = targetSnap[iS]*nMaxZ*MAGS_N_BANDS;
        memcpy(working + offsetWorking, spectra[iS].working, nSize*sizeof(double));
        offsetWorking += nSize;
    }
    //// Copy special templates
    offsetInBC = offsetWorking;
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nSize = nMaxZ*MAGS_N_BANDS;
        memcpy(working + offsetInBC, spectra[iS].inBC, nSize*sizeof(double));
        offsetInBC += nSize;
    }
    offsetOutBC = offsetInBC;
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        nSize = nMaxZ*MAGS_N_BANDS;
        memcpy(working + offsetOutBC, spectra[iS].outBC, nSize*sizeof(double));
        offsetOutBC += nSize;
    }
    //// Copy wavelengths (same at each target snapshot)
    offsetWaves = offsetOutBC;
    memcpy(working + offsetWaves, spectra->centreWaves, MAGS_N_BANDS*sizeof(double));
    offsetWaves += MAGS_N_BANDS;
    memcpy(working + offsetWaves, spectra->logWaves, MAGS_N_BANDS*sizeof(double));
    ////
    memcpy(miniSpectra->targetSnap, targetSnap, MAGS_N_SNAPS*sizeof(int));
    miniSpectra->minZ = spectra->minZ;
    miniSpectra->maxZ = spectra->maxZ;
    miniSpectra->nMaxZ = nMaxZ;
    ////   -Find the interval for birth cloud
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS)
        miniSpectra->iAgeBC[iS] = \
        birth_cloud_interval(tBC, spectra[iS].ageStep, spectra[iS].nAgeStep);
    miniSpectra->totalSize = totalSize;
    miniSpectra->working = working;
    miniSpectra->inBC = working + offsetWorking;
    miniSpectra->outBC = working + offsetInBC;
    miniSpectra->centreWaves = working + offsetOutBC;
    miniSpectra->logWaves = working + offsetWaves;

    // Free full templates
    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
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

void init_magnitudes(void) {
    int mpi_rank = run_globals.mpi_rank;
    mag_params_t *mags_params = &run_globals.mags_params;

    // Initalise all relevant parameters at the master core
    if (mpi_rank == MASTER) {
        printf("#***********************************************************\n");
        printf("# Compute magnitudes\n");

        // Read target snapshots
        run_params_t *params = &run_globals.params;
        char str[STRLEN];
        char delim[] = ",";
        char *token;
        int target_snaps[MAGS_N_SNAPS];

        memcpy(str, params->TargetSnaps, sizeof(str));
        token = strtok(str, delim);
        for(int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap) {
            if (token != NULL) {
                target_snaps[i_snap] = atoi(token);
                token = strtok(NULL, delim);
            }
            else if (i_snap != MAGS_N_SNAPS - 1) {
                mlog_error("TargetSnaps does not match MAGS_N_SNAPS!");
                ABORT(EXIT_FAILURE);
            }
        }
        printf("# Target snapshots: ");
        for(int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap)
            printf("%d ", target_snaps[i_snap]);
        printf("\n");
        // Read beta filters
        double beta_bands[2*MAGS_N_BANDS];
        int n_beta = 0;

        memcpy(str, params->BetaBands, sizeof(str));
        token = strtok(str, delim);
        for(int i_band = 0; i_band < 2*MAGS_N_BANDS; ++i_band) {
            if (token != NULL) {
                beta_bands[i_band] = atof(token);
                token = strtok(NULL, delim);
                ++n_beta;
            }
            else
                break;
        }
        if (n_beta%2 == 0)
            n_beta /= 2;
        else {
            mlog_error("Wrong BetaBands!");
            ABORT(EXIT_FAILURE);
        }
        printf("# Beta filters:\n");
        for(int i_band = 0; i_band < n_beta; ++i_band)
            printf("#\t%.1f AA to %.1f\n", beta_bands[2*i_band], beta_bands[2*i_band + 1]);
        // Read rest-frame filters
        double rest_bands[2*MAGS_N_BANDS];
        int n_rest = 0;

        memcpy(str, params->RestBands, sizeof(str));
        token = strtok(str, delim);
        for(int i_band = 0; i_band < 2*MAGS_N_BANDS; ++i_band) {
            if (token != NULL) {
                rest_bands[i_band] = atof(token);
                token = strtok(NULL, delim);
                ++n_rest;
            }
            else
                break;
        }
        if (n_rest%2 == 0)
            n_rest /= 2;
        else {
            mlog_error("Wrong RestBands!");
            ABORT(EXIT_FAILURE);
        }
        printf("# Rest-frame filters:\n");
        for(int i_band = 0; i_band < n_rest; ++i_band)
            printf("#\t%.1f AA to %.1f\n", rest_bands[2*i_band], rest_bands[2*i_band + 1]);
        //
        if (n_beta + n_rest != MAGS_N_BANDS) {
            mlog_error("Number of beta and rest-frame filters do not match MAGS_N_BANDS!");
            ABORT(EXIT_FAILURE);
        }
        printf("#***********************************************************\n\n");

        // Initialise SED templates
        ////
        char *fname = params->PhotometricTablesDir;
        strcat(fname, "/sed_library.hdf5");
        ////Convert time unit to yr
        int snaplist_len = params->SnaplistLength;
        double *LTTime = malloc(snaplist_len*sizeof(double));
        double time_unit = run_globals.units.UnitTime_in_Megayears/params->Hubble_h*1e6;

        memcpy(LTTime, run_globals.LTTime, snaplist_len*sizeof(double));
        for(int i_time = 0; i_time < snaplist_len; ++i_time)
            LTTime[i_time] *= time_unit;
        ////
        init_templates_mini(mags_params, fname, LTTime, target_snaps, run_globals.ZZ,
                            beta_bands, n_beta, rest_bands, n_rest, params->BirthCloudLifetime);
    }

    // Broadcast parameters to all cores
    MPI_Comm mpi_comm = run_globals.mpi_comm;
    double *working;
    ptrdiff_t offset_inBC;
    ptrdiff_t offset_outBC;
    ptrdiff_t offset_waves;
    ptrdiff_t offset_logWaves;

    if (mpi_rank == MASTER) {
        working = mags_params->working;
        offset_inBC = mags_params->inBC - working;
        offset_outBC = mags_params->outBC - working;
        offset_waves = mags_params->centreWaves - working;
        offset_logWaves = mags_params->logWaves - working;

        mags_params->working = NULL;
        mags_params->inBC = NULL;
        mags_params->outBC = NULL;
        mags_params->centreWaves = NULL;
        mags_params->logWaves = NULL;
    }

    MPI_Bcast(mags_params, sizeof(mag_params_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_inBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_outBC, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_waves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    MPI_Bcast(&offset_logWaves, sizeof(ptrdiff_t), MPI_BYTE, MASTER, mpi_comm);
    if (mpi_rank != MASTER)
        working = (double*)malloc(mags_params->totalSize);
    MPI_Bcast(working, mags_params->totalSize, MPI_BYTE, MASTER, mpi_comm);

    mags_params->working = working;
    mags_params->inBC = working + offset_inBC;
    mags_params->outBC = working + offset_outBC;
    mags_params->centreWaves = working + offset_waves;
    mags_params->logWaves = working + offset_logWaves;
}

void cleanup_mags(void) {
    H5Tclose(run_globals.hdf5props.array_nmag_f_tid);
    free(run_globals.mags_params.working);
}

void get_output_magnitudes(float *target, galaxy_t *gal, int snapshot) {
    // Check if ``snapshot`` is a target snapshot
    int iS;
    int *targetSnap = run_globals.mags_params.targetSnap;
    double *pInBCFlux = gal->inBCFlux;
    double *pOutBCFlux = gal->outBCFlux;

    for(iS = 0; iS < MAGS_N_SNAPS; ++iS) {
        if (snapshot == targetSnap[iS])
            break;
        else {
            pInBCFlux += MAGS_N_BANDS;
            pOutBCFlux += MAGS_N_BANDS;
        }
    }

    if (iS != MAGS_N_SNAPS) {
        // Convert fluxes to magnitudes
        //   -Convert the unit of SFRs
        double mags[MAGS_N_BANDS];
        double sfr_unit = -2.5*log10(run_globals.units.UnitMass_in_g
                                     /run_globals.units.UnitTime_in_s
                                     *SEC_PER_YEAR/SOLAR_MASS);
        for(int i_band = 0; i_band < MAGS_N_BANDS; ++i_band)
            target[i_band] = (float)(-2.5*log10(pInBCFlux[i_band] + pOutBCFlux[i_band]) + 8.9 \
                                     + sfr_unit);
    }
}
#endif
