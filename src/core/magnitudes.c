#include "meraxes.h"


void init_magnitudes(void) {
    // Initalise all relevant parameters at the first core
    if (run_globals.mpi_rank == 0) {
        printf("#***********************************************************\n");
        printf("# Compute magnitudes\n");

        // Read target snapshots
        char str[STRLEN];
        char delim[] = ",";
        char *token;
        int target_snaps[MAGS_N_SNAPS];

        memcpy(str, run_globals.params.TargetSnaps, sizeof(str));
        token = strtok(str, delim);
        for(int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap) {
            target_snaps[i_snap] = atoi(token);
            token = strtok(NULL, delim);
        }

        printf("# Target snapshots: ");
        for(int i_snap = 0; i_snap < MAGS_N_SNAPS; ++i_snap)
            printf("%d ", target_snaps[i_snap]);
        printf("\n");

        // Read rest-frame filters
        double rest_bands[2*MAGS_N_BANDS];

        memcpy(str, run_globals.params.RestBands, sizeof(str));
        token = strtok(str, delim);
        for(int i_bound = 0; i_bound < 2*MAGS_N_BANDS; ++i_bound) {
            rest_bands[i_bound] = atof(token);
            token = strtok(NULL, delim);
        }

        printf("# Rest-frame filters:\n");
        for(int i_band = 0; i_band < MAGS_N_BANDS; ++i_band)
            printf("#\t%.1f AA to %.1f\n", rest_bands[2*i_band], rest_bands[2*i_band + 1]);
        printf("#***********************************************************\n\n");
    }
}

void cleanup_mags(void) {
    H5Tclose(run_globals.hdf5props.array_nmag_f_tid);
}
