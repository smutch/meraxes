#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>

#define NMETAL 40
#define MIN_Z 0
#define MAx_Z 39
#define NAGE 2000

static double age[NAGE];
static double yield_tables[Y_NELEMENT][NMETAL*NAGE];
static double yield_tables_working[N_HISTORY_SNAPS][NMETAL][Y_NELEMENT];

void read_yield_tables(void) {
    if (run_globals.mpi_rank == 0) {
        hid_t fd;
        char fname[STRLEN];

        sprintf(fname, "%s/yield_tables.hdf5", run_globals.params.YieldDir);
        fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        // Read age [Myr]
        H5LTread_dataset_double(fd, "age", age);
        // Read total yield [1/Myr]
        H5LTread_dataset_double(fd, "total_yield", yield_tables[Y_TOTAL]);
        // Read total metal yield [1/Myr]
        H5LTread_dataset_double(fd, "total_metal_yield", yield_tables[Y_TOTAL_METAL]);

        H5Fclose(fd);
    }
    // Broadcast the values to all cores
    MPI_Bcast(age, sizeof(age), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(yield_tables, sizeof(yield_tables), MPI_BYTE, 0, run_globals.mpi_comm);
}

void compute_yield_tables(int snapshot) {
    int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;
    double *LTTime = run_globals.LTTime;
    double time_unit = run_globals.units.UnitTime_in_Megayears/run_globals.params.Hubble_h;
    double low, high;
    double *pData;

    for(int i_burst = 0; i_burst < n_bursts; ++i_burst) {
        if (i_burst > 0)
            low = ((LTTime[snapshot - i_burst - 1] + LTTime[snapshot - i_burst])/2.
                  - LTTime[snapshot - 1])*time_unit;
        else
            low = age[0];
        high = ((LTTime[snapshot - i_burst - 1] + LTTime[snapshot - i_burst])/2.
               - LTTime[snapshot])*time_unit;
        if (high < age[NAGE - 1]) {
            for(int i_element = 0; i_element < Y_NELEMENT; ++i_element) {
                pData = yield_tables[i_element];
                for(int i_metal = 0; i_metal < NMETAL; ++i_metal) {
                    yield_tables_working[i_burst][i_metal][i_element] = \
                    trapz_table(pData, age, NAGE, low, high);
                    pData += NAGE;
                }
            }
        }
        else {
            for(int i_element = 0; i_element < Y_NELEMENT; ++i_element)
                for(int i_metal = 0; i_metal < NMETAL; ++i_metal)
                    yield_tables_working[i_burst][i_metal][i_element] = 0.;
        }
    }
    
    /*
    for(int i_burst = 0; i_burst < n_bursts; ++i_burst) {
        low = yield_tables_working[i_burst][0][0];
        high = yield_tables_working[i_burst][0][1];
        printf("rank = %d, nB = %d, %.6f, %.6f, %.6f\n", run_globals.mpi_rank, 
                                                         i_burst, low, high, high/low);
    }
    printf("\n");

    for(int i_burst = 0; i_burst < n_bursts; ++i_burst) {
        low = yield_tables_working[i_burst][12][0];
        high = yield_tables_working[i_burst][12][1];
        printf("rank = %d, nB = %d, %.6f, %.6f, %.6f\n", run_globals.mpi_rank, 
                                                         i_burst, low, high, high/low);
    }
    printf("\n");

    for(int i_burst = 0; i_burst < n_bursts; ++i_burst) {
        low = yield_tables_working[i_burst][39][0];
        high = yield_tables_working[i_burst][39][1];
        printf("rank = %d, nB = %d, %.6f, %.6f, %.6f\n", run_globals.mpi_rank, 
                                                         i_burst, low, high, high/low);
    }
    printf("\n");
    */
}


double get_yield(int i_burst, double metals, int element) {
    int Z = (int)(metals*1000 - .5);
    if (Z < MIN_Z)
        Z = MIN_Z;
    else if (Z > MAx_Z)
        Z = MAx_Z;
    return yield_tables_working[i_burst][Z][element];
}
