#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>


// These numbers are based on pre-knowledge on the yield and energy tables
#define NMETAL 40
#define MIN_Z 0
#define MAX_Z 39
#define NAGE 2000


static double age[NAGE];
static double yield_tables[Y_NELEMENT][NMETAL*NAGE];
static double yield_tables_working[N_HISTORY_SNAPS][NMETAL][Y_NELEMENT];
static double energy_tables[NMETAL*NAGE];
static double energy_tables_working[N_HISTORY_SNAPS][NMETAL];

void read_stellar_feedback_tables(void) {
    if (run_globals.mpi_rank == 0) {
        hid_t fd;
        char fname[STRLEN];
        double energy_unit = run_globals.units.UnitEnergy_in_cgs;

        sprintf(fname, "%s/stellar_feedback_tables.hdf5", run_globals.params.StellarFeedbackDir);
        fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        // Read age [Myr]
        H5LTread_dataset_double(fd, "age", age);
        // Read total yield [1/Myr]
        H5LTread_dataset_double(fd, "total_yield", yield_tables[Y_TOTAL]);
        // Read total metal yield [1/Myr]
        H5LTread_dataset_double(fd, "total_metal_yield", yield_tables[Y_TOTAL_METAL]);
        // Read energy [1/(10^10 M_solar)]
        H5LTread_dataset_double(fd, "energy", energy_tables);
        H5Fclose(fd);

        // Convert unit
        for(int i = 0; i < NMETAL*NAGE; ++i)
            energy_tables[i] /= energy_unit;
    }

    // Broadcast the values to all cores
    MPI_Bcast(age, sizeof(age), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(yield_tables, sizeof(yield_tables), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(energy_tables, sizeof(energy_tables), MPI_BYTE, 0, run_globals.mpi_comm);
}


void compute_stellar_feedback_tables(int snapshot) {
    int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;
    double *LTTime = run_globals.LTTime;
    double time_unit = run_globals.units.UnitTime_in_Megayears/run_globals.params.Hubble_h;
    double low, high;
    double *pData;

    for(int i_burst = 0; i_burst < n_bursts; ++i_burst) {
        // Assume that SN happens at the middle time of the given snapshot.
        // As an approximation adopt the value at the middle time of each 
        // snapshot for yields and energy injection.
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
            pData = energy_tables;
            for(int i_metal = 0; i_metal < NMETAL; ++i_metal) {
                energy_tables_working[i_burst][i_metal] = \
                interp(high, age, pData, NAGE) - interp(low, age, pData, NAGE);
                pData += NAGE;
            }
        }
        else {
            // When the stellar age is older than the time last grid,
            // yields and energy injection are negligible.
            for(int i_element = 0; i_element < Y_NELEMENT; ++i_element)
                for(int i_metal = 0; i_metal < NMETAL; ++i_metal)
                    yield_tables_working[i_burst][i_metal][i_element] = 0.;
            for(int i_metal = 0; i_metal < NMETAL; ++i_metal)
                energy_tables_working[i_burst][i_metal] = 0.;
        }
    }
}


double get_yield(int i_burst, double metals, int element) {
    // Convert the metallicity to an integer
    int Z = (int)(metals*1000 - .5);
    if (Z < MIN_Z)
        Z = MIN_Z;
    else if (Z > MAX_Z)
        Z = MAX_Z;
    return yield_tables_working[i_burst][Z][element];
}


double get_SN_energy(int i_burst, double metals) {
    // Convert the metallicity to an integer
    int Z = (int)(metals*1000 - .5);
    if (Z < MIN_Z)
        Z = MIN_Z;
    else if (Z > MAX_Z)
        Z = MAX_Z;
    return energy_tables_working[i_burst][Z];
}


double get_total_SN_energy(void) {
    // The last grid of the energy table is the total SNII energy,
    // and is independent to metallicity
    return energy_tables[NAGE - 1];
}

