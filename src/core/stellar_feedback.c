#include <hdf5.h>
#include <hdf5_hl.h>

#include "meraxes.h"
#include "misc_tools.h"
#include "stellar_feedback.h"
#if USE_MINI_HALOS
#include "PopIII.h"
#endif

static double age[NAGE];
static double yield_tables[NELEMENT][NMETAL * NAGE];
static double yield_tables_working[N_HISTORY_SNAPS][NMETAL][NELEMENT];
static double energy_tables[NMETAL * NAGE];
static double energy_tables_working[N_HISTORY_SNAPS][NMETAL];

static void check_n_history_snaps(void)
{
  int last_snap = run_globals.ListOutputSnaps[run_globals.NOutputSnaps - 1];

  if (last_snap <= 1) {
    mlog_error("Choose larger output snapshots");
    ABORT(EXIT_FAILURE);
  }

  // Calculate the minimum number of snapshots to fully track the SN feedback
  int i_burst;
  int isconverge;
  int n_history_snaps = 0;
  double* pData;
  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h;
  double t_begin;

  for (int i_snap = 1; i_snap < last_snap; ++i_snap) {
    for (i_burst = 0; i_burst < i_snap; ++i_burst) {
      if (i_burst > 0)
        t_begin = ((LTTime[i_snap - i_burst - 1] + LTTime[i_snap - i_burst]) / 2. - LTTime[i_snap - 1]) * time_unit;
      else
        t_begin = age[0];
      if (t_begin > age[NAGE - 1])
        break;

      isconverge = 1;
      pData = energy_tables;
      for (int i_metal = 0; i_metal < NMETAL; ++i_metal) {
        // If the energy released by SNs between the time at the beginning of
        // a snapshot and when the startburst happens is equal to the total
        // energy that a SSP can release, no energy will be released during the
        // snapshot. Therefore, this snapshot does not need to be included to
        // track the SN feedback. Otherwise, more snapshots is needed.
        if (interp(t_begin, age, pData, NAGE) != pData[NAGE - 1]) {
          isconverge = 0;
          break;
        } else
          pData += NAGE;
      }
      if (isconverge)
        break;
    }
    if (i_burst > n_history_snaps)
      n_history_snaps = i_burst;
  }

  if (n_history_snaps > N_HISTORY_SNAPS) {
    mlog_error("N_HISTORY_SNAPS is expected to be %d", n_history_snaps);
    ABORT(EXIT_FAILURE);
  }
}

void read_stellar_feedback_tables(void)
{
  if (run_globals.mpi_rank == 0) {
    hid_t fd;
    char fname[STRLEN];
    double energy_unit = run_globals.units.UnitEnergy_in_cgs;

    sprintf(fname, "%s/stellar_feedback_tables.hdf5", run_globals.params.StellarFeedbackDir);
    fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    // Read age [Myr]
    H5LTread_dataset_double(fd, "age", age);
    // Read total yield [1/Myr]
    H5LTread_dataset_double(fd, "total_yield", yield_tables[RECYCLING_FRACTION]);
    // Read total metal yield [1/Myr]
    H5LTread_dataset_double(fd, "total_metal_yield", yield_tables[TOTAL_METAL]);
    // Read energy [1/(10^10 M_solar)]
    H5LTread_dataset_double(fd, "energy", energy_tables);
    H5Fclose(fd);

    // Convert unit
    for (int i = 0; i < NMETAL * NAGE; ++i)
      energy_tables[i] /= energy_unit;

    check_n_history_snaps();
  }

  // Broadcast the values to all cores
  MPI_Bcast(age, sizeof(age), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(yield_tables, sizeof(yield_tables), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(energy_tables, sizeof(energy_tables), MPI_BYTE, 0, run_globals.mpi_comm);
}

void compute_stellar_feedback_tables(int snapshot)
{
  int n_bursts = (snapshot >= N_HISTORY_SNAPS) ? N_HISTORY_SNAPS : snapshot;
  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h;
  double t_begin, t_end;
  double* pData;

  for (int i_burst = 0; i_burst < n_bursts; ++i_burst) {
    // Assume that SN happens at the middle time of the given snapshot.
    // As an approximation adopt the value at the middle time of each
    // snapshot for yields and energy injection.
    if (i_burst > 0)
      t_begin = ((LTTime[snapshot - i_burst - 1] + LTTime[snapshot - i_burst]) / 2. - LTTime[snapshot - 1]) * time_unit;
    else
      t_begin = age[0];
    t_end = ((LTTime[snapshot - i_burst - 1] + LTTime[snapshot - i_burst]) / 2. - LTTime[snapshot]) * time_unit;
    if (t_end < age[NAGE - 1]) {
      for (int i_element = 0; i_element < NELEMENT; ++i_element) {
        pData = yield_tables[i_element];
        for (int i_metal = 0; i_metal < NMETAL; ++i_metal) {
          yield_tables_working[i_burst][i_metal][i_element] = trapz_table(pData, age, NAGE, t_begin, t_end);
          pData += NAGE;
        }
      }
      pData = energy_tables;
      for (int i_metal = 0; i_metal < NMETAL; ++i_metal) {
        energy_tables_working[i_burst][i_metal] = interp(t_end, age, pData, NAGE) - interp(t_begin, age, pData, NAGE);
        pData += NAGE;
      }
    } else {
      // When the stellar age is older than the time last grid,
      // yields and energy injection are negligible.
      for (int i_element = 0; i_element < NELEMENT; ++i_element)
        for (int i_metal = 0; i_metal < NMETAL; ++i_metal)
          yield_tables_working[i_burst][i_metal][i_element] = 0.;
      for (int i_metal = 0; i_metal < NMETAL; ++i_metal)
        energy_tables_working[i_burst][i_metal] = 0.;
    }
  }
}

static inline int get_integer_metallicity(double metals)
{
  int Z = (int)(metals * 1000 - .5);
  if (Z < MIN_Z)
    Z = MIN_Z;
  else if (Z > MAX_Z)
    Z = MAX_Z;
  return Z;
}

double get_recycling_fraction(int i_burst, double metals)
{
  // The recycling fraction equals to the yield of all elements including H & He
  return yield_tables_working[i_burst][get_integer_metallicity(metals)][RECYCLING_FRACTION];
}

double get_metal_yield(int i_burst, double metals)
{
  // The metal yield includes all elements execpt H & He
  return yield_tables_working[i_burst][get_integer_metallicity(metals)][TOTAL_METAL];
}

double get_SN_energy(int i_burst, double metals)
{
  // Convert the metallicity to an integer
  return energy_tables_working[i_burst][get_integer_metallicity(metals)];
}

double get_total_SN_energy(void)
{
  // The last grid of the energy table is the total SNII energy,
  // and is independent to metallicity
  return energy_tables[NAGE - 1];
}

#if USE_MINI_HALOS
// Stuff for Pop. III feedback

double get_SN_energy_PopIII(int i_burst,
                            int snapshot,
                            int SN_type) // SN_type = 0 -> CC, 1 -> PISN
                                         // Pop. III have higher masses so we need to account also for PISN!
{
  double NumberPISN = run_globals.NumberPISN;
  double NumberSNII = run_globals.NumberSNII;
  double Enova;
  // Core Collapse SN
  if (SN_type == 0) {
    Enova = ENOVA_CC;
    double CC_Fraction = CCSN_PopIII_Fraction(i_burst, snapshot, 0);
    return Enova * CC_Fraction * NumberSNII * 1e10 / run_globals.params.Hubble_h; // result in erg * (1e10 Msol / h)
  }
  // PISN (feedback here is contemporaneous).
  else if (SN_type == 1) {
    if (i_burst != 0) {
      mlog_error("PISN feedback is instantaneous!");
      return 0;
    }
    Enova = ENOVA_PISN;
    return Enova * NumberPISN / (NumberPISN + NumberSNII) * NumberPISN * 1e10 /
           run_globals.params.Hubble_h; // same as above
  } else {
    mlog_error("SN_type must be 0 or 1!");
    return 0;
  }
}
#endif
