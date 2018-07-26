#include "meraxes.h"

void init_magnitudes(void) {
    struct sed_params spectra;
    init_templates_raw(&spectra, "./");
}

void cleanup_mags(void) {
    H5Tclose(run_globals.hdf5props.array_nmag_f_tid);
}
