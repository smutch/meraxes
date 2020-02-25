#ifndef MODIFIERS_H
#define MODIFIERS_H

// Qin, Y. et al., 2017. Dark-ages Reionization and Galaxy Formation Simulation VIII.
// Suppressed growth of dark matter halos during the Epoch of Reionization.
// Monthly Notices of the Royal Astronomical Society, 9(November), p.stx083.
// Available at: http://mnras.oxfordjournals.org/lookup/doi/10.1093/mnras/stx083.
// DO NOT CHANGE THESE UNLESS YOU KNOW WHAT THESE MEAN
#define NFIELDS 8
#define N_START 0
#define N_LOGMS 31
#define DELTA_M 0.1
#define MIN_LOGM 7.5
#define MAX_LOGM 11.5
#define M_OFFSET 0.5

typedef struct Modifier {
    float logMmin;
    float logMmax;
    float mass_mean;
    float mass_errl;
    float mass_erru;
    float ratio;
    float ratio_errl;
    float ratio_erru;
} Modifier;

#ifdef __cplusplus
extern "C" {
#endif

void read_mass_ratio_modifiers(int snapshot);
void read_baryon_frac_modifiers(int snapshot);
double interpolate_modifier(Modifier* modifier_data, double logM);

#ifdef __cplusplus
}
#endif

#endif
