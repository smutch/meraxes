#ifndef COOLING_H
#define COOLING_H

#define N_METALLICITIES 8
#define N_TEMPS 91
#define MIN_TEMP 4.0 // log10(T/Kelvin)
#define MAX_TEMP 8.5 // log10(T/Kelvin)

// metallicies with repect to solar.
// These will be converted to absolute metallicities by adding log10(Z_sun), Zsun=0.02
static double metallicities[N_METALLICITIES] = {
    -5.0, // This is actually primordial in the tables but that is log10(0) = -infinity
    -3.0,
    -2.0,
    -1.5,
    -1.0,
    -0.5,
    +0.0,
    +0.5
};

static char group_name[N_METALLICITIES][6] = {
    "mzero",
    "m-30",
    "m-20",
    "m-15",
    "m-10",
    "m-05",
    "m-00",
    "m+05"
};

static double cooling_rate[N_METALLICITIES][N_TEMPS];


#ifdef __cplusplus
extern "C" {
#endif

void read_cooling_functions(void);
double interpolate_cooling_rate(double logTemp, double logZ);

#ifdef __cplusplus
}
#endif

#endif
