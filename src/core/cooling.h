#ifndef COOLING_H
#define COOLING_H

#define N_METALLICITIES 8
#define N_TEMPS 91
#define MIN_TEMP 4.0 // log10(T/Kelvin)
#define MAX_TEMP 8.5 // log10(T/Kelvin)

#ifdef __cplusplus
extern "C"
{
#endif

  void read_cooling_functions(void);
  double interpolate_cooling_rate(double logTemp, double logZ);
  double LTE_Mcool(double Temp, double nH);

#ifdef __cplusplus
}
#endif

#endif
