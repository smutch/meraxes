#ifndef VIRIAL_PROPERTIES_H
#define VIRIAL_PROPERTIES_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  double Tvir_to_Mvir(double T, double z);
  double calculate_Mvir(double Mvir, int len);
  double hubble_at_snapshot(int snapshot);
  double hubble_time(int snapshot);
  double calculate_Rvir(double Mvir, int snapshot);
  double calculate_Vvir(double Mvir, double Rvir);
  double calculate_spin_param(halo_t* halo);
  double Vvir_to_Tvir(double Vvir, int halo_type);
  double Vvir_to_Mvir(double Vvir, double redshift, int halo_type);
  
#ifdef __cplusplus
}
#endif

#endif
