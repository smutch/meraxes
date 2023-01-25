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
  double calculate_gasMass(int snapshot, double length);
  double calculate_zeq(double OmegaM);
  double Transfer_function(double k);
  double Growth_Factor(double redshift);
  double integrand_GF(double redshift);
  double PowerSpectrum(double redshift, double scale);
  double integrand_S2(double redshift, double HaloMass, double k);
  double Sigma(double redshift, double HaloMass);
  double nuc(double redshift, double HaloMass);
  double R0(double redshift, double HaloMass);
  

#ifdef __cplusplus
}
#endif

#endif
