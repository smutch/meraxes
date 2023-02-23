#ifndef VIRIAL_PROPERTIES_H
#define VIRIAL_PROPERTIES_H

#include "meraxes.h"

// Below gives grid sizes for the interpolation arrays
#define x_int_NCFVALS 36000 //Number of CF vals in your table
#define MAX_RAD 500 //Largest Radius value in your table

#ifdef __cplusplus
extern "C"
{
#endif

  double Tvir_to_Mvir(double T, double z);
  double calculate_Mvir(double Mvir, int len);
  double calculate_Mvir_2(double Rvir, double redshift);
  double hubble_at_snapshot(int snapshot);
  double hubble_time(int snapshot);
  double calculate_Rvir(double Mvir, int snapshot);
  double calculate_Rvir_2(double Mvir, double redshift);
  double calculate_Vvir(double Mvir, double Rvir);
  double calculate_spin_param(halo_t* halo);
  double calculate_gasMass(int snapshot, double length);
  double calculate_zeq(double OmegaM);
  double Transfer_function(double k);
  double Growth_Factor(double redshift);
  double integrand_GF(double redshift);
  double GF_norm();
  double PowerSpectrum(double redshift, double scale);
  //double integrand_S2(double redshift, double HaloMass, double k);
  double Sigma(double redshift, double Halo_Mass);
  double SigmaNorm(double redshift);
  double nuc(double redshift, double Halo_Mass);
  double nuc_2(double redshift, double Halo_Mass);
  double R0(double redshift, double Halo_Mass);
  double TwoPointCF(double Radius, double Corr_length);
  //double integrand_2pointCF(double k);
  double TwoPointCF_2(double redshift, double Halo_Mass, double SpatialCFval);
  void initialize_interpCF_arrays(void);
  double read_SpatialCF(double redshift, double Radius);
  

#ifdef __cplusplus
}
#endif

#endif
