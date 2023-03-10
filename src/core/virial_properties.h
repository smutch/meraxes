#ifndef VIRIAL_PROPERTIES_H
#define VIRIAL_PROPERTIES_H

#include "meraxes.h"

// Below gives grid sizes for the interpolation arrays
#define x_int_NCFVALS 9600 //Number of CF vals in your table
#define MAX_RAD 0.51 //Largest Metal Radius value in your table (cMpc / h)
#define MAX_Rvir 0.32 // Largest Rvir value in your table (cMpc / h)

//Parameters for fitting function, later these will be put in the input parameter file
#define Alpha_ind -1.4 
#define Beta_ind 0.8
#define Gamma_ind 2.8
#define Psi_Norm 200.0

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
  double Growth_Factor(double redshift);
  double integrand_GF(double redshift);
  double GF_norm();
  double nuc(double redshift, double Halo_Mass);
  double TwoPointCF_2(double redshift, double Halo_Radius, double Radius);
  void initialize_interpCF_arrays(void);
  void initialize_interpSigma_arrays(void);
  double read_SpatialCF(double redshift, double Radius);
  double read_Sigma(double redshift, double RvirVal);
  double NLBias(double Dist_Radius, double Halo_Mass, double redshift);
  

#ifdef __cplusplus
}
#endif

#endif
