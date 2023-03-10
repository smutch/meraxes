#ifndef VIRIAL_PROPERTIES_H
#define VIRIAL_PROPERTIES_H

#include "meraxes.h"

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
  double hubble_at_snapshot(int snapshot);
  double hubble_time(int snapshot);
  double calculate_Rvir(double Mvir, int snapshot);
  double calculate_Vvir(double Mvir, double Rvir);
  double calculate_spin_param(halo_t* halo);
  double calculate_gasMass(int snapshot, double length);
  double NLBias(double Dist_Radius, double Halo_Mass, double redshift);
  

#ifdef __cplusplus
}
#endif

#endif
