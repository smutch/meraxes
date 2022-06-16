/*
 * This code is an amalgamation of the requisite functions for LW background using
 * formulation of Qin2020a.
 * Written to work within Meraxes by Emanuele M. Ventura.
 */
 
#ifndef COMPUTE_JLW_H
#define COMPUTE_JLW_H

#include "utils.h"
#include "meraxes.h"

#define R_XLy_MAX (float)(500)
#define NU_LL (double)(3.29e15) 
#define NU_LW (double)(2.71e15) 
#define PLANCK_EV (double)(4.1357e-15)

#ifdef __cplusplus
extern "C"
{
#endif

  void ComputeJLW(int snapshot, timer_info* timer_total);

  double* zpp_edgee; 
  double* sum_lyn_LW;
  extern double* sum_lyn_LW;
 
  double spectral_emissivity_LW(double nu_norm, int Population);
 
  void evolveLW(float zp, const double integrand_POP2[], const double StarF_POP2[], double deriv[]);
 
  void destruct_LW();
  int init_LW();

#ifdef __cplusplus
}
#endif

#endif
