/*
 * This code is an amalgamation of the requisite functions for LW background using
 * formulation of Qin2020a.
 * Written to work within Meraxes by Emanuele M. Ventura.
 */
 
 //>>>>>>Here you need to put all the includes>>>>>>//
 
 //#define NU_LL (double)(3.29e15) // frequency of Lyman-limit band in Hz
 //#define NU_LW (double)(2.71e15) // frequency of lower-limit of LW band in Hz
 //#define PLANCK_EV (double)(4.1357e-15) // ev s
 //#define R_XLy_MAX (float)(500)
 
 //#include "meraxes.h"
 //#include "utils.h"
 //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
 
// double* zpp_edgee; // Not a good practice as pointed by Brad but I am dumb
 
// double nu_integral(int nn, float redd, float redprime, double SFRD);
 //double nu_integrand(double nu, void* params);
 
// void evolveLW(float zp, const double integrand_POP2[], double deriv[]);
 
// void destruct_LW();
 //int init_LW();

// void ComputeJLW(int snapshot, timer_info* timer_total);
 
#ifndef COMPUTE_JLW_H
#define COMPUTE_JLW_H

#include "utils.h"
#include "meraxes.h"

#define R_XLy_MAX (float)(500)
#define NU_LL (double)(3.29e15) // frequency of Lyman-limit band in Hz
#define NU_LW (double)(2.71e15) // frequency of lower-limit of LW band in Hz
#define PLANCK_EV (double)(4.1357e-15) // ev s

#ifdef __cplusplus
extern "C"
{
#endif

  void ComputeJLW(int snapshot, timer_info* timer_total);

  double* zpp_edgee; // Not a good practice as pointed by Brad but I am dumb
  double* sum_lyn_LW;
 
  //double nu_integral(int nn, float redd, float redprime, double SFRD);
  //double nu_integrand(double nu, void* params);
  double spectral_emissivity_LW(double nu_norm, int flag, int Population);
 
  void evolveLW(float zp, const double integrand_POP2[], double deriv[]);
 
  void destruct_LW();
  int init_LW();

#ifdef __cplusplus
}
#endif

#endif
