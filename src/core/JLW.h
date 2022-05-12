/*
 * This code is an amalgamation of the requisite functions for LW background using
 * formulation of Qin2020a.
 * Written to work within Meraxes by Emanuele M. Ventura.
 */
 
 //>>>>>>Here you need to put all the includes>>>>>>//
 
 #define NU_LL (double)(3.29e15) // frequency of Lyman-limit band in Hz
 #define NU_LW (double)(2.71e15) // frequency of lower-limit of LW band in Hz
 #define PLANCK_EV (double)(4.1357e-15) // ev s
 #define R_XLy_MAX (float)(500)
 
 #include "meraxes.h"
 //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
 
 double nu_integral(int nn, float redd, float redprime, double SFRD); // Finished! 
 double nu_integrand(double nu, void* params); // Finished!
 
 void evolveLW(float zp, const double integrand_POP2[], double deriv[]);
 
 void destruct_LW();
 int init_LW();
