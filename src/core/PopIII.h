#ifndef POPIII_H
#define POPIII_H

#include "meraxes.h"

#define MASS_BINS 1000 //Must be larger or equal than Mmax - Mmin of the IMF otherwise the thing is gonna fail!

#define MminSnII 8
#define MmaxSnII 40 

#define MminPISN 140
#define MmaxPISN 260

#ifdef __cplusplus
extern "C"
{
#endif
  
  void initialize_time_interp_arrays();
  double CCSN_PopIII_Fraction(int snapshot); 
  float interp_mass(float lifetime);
  double IMFnorm(double MminIMF, double MmaxIMF);
  double get_StellarAge(double StarMass);
  double getIMF(double StarMass);
  double Number_SNII(void);
  double Number_PISN(void);

#ifdef __cplusplus
}
#endif

#endif
