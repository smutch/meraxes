#ifndef POPIII_H
#define POPIII_H

#include "meraxes.h"

#define MminSnII 8
#define MmaxSnII 40 

#define MminPISN 140
#define MmaxPISN 260

#define IMF_MASS_STEP 0.5

#ifdef __cplusplus
extern "C"
{
#endif
  
  void initialize_time_interp_arrays(double MminIMF, double MmaxIMF);
  void initialize_PopIII(void);
  double CCSN_PopIII_Fraction(int i_burst, int curr_snap, int flagMW); 
  double interp_mass(double lifetime);
  double IMFnorm(double MminIMF, double MmaxIMF);
  double get_StellarAge(double StarMass);
  double getIMF(double StarMass);
  double getIMF_massweighted(double StarMass);
  double Number_SNII(void);
  double Mass_SNII(void);
  double Number_PISN(void);
  double Mass_PISN(void);
  double Mass_BHs(void); 
  double CCSN_PopIII_Yield(int i_burst, int curr_snap, int yield_type);
  double PISN_PopIII_Yield(int yield_type);
  double NLBias(double Dist_Radius, double Halo_Mass, double redshift);

#ifdef __cplusplus
}
#endif

#endif
