#ifndef POPIII_H
#define POPIII_H

#include "meraxes.h"

#define MASS_BINS 1000 //Must be larger or equal than Mmax - Mmin of the IMF otherwise the thing is gonna fail!

#define MminSnII 8
#define MmaxSnII 40 

#define MminPISN 140
#define MmaxPISN 260

#ifdef _POPIII_C

static double NumberPISN;
static double MassPISN;
static double NumberSNII;
static double MassSNII;
static double MassBHs;

#else

extern double NumberPISN;
extern double MassPISN;
extern double NumberSNII;
extern double MassSNII;
extern double MassBHs;

#ifdef __cplusplus
extern "C"
{
#endif
  
  void initialize_time_interp_arrays(void);
  void initialize_PopIII_stuff(void);
  double CCSN_PopIII_Fraction(int i_burst, int curr_snap); 
  double CCSN_PopIII_MassFraction(int i_burst, int curr_snap);
  double interp_mass(double lifetime);
  double IMFnorm(double MminIMF, double MmaxIMF);
  double get_StellarAge(double StarMass);
  double getIMF(double StarMass);
  double getIMF_2(double StarMass);
  double Number_SNII(void);
  double Mass_SNII(void);
  double Number_PISN(void);
  double Mass_PISN(void);
  double Mass_BHs(void); 
  double CCSN_PopIII_Yield(int i_burst, int curr_snap, int yield_type);

#ifdef __cplusplus
}
#endif

#endif
