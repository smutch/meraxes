#include "meraxes.h"
#include "misc_tools.h"
#include "stellar_feedback.h"
#include "PopIII.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>

static float Mass_Values[MASS_BINS];
static float Time_Values[MASS_BINS];

void initialize_time_interp_arrays()
{
  double MminIMF = run_globals.params.physics.MminIMF; //remember to put these new parameters for IMF!
  double MmaxIMF = run_globals.params.physics.MmaxIMF;

  int mass_bins = MmaxIMF - MminIMF;
  double mass_val = log10(MminIMF);
  int i;
  
  double time_val;

  if (run_globals.mpi_rank == 0) {

    for (i = 0; i < mass_bins; i++) {
      mass_val = log10(MminIMF + i); //Summing double + int! Is it a problem? (Maybe double check it later)
      Mass_Values[i] = mass_val; //&Mass_Values[i] ?
      Time_Values[i] = get_StellarAge(mass_val); //&Time_Values[i] ?
    }
  }

  // broadcast the values to all cores
  MPI_Bcast(&Mass_Values, sizeof(Mass_Values), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&Time_Values, sizeof(Time_Values), MPI_BYTE, 0, run_globals.mpi_comm);
}

float interp_mass(float lifetime) // CHECK THIS!!! 
{

  int n_low, n_high;
  double MminIMF = run_globals.params.physics.MminIMF;
  double log10MminIMF = log10(MminIMF);
  double MinLifetime = get_StellarAge(log10MminIMF); // Should be equivalent of using Time_Values[0]

  float massfinal_result;

  // Check if Mass is inside interpolation boundaries (That shouldn't happen, so maybe put an error message or a print
  if (lifetime > 0.999 * Time_Values[MASS_BINS - 1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    lifetime = (float)(Time_Values[MASS_BINS - 1] * 0.999);
  } else if (lifetime < Time_Values[0]) {
    return 0.0;
  }
  for (int i = 0; i < MASS_BINS; i++) { //find index. You could add a safety condition here
   if ((lifetime > Time_Values[i]) && (lifetime < Time_Values[i+1])){
     n_low = i;
     break;
     }
   }
  n_high = n_low + 1;


  // Linear interpolation (you can do that because input values are in log10 units! STOP HERE! You need to figure out this!
  
  massfinal_result = Mass_Values[n_low] + ((lifetime - Time_Values[n_low]) * (Mass_Values[n_high] - Mass_Values[n_low])) / (Time_Values[n_high] - Time_Values[n_low]);

  return pow(10, massfinal_result); //Return result in SolarMass
}

double integrand_IMFnorm(double StarMass) // You might want a case 2 for a different type of IMF (with Characteristic mass)
{
  double AlphaIMF = run_globals.params.physics.AlphaIMF; //remember to put these new parameters for IMF!
  
  return StarMass * pow(StarMass, AlphaIMF);
}

double IMFnorm(double Mmin_IMF, double Mmax_IMF) //get normalization of Pop III IMF
{

  double norma, normaerr;
  double rel_tol = 0.01; //<- relative tolerance
  gsl_function F;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
  F.function = &integrand_IMFnorm;
  
  gsl_integration_qag(&F,
                      Mmin_IMF,
                      Mmax_IMF,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &norma,
                      &normaerr);
  gsl_integration_workspace_free(w);

  return 1.0 / norma;
}

double getIMF(double StarMass)
{
  double MminIMF = run_globals.params.physics.MminIMF; //remember to put these new parameters for IMF!
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double AlphaIMF = run_globals.params.physics.AlphaIMF;
  
  double Anorm = IMFnorm(MminIMF, MmaxIMF);
  
  return Anorm * pow(StarMass, AlphaIMF);
}

/*double get_StellarAge(double StarMass) //Stellar lifetime in Myr 
{
  int PopIIIAgePrescription = run_globals.params.physics.PopIIIAgePrescription; //Remember to put this as a parameter!
  double logStarMass = log10(StarMass);
  switch (PopIIIAgePrescription) {
      case 1: //stellar lifetimes by Schaerer 2002 (popIII) model at met=0 with strong mass loss, mass range [80-1000] Msun
        double a0 = 8.795;
        double a1 = -1.797;
        double a2 = 0.332;
        double a3 = 0;
        break;
      case 2: //stellar lifetimes by Schaerer 2002 (popIII) model at met=0 with no mass loss, mass range [9-500] Msun
        double a0 = 9.785;
        double a1 = -3.759;
        double a2 = 1.413;
        double a3 = -0.186;
        break;
      default:
        mlog_error("Unrecognised value for PopIIIAgePrescription parameter! Defaulting to Strong Mass Loss case");
        double a0 = 8.795;
        double a1 = -1.797;
        double a2 = 0.332;
        double a3 = 0;
        break;
    }
  return pow(10, a0 + a1 * logStarMass + a2 * logStarMass * logStarMass + a3 * pow(logStarMass, 3)) / 1e6;   
}*/

double get_StellarAge(double StarMass) //Star Mass in log10(Msol) to get tstar in log10(yr). Use this so you can do a linear interpolation!!! 
{
  int PopIIIAgePrescription = run_globals.params.physics.PopIIIAgePrescription; //Remember to put this as a parameter!
  double a0, a1, a2, a3;
  switch (PopIIIAgePrescription) {
      case 1: //stellar lifetimes by Schaerer 2002 (popIII) model at met=0 with strong mass loss, mass range [80-1000] Msun
        a0 = 8.795;
        a1 = -1.797;
        a2 = 0.332;
        a3 = 0;
        break;
      case 2: //stellar lifetimes by Schaerer 2002 (popIII) model at met=0 with no mass loss, mass range [9-500] Msun
        a0 = 9.785;
        a1 = -3.759;
        a2 = 1.413;
        a3 = -0.186;
        break;
      default:
        mlog_error("Unrecognised value for PopIIIAgePrescription parameter! Defaulting to Strong Mass Loss case");
        a0 = 8.795;
        a1 = -1.797;
        a2 = 0.332;
        a3 = 0;
        break;
    }
  return a0 + a1 * StarMass + a2 * StarMass * StarMass + a3 * pow(StarMass, 3);   
}

double Number_SNII(void)
{
  double result, err;
  double rel_tol = 0.01; //<- relative tolerance
  gsl_function F;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
  F.function = &getIMF;
  
  gsl_integration_qag(&F,
                      MminSnII,
                      MmaxSnII,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &err);
  gsl_integration_workspace_free(w);

  return result;  
}

double Number_PISN(void)
{
  double result, err;
  double rel_tol = 0.01; //<- relative tolerance
  gsl_function F;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
  F.function = &getIMF;
  
  gsl_integration_qag(&F,
                      MminPISN,
                      MmaxPISN,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &err);
  gsl_integration_workspace_free(w);

  return result;  
}

double CCSN_PopIII_Fraction(int Snapshot) //Eq. 17 from Mutch et al. 2016 YOU ARE HERE! NEED TO WRITE THIS FUNCTION!!!
{

  int last_snap = run_globals.ListOutputSnaps[run_globals.NOutputSnaps - 1];
  
  mlog("last_snap = %d", MLOG_MESG, last_snap);
  
  if (last_snap <= 1) {
    mlog_error("Choose larger output snapshots");
    ABORT(EXIT_FAILURE);
  }
  
  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6; //You need result in yrs
  double m_min;
  double m_max;
  
  double TotalCCSN;
  
  double DeltaTime = (LTTime[last_snap] - LTTime[snapshot]) * time_unit;
  double DeltaTimeSnap = (LTTime[Snapshot] - LTTime[Snapshot - 1]) * time_unit; //Should be correct, might need to double check that!
  
  mlog("DeltaT Snap, DeltaTSF (Myr) = %f", MLOG_MESG, DeltaTimeSnap / 1e6, DeltaTime / 1e6);

  m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
  m_max = interp_mass(DeltaTime - DeltaTimeSnap / 2);
  
  mlog("m_min = %f, m_max = %f", MLOG_MESG, m_min, m_max);
  
  double result, err;
  double rel_tol = 0.01; //<- relative tolerance
  gsl_function F;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
  F.function = &getIMF;
  
  gsl_integration_qag(&F,
                      m_min,
                      m_max,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &err);
  gsl_integration_workspace_free(w);
  
  TotalCCSN = Number_SNII();
  
  mlog("TotCCSN, Frac = %f", MLOG_MESG, TotalCCSN, result / TotalCCSN);

  return result / TotalCCSN;  
} 
