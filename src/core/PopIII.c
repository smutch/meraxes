#include "meraxes.h"
#include "misc_tools.h"
#include "stellar_feedback.h"
#include "PopIII.h"

#include <assert.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>

static float Mass_Values[MASS_BINS]; //Is there a smart way to define MASS_BINS depending on the limits of the IMF that we choose?
static float Time_Values[MASS_BINS];
static double NumberPISN;
static double MassPISN;
static double NumberSNII;
static double MassSNII;
static double MassBHs;

void initialize_time_interp_arrays()
{
  double MminIMF = run_globals.params.physics.MminIMF; 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;

  int mass_bins = 2 * (MmaxIMF - MminIMF);
  double mass_step = (MmaxIMF - MminIMF) / mass_bins;
  double mass_val = log10(MminIMF);
  int i;

  if (run_globals.mpi_rank == 0) {

    for (i = 0; i < mass_bins; i++) {
      mass_val = log10(MminIMF + (double)i * mass_step);
      Mass_Values[i] = mass_val;
      Time_Values[i] = get_StellarAge(mass_val); 
      //mlog("Mass = %f, Time = %f", MLOG_MESG, pow(10, Mass_Values[i]), pow(10, Time_Values[i]) / 1e6);
    }
  }

  // broadcast the values to all cores
  MPI_Bcast(&Mass_Values, sizeof(Mass_Values), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&Time_Values, sizeof(Time_Values), MPI_BYTE, 0, run_globals.mpi_comm);
}

void initialize_PopIII_stuff() //Initialize PopIII quantities that are easily computed just from the IMF.
{
  if (run_globals.mpi_rank == 0) {
    NumberPISN = Number_PISN();
    MassPISN = Mass_PISN();
    NumberSNII = Number_SNII();
    MassSNII = Mass_SNII();
    MassBHs = Mass_BHs();   
  }
  MPI_Bcast(&NumberPISN, sizeof(NumberPISN), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&MassPISN, sizeof(MassPISN), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&NumberSNII, sizeof(NumberSNII), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&MassSNII, sizeof(MassSNII), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&MassBHs, sizeof(MassBHs), MPI_BYTE, 0, run_globals.mpi_comm);
}

double interp_mass(double lifetime) // Lifetime in yr units!!
{
  double MminIMF = run_globals.params.physics.MminIMF; 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;

  int mass_bins = 2 * (MmaxIMF - MminIMF);
  int n_low, n_high;

  double massfinal_result;
  
  double log10lifetime = log10(lifetime);
 
  // Check if Mass is inside interpolation boundaries (That shouldn't happen, so maybe put an error message or a print
  if (log10lifetime < Time_Values[mass_bins - 1]) {
    // If it is above the upper limit, we just assume that it is the upper limit
    log10lifetime = (Time_Values[mass_bins - 1]);
    massfinal_result = (Mass_Values[mass_bins - 1]);
    } 
  else if (log10lifetime > Time_Values[0]) {
    log10lifetime = (Time_Values[0]);
    massfinal_result = (Mass_Values[0]);
    }
  else {
    for (int i = 0; i < mass_bins; i++) { //find index. You could add a safety condition here
      if ((log10lifetime <= Time_Values[i]) && (log10lifetime >= Time_Values[i+1])){
        n_low = i;
        break;
        }
      }
    
    n_high = n_low + 1;

    // Linear interpolation 
  
    massfinal_result = Mass_Values[n_low] + ((log10lifetime - Time_Values[n_low]) * (Mass_Values[n_high] - Mass_Values[n_low])) / (Time_Values[n_high] - Time_Values[n_low]);
    }

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
  double MminIMF = run_globals.params.physics.MminIMF; 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double AlphaIMF = run_globals.params.physics.AlphaIMF;
  
  double Anorm = IMFnorm(MminIMF, MmaxIMF);
  
  return Anorm * pow(StarMass, AlphaIMF);
}

double getIMF_2(double StarMass) // Don't need a second function, you could put this inside getIMF with a flag
{
  double MminIMF = run_globals.params.physics.MminIMF;
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double AlphaIMF = run_globals.params.physics.AlphaIMF;
  
  double Anorm = IMFnorm(MminIMF, MmaxIMF);
  
  return Anorm * StarMass * pow(StarMass, AlphaIMF);
}

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

double Mass_SNII(void)
{
  double result, err;
  double rel_tol = 0.01; //<- relative tolerance
  gsl_function F;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
  F.function = &getIMF_2;
  
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

double Mass_BHs(void) // Add BHs for Pop III with M>40Msol. Atm they don't do anything, it's just because we don't want Stellar population surviving!
{ 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  
  // First check if your IMF allows PISN!
  if (MmaxIMF <= MmaxSnII) //NO BHs
    return 0.0;
  
  else {
    double result_1, result_2, err_1, err_2; // You have 2 limits of integration
    double rel_tol = 0.01; //<- relative tolerance
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    F.function = &getIMF_2;
  
    if (MmaxIMF > MmaxPISN) {
      gsl_integration_qag(&F,
                      MmaxSnII,
                      MminPISN,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result_1,
                      &err_1);
      
      gsl_integration_qag(&F,
                      MmaxPISN,
                      MmaxIMF,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result_2,
                      &err_2);
      }
    else if (MmaxIMF <= MminPISN) {
      gsl_integration_qag(&F,
                      MmaxSnII,
                      MmaxIMF,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result_1,
                      &err_1);
      result_2 = 0.0;
      }
    gsl_integration_workspace_free(w);
    return (result_1 + result_2);
    }  
}

double Number_PISN(void)
{ 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  
  // First check if your IMF allows PISN!
  if (MmaxIMF < MminPISN)
    return 0.0;
  
  else {
    double result, err;
    double rel_tol = 0.01; //<- relative tolerance
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
  F.function = &getIMF;
    if (MmaxIMF >= MmaxPISN) {
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
      }
    else {
      gsl_integration_qag(&F,
                      MminPISN,
                      MmaxIMF,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &err);
      }
    gsl_integration_workspace_free(w);
    return result;
    }  
}

double Mass_PISN(void)
{
 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  
  // First check if your IMF allows PISN!
  if (MmaxIMF < MminPISN)
    return 0.0;
  
  else {
    double result, err;
    double rel_tol = 0.01; //<- relative tolerance
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF_2;
    if (MmaxIMF >= MmaxPISN) {
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
      }
    else {
      gsl_integration_qag(&F,
                      MminPISN,
                      MmaxIMF,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &err);
      }
    gsl_integration_workspace_free(w);
    return result;
    }    
}

double CCSN_PopIII_Fraction(int i_burst, int curr_snap) //Eq. 17 from Mutch et al. 2016 Result in adimensional number
{

  //mlog("curr_snap = %d, snap_SF = %d", MLOG_MESG, curr_snap, curr_snap - i_burst);
  
  /*if (last_snap <= 1) {
    mlog_error("Choose larger output snapshots");
    ABORT(EXIT_FAILURE);
  }*/
  
  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6; //You need result in yrs
  double m_min;
  double m_max;
  
  double TotalCCSN;
  
  double DeltaTime = (LTTime[curr_snap - i_burst] - LTTime[curr_snap]) * time_unit;
  double DeltaTimeSnap = (LTTime[curr_snap - i_burst - 1] - LTTime[curr_snap - i_burst]) * time_unit; //Should be correct, might need to double check that!
    
  
  if (i_burst != 0) {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = interp_mass(DeltaTime - DeltaTimeSnap / 2);
    }
    
  else {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = MmaxSnII;
    }
  
  //if (i_burst < 10)
    //mlog("m_min = %f, m_max = %f", MLOG_MESG, m_min, m_max);
    
  if (m_min > MmaxSnII) { //There are no CCSN in this snapshot!
    //mlog("m_min = %f, there are no CCSN in this snapshot", MLOG_MESG, m_min);
    return 0.0;
  }
  
  else if (m_max < MminSnII) { //There are no CCSN in this snapshot!
    //mlog("m_max = %f, there are no CCSN in this snapshot", MLOG_MESG, m_max);
    return 0.0; //Maybe put -1 or a break because this condition should stop your while loop!
  }
    
  else { 
    if (m_max > MmaxSnII) // Firstly, you are only interested in stars in the CCSN mass range
      m_max = MmaxSnII;
      
    if (m_min < MminSnII) 
      m_min = MminSnII;
      
    //mlog("m_min = %f, m_max = %f", MLOG_MESG, m_min, m_max);
  
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
  
    //TotalCCSN = Number_SNII();
    //TotalCCSN = Number_SNII() + Number_PISN(); //I am still not 100% sure if I have to consider only SNII (I believe so)
    TotalCCSN = NumberSNII + Number_PISN;
    
    //mlog("TotCCSN = %f, Frac = %f", MLOG_MESG, TotalCCSN, result / TotalCCSN);

    return result / TotalCCSN;
    //return result;
  }  
} 

double CCSN_PopIII_MassFraction(int i_burst, int curr_snap) //Eq. 22 from Mutch et al. 2016 Result in adimensional number
{

  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6; //You need result in yrs
  double m_min;
  double m_max;
  
  double TotalMassCCSN;
  
  double DeltaTime = (LTTime[curr_snap - i_burst] - LTTime[curr_snap]) * time_unit;
  double DeltaTimeSnap = (LTTime[curr_snap - i_burst - 1] - LTTime[curr_snap - i_burst]) * time_unit; //Should be correct, might need to double check that!
    
  
  if (i_burst != 0) {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = interp_mass(DeltaTime - DeltaTimeSnap / 2);
    }
    
  else {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = MmaxSnII;
    }
  
  if (m_min > MmaxSnII) { //There are no CCSN in this snapshot!
    //mlog("m_min = %f, there are no CCSN in this snapshot", MLOG_MESG, m_min);
    return 0.0;
  }
  
  else if (m_max < MminSnII) { //There are no CCSN in this snapshot!
    //mlog("m_max = %f, there are no CCSN in this snapshot", MLOG_MESG, m_max);
    return 0.0; //Maybe put -1 or a break because this condition should stop your while loop!
  }
    
  else { 
    if (m_max > MmaxSnII) // Firstly, you are only interested in stars in the CCSN mass range
      m_max = MmaxSnII;
      
    if (m_min < MminSnII) 
      m_min = MminSnII;
      
    double result, err;
    double rel_tol = 0.01; //<- relative tolerance
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF_2;
  
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
  
    //TotalCCSN = Number_SNII();
    //TotalMassCCSN = Mass_SNII() + Mass_PISN(); //I am still not 100% sure if I have to consider only SNII (I believe so)
    TotalMassCCSN = MassSNII + MassPISN;
  
    return result / TotalMassCCSN;
  }  
} 

double CCSN_PopIII_Yield(int i_burst, int curr_snap, int yield_type) //0 = Tot, 1 = Metals, 2 = Remnant. Based on Heger & Woosley 2010, No Mixing, S4 Model.
{

  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6; //You need result in yrs
  double m_min;
  double m_max;
  
  double TotalMassCCSN;
  
  double DeltaTime = (LTTime[curr_snap - i_burst] - LTTime[curr_snap]) * time_unit;
  double DeltaTimeSnap = (LTTime[curr_snap - i_burst - 1] - LTTime[curr_snap - i_burst]) * time_unit; //Should be correct, might need to double check that!
    
  
  if (i_burst != 0) {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = interp_mass(DeltaTime - DeltaTimeSnap / 2);
    }
    
  else {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = MmaxSnII;
    }
  
  if (m_min > MmaxSnII) { //There are no CCSN in this snapshot!
    //mlog("m_min = %f, there are no CCSN in this snapshot", MLOG_MESG, m_min);
    return 0.0;
  }
  
  else if (m_max < MminSnII) { //There are no CCSN in this snapshot!
    //mlog("m_max = %f, there are no CCSN in this snapshot", MLOG_MESG, m_max);
    return 0.0; //Maybe put -1 or a break because this condition should stop your while loop!
  }
    
  else {
    double Y; 
    if (m_max > MmaxSnII) // Firstly, you are only interested in stars in the CCSN mass range
      m_max = MmaxSnII;
      
    if (m_min < MminSnII) 
      m_min = MminSnII;
      
    double result, err;
    double rel_tol = 0.01; //<- relative tolerance
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF_2;
  
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
  
    //TotalCCSN = Number_SNII();
    //TotalMassCCSN = Mass_SNII() + Mass_PISN(); //I am still not 100% sure if I have to consider only SNII (I believe so)
    TotalMassCCSN = MassSNII + MassPISN;
    
    if (yield_type == 0) { // All (recycling mass) 
      if (m_max <= 30.0)
        Y = 0.88;
      else 
        Y = 0.6;
    }
    
    if (yield_type == 1) { // Metals 
      if (m_max <= 15.0)
        Y = 0.05;
      else if (m_max <= 25.0)
        Y = 0.09;
      else if (m_max <= 30.0)
        Y = 0.15;
      else
        Y = 0.0;
    }
    
    if (yield_type == 2) { //Remnant 
      if (m_max <= 30.0)
        Y = 0.12;
      else 
        Y = 0.4;
    }
  
    return result / TotalMassCCSN * Y;
  }  
} 
