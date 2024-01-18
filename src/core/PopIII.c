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

void initialize_time_interp_arrays(double MminIMF, double MmaxIMF)
{
  int mass_bins = (MmaxIMF - MminIMF) / IMF_MASS_STEP;

  run_globals.Mass_Values = malloc(sizeof(float) * mass_bins);
  run_globals.Time_Values = malloc(sizeof(float) * mass_bins);

  if (run_globals.mpi_rank == 0) {
    for (int i = 0; i < mass_bins; i++) {
      run_globals.Mass_Values[i] = log10(MminIMF + (double)i * IMF_MASS_STEP);
      run_globals.Time_Values[i] = get_StellarAge(run_globals.Mass_Values[i]); 
    }
  }

  // broadcast the values to all cores
  MPI_Bcast(run_globals.Mass_Values, mass_bins, MPI_FLOAT, 0, run_globals.mpi_comm);
  MPI_Bcast(run_globals.Time_Values, mass_bins, MPI_FLOAT, 0, run_globals.mpi_comm);
}

void initialize_PopIII() //Initialize PopIII quantities that are easily computed just from the IMF.
{

  int IMF_Type = run_globals.params.physics.PopIII_IMF;
  double MminIMF;
  double MmaxIMF;
  double AlphaIMF;
  double McharIMF;
  double SigmaIMF;
  double NionBaryIII;
  
  switch (IMF_Type) { //Use the one for which you have the spectra from Raiter et al. (2010).
    case 1:
      MminIMF = 1.0;
      MmaxIMF = 500.0;
      AlphaIMF = -2.35;
      NionBaryIII = 22000.0;
      break;
    case 2:
      MminIMF = 50.0;
      MmaxIMF = 500.0;
      AlphaIMF = -2.35;
      NionBaryIII = 72000.0;
      break;
    case 3:
      MminIMF = 1.0;
      MmaxIMF = 500.0;
      McharIMF = 10.0;
      SigmaIMF = 1.0;
      NionBaryIII = 47600.0;
      break;
    case 4:
      MminIMF = 1.0;
      MmaxIMF = 500.0;
      McharIMF = 60.0;
      SigmaIMF = 1.0;
      NionBaryIII = 71000.0;
      break;
    default:
      mlog_error("Unrecognised value for PopIII_IMF! Defaulting to Salpeter 1.");
      MminIMF = 1.0;
      MmaxIMF = 100.0;
      AlphaIMF = -2.35;
      NionBaryIII = 22000.0;
      break;
  }
  
  //double MminIMF = run_globals.params.physics.MminIMF; 
  //double MmaxIMF = run_globals.params.physics.MmaxIMF;
  assert(MminIMF<MmaxIMF);
  
  run_globals.params.physics.MminIMF = MminIMF;
  run_globals.params.physics.MmaxIMF = MmaxIMF;
  
  if (IMF_Type < 3) //Salpeter IMF
    run_globals.params.physics.AlphaIMF = AlphaIMF;
  else { //logNormal IMF
    run_globals.params.physics.McharIMF = McharIMF;
    run_globals.params.physics.SigmaIMF = SigmaIMF;
    }
  run_globals.params.physics.ReionNionPhotPerBaryIII = NionBaryIII; // Assign this value because this depends on the choice of the IMF

  initialize_time_interp_arrays(MminIMF, MmaxIMF);
  
  run_globals.IMFnorm = IMFnorm(MminIMF, MmaxIMF);
  run_globals.NumberPISN = Number_PISN();
  run_globals.MassPISN = Mass_PISN();
  run_globals.NumberSNII = Number_SNII();
  run_globals.MassSNII = Mass_SNII();
  run_globals.MassBHs = Mass_BHs();
  mlog("Init IMF quantities: NPISN = %f, MPISN = %f, NSNII = %f, MSNII = %f, MBHs = %f", MLOG_MESG, run_globals.NumberPISN, run_globals.MassPISN, run_globals.NumberSNII, run_globals.MassSNII, run_globals.MassBHs);    
}

double interp_mass(double lifetime) // Lifetime in yr units!!
{
  double MminIMF = run_globals.params.physics.MminIMF; 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;

  int mass_bins = (MmaxIMF - MminIMF) / IMF_MASS_STEP;
  int n_low, n_high;

  double massfinal_result;
  
  double log10lifetime = log10(lifetime);
 
  // Check if Mass is inside interpolation boundaries
  if (log10lifetime < run_globals.Time_Values[mass_bins - 1]) {
    // If it is above the upper limit, we just assume that it is the upper limit
    log10lifetime = (run_globals.Time_Values[mass_bins - 1]);
    massfinal_result = (run_globals.Mass_Values[mass_bins - 1]);
    } 
  else if (log10lifetime > run_globals.Time_Values[0]) {
    log10lifetime = (run_globals.Time_Values[0]);
    massfinal_result = (run_globals.Mass_Values[0]);
    }
  else {
    for (int i = 0; i < mass_bins; i++) { //find index.
      if ((log10lifetime <= run_globals.Time_Values[i]) && (log10lifetime > run_globals.Time_Values[i+1])){
        n_low = i;
        break;
        }
      }   
    n_high = n_low + 1;
    // Linear interpolation 
    massfinal_result = run_globals.Mass_Values[n_low] + ((log10lifetime - run_globals.Time_Values[n_low]) * (run_globals.Mass_Values[n_high] - run_globals.Mass_Values[n_low])) / (run_globals.Time_Values[n_high] - run_globals.Time_Values[n_low]);
    }
  return pow(10, massfinal_result); //Return result in SolarMass
}

double integrand_IMFnorm(double StarMass) // You might want a case 2 for a different type of IMF (with Characteristic mass)
{
  int IMF_Type = run_globals.params.physics.PopIII_IMF;
  
  if (IMF_Type < 3) { //Salpeter IMF
    double AlphaIMF = run_globals.params.physics.AlphaIMF;
    return StarMass * pow(StarMass, AlphaIMF);
    }
  else { //logN
    double McharIMF = run_globals.params.physics.McharIMF;
    double SigmaIMF = run_globals.params.physics.SigmaIMF;
    return exp((-0.5 / (SigmaIMF * SigmaIMF) * log(StarMass / McharIMF) * log(StarMass / McharIMF)));
    }
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
  int IMF_Type = run_globals.params.physics.PopIII_IMF;
  double Anorm = run_globals.IMFnorm;
  
  if (IMF_Type < 3) { //Salpeter IMF
    double AlphaIMF = run_globals.params.physics.AlphaIMF;
    return Anorm * pow(StarMass, AlphaIMF);
    }
  else { //logNormal
    double McharIMF = run_globals.params.physics.McharIMF;
    double SigmaIMF = run_globals.params.physics.SigmaIMF;
    return 1.0 / StarMass * exp(log(Anorm) - 0.5 / (SigmaIMF * SigmaIMF) * log(StarMass / McharIMF) * log(StarMass / McharIMF));
    }
}

double getIMF_massweighted(double StarMass) 
{
  return StarMass * getIMF(StarMass);
}

double get_StellarAge(double StarMass) //Star Mass in log10(Msol) to get tstar in log10(yr).  
{
  int PopIIIAgePrescription = run_globals.params.physics.PopIIIAgePrescription; 
  if ((PopIIIAgePrescription == 1) && (run_globals.params.physics.MmaxIMF > 500.0)) {
    mlog_error("MmaxIMF must be smaller than 500 if PopIIIAgePrescription == 1");
    mlog("MmaxIMF = %f", MLOG_MESG, run_globals.params.physics.MmaxIMF);
    ABORT(EXIT_FAILURE); 
    }
  double a0, a1, a2, a3;
  switch (PopIIIAgePrescription) {
      case 1: //stellar lifetimes by Schaerer 2002 (popIII) model at met=0 with strong mass loss, mass range [80-1000] Msun
        a0 = 8.795;
        a1 = -1.797;
        a2 = 0.332;
        a3 = 0;
        break;
      case 2: //stellar lifetimes by Schaerer 2002 (popIII) model at met=0 with no mass loss, mass range [5-500] Msun
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
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double MminIMF = run_globals.params.physics.MminIMF;
  if ((MmaxIMF < MminSnII) || (MminIMF > MmaxSnII))
    return 0.0;
  
  else {
    double result, err;
    double rel_tol = 0.01; 
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF;
    
    gsl_integration_qag(&F,
                        MminIMF < MminSnII ? MminSnII : MminIMF,
                        MmaxIMF > MmaxSnII ? MmaxSnII : MmaxIMF,
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
}

double Mass_SNII(void)
{
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double MminIMF = run_globals.params.physics.MminIMF;
  if ((MmaxIMF < MminSnII) || (MminIMF > MmaxSnII))
    return 0.0;
  
  else {
    double result, err;
    double rel_tol = 0.01; 
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF_massweighted;
    
    gsl_integration_qag(&F,
                        MminIMF < MminSnII ? MminSnII : MminIMF,
                        MmaxIMF > MmaxSnII ? MmaxSnII : MmaxIMF,
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
}

double Mass_BHs(void) // Add BHs for Pop III with M>40Msol. Atm they don't do anything, it's just because we don't want Stellar population surviving!
{ 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double MminIMF = run_globals.params.physics.MminIMF;
  
  // First check if your IMF allows BHs as failed SN!
  if (MmaxIMF <= MmaxSnII) //NO BHs
    return 0.0;
  
  else {
    double result_1, result_2, err_1, err_2; // You have 2 limits of integration
    double rel_tol = 0.01;
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    F.function = &getIMF_massweighted;
  
    if (MminIMF < MminPISN) {
        gsl_integration_qag(&F,
                    MminIMF < MmaxSnII ? MmaxSnII : MminIMF,
                    MmaxIMF > MminPISN ? MminPISN : MmaxIMF,
                    0,
                    rel_tol,
                    1000,
                    GSL_INTEG_GAUSS61,
                    w,
                    &result_1,
                    &err_1);
        }
    else
	result_1 = 0.;

    if (MmaxIMF > MmaxPISN) {
         gsl_integration_qag(&F,
                    MminIMF < MmaxPISN ? MmaxPISN : MminIMF,
                    MmaxIMF,
                    0,
                    rel_tol,
                    1000,
                    GSL_INTEG_GAUSS61,
                    w,
                    &result_2,
                    &err_2);
        }
    else
        result_2 = 0.;

    gsl_integration_workspace_free(w);
    return (result_1 + result_2);
    }
}

double Number_PISN(void)
{ 
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double MminIMF = run_globals.params.physics.MminIMF;
  
  // First check if your IMF allows PISN!
  if (MmaxIMF < MminPISN)
    return 0.0;
  
  else {
    double result, err;
    double rel_tol = 0.01;
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF;
    gsl_integration_qag(&F,
                    MminIMF < MminPISN ? MminPISN : MminIMF,
                    MmaxIMF > MmaxPISN ? MmaxPISN : MmaxIMF,
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
}

double Mass_PISN(void)
{
  double MmaxIMF = run_globals.params.physics.MmaxIMF;
  double MminIMF = run_globals.params.physics.MminIMF;
  
  // First check if your IMF allows PISN!
  if (MmaxIMF < MminPISN)
    return 0.0;
  
  else {
    double result, err;
    double rel_tol = 0.01; //<- relative tolerance
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  
    F.function = &getIMF_massweighted;
    gsl_integration_qag(&F,
                    MminIMF < MminPISN ? MminPISN : MminIMF,
                    MmaxIMF > MmaxPISN ? MmaxPISN : MmaxIMF,
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
}

double CCSN_PopIII_Number(int i_burst, int curr_snap, int flagMW)
{
  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6; //You need result in yrs
  double m_min;
  double m_max;
  
  double DeltaTime = (LTTime[curr_snap - i_burst] - LTTime[curr_snap]) * time_unit;
  double DeltaTimeSnap = (LTTime[curr_snap - i_burst - 1] - LTTime[curr_snap - i_burst]) * time_unit;
  double DeltaTimeSnap2 = (LTTime[curr_snap - i_burst] - LTTime[curr_snap - i_burst + 1]) * time_unit;
  
  if (i_burst != 0) {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = interp_mass(DeltaTime - DeltaTimeSnap2 / 2);
    }
    
  else {
    m_min = interp_mass(DeltaTime + DeltaTimeSnap / 2);
    m_max = MmaxSnII;
    }
  
  if (m_min > MmaxSnII) // There are no CCSN in this snapshot!
    return 0.0;
  
  else if (m_max < MminSnII) // There are no CCSN in this snapshot!
    return 0.0; // Maybe put -1 or a break because this condition should stop your while loop!
    
  else {
    if (m_max > MmaxSnII) // Firstly, you are only interested in stars in the CCSN mass range
      m_max = MmaxSnII;
      
    if (m_min < MminSnII) 
      m_min = MminSnII;
    
    double result, err;
    double rel_tol = 0.01;
    gsl_function F;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    
    if (flagMW == 0)
      F.function = &getIMF;
    
    else if (flagMW == 1)
      F.function = &getIMF_massweighted;
  
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

    return result;
  }  
}

double CCSN_PopIII_Fraction(int i_burst, int curr_snap, int flagMW) //from Mutch et al. 2016 flag == 0 Eq. 17, flag == 1 Eq. 22  Result in adimensional number
{

  double result = CCSN_PopIII_Number(i_burst, curr_snap, flagMW);
  if (result < REL_TOL)
      return 0.0;

  double TotalSN;

  if (flagMW == 0) { 
    double NumberPISN = run_globals.NumberPISN;
    double NumberSNII = run_globals.NumberSNII;
    TotalSN = NumberSNII + NumberPISN;
  }
  else if (flagMW == 1) {
    double MassPISN = run_globals.MassPISN;
    double MassSNII = run_globals.MassSNII;
    TotalSN = MassPISN + MassSNII;
  }

  return result / TotalSN;
} 

double CCSN_PopIII_Yield(int i_burst, int curr_snap, int yield_type) //0 = Tot, 1 = Metals, 2 = Remnant. Based on Heger & Woosley 2010, No Mixing, S4 Model.
{

  double result = CCSN_PopIII_Number(i_burst, curr_snap, 1);
  if (result < REL_TOL)
      return 0.0;

  double* LTTime = run_globals.LTTime;
  double time_unit = run_globals.units.UnitTime_in_Megayears / run_globals.params.Hubble_h * 1e6; //You need result in yrs
  double m_max;
  
  double TotalMassSN = run_globals.MassSNII;
  
  double DeltaTime = (LTTime[curr_snap - i_burst] - LTTime[curr_snap]) * time_unit;
  double DeltaTimeSnap2 = (LTTime[curr_snap - i_burst] - LTTime[curr_snap - i_burst + 1]) * time_unit;
  
  if (i_burst != 0)
    m_max = interp_mass(DeltaTime - DeltaTimeSnap2 / 2);
  else
    m_max = MmaxSnII;
  if (m_max > MmaxSnII) // Firstly, you are only interested in stars in the CCSN mass range
    m_max = MmaxSnII;
  
  double Y; 
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
  return result / TotalMassSN * Y;
} 

double PISN_PopIII_Yield(int yield_type) //Yield_type = 0 -> Recycling, 1 -> Metals
{
  //Remember that PISN feedback is always contemporaneous and they leave no remnants!
  double MassPISN = run_globals.MassPISN;
  double MassSNII = run_globals.MassSNII;
  double PISN_MassFraction = MassPISN / MassPISN; 
  if (yield_type == 0) 
    return PISN_MassFraction; 
  if (yield_type == 1)
    return PISN_MassFraction / 2.0;      
}

double NLBias(double Dist_Radius, double Halo_Mass, double redshift) 
// Fitting function written to match the results of Iliev+03 and Dijkstra+08. 
// Input in internal units. We don't use this in the small box where we have a grid resolution of about 200ckpc
// We want to use it for smaller grid resolution (about 1Mpc)
{
  double Alpha_ind = run_globals.params.physics.AlphaCluster;
  double Beta_ind = run_globals.params.physics.BetaCluster;
  double Gamma_ind = run_globals.params.physics.GammaCluster;
  double Psi_Norm = run_globals.params.physics.NormCluster;
  
  double little_h = run_globals.params.Hubble_h;
  Dist_Radius /= little_h;  
  
  Halo_Mass = Halo_Mass * 1e10 / little_h;
  
  return (Psi_Norm * pow(Dist_Radius / 0.01, Alpha_ind) * pow(Halo_Mass / 1e6, Beta_ind) * pow(redshift / 20.0, Gamma_ind));
}
