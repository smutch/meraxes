/*
 * This code is an amalgamation of the requisite functions for X-ray heating taken from
 * 21cmFAST. Specifically, from heating_helper_progs.c and elec_interp.c.
 * The original code was written by Andrei Mesinger.
 * Written to work within Meraxes by Bradley Greig.
 */

#include <assert.h>
#include <fftw3-mpi.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

// DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>

#include "reionization.h"

// TODO: Remove global variables here!
#define _XRAY_HEATING_FUNCTIONS_C
#include "XRayHeatingFunctions.h"

static double *sigma_atR, *sigma_Tmin, *ST_over_PS;
static int NO_LIGHT;

static float x_int_Energy[x_int_NENERGY];
static float x_int_fheat[x_int_NXHII][x_int_NENERGY];
static float x_int_n_Lya[x_int_NXHII][x_int_NENERGY];
static float x_int_nion_HI[x_int_NXHII][x_int_NENERGY];
static float x_int_nion_HeI[x_int_NXHII][x_int_NENERGY];
static float x_int_nion_HeII[x_int_NXHII][x_int_NENERGY];

int init_heat()
{

  size_t TsNumFilterSteps = (size_t)run_globals.params.TsNumFilterSteps;
  zpp_edge = calloc(TsNumFilterSteps, sizeof(double)); // calloc allocates the requested memory and returns a pointer to
                                                       // it. callocs(number of elem to be allocated, size)
  sigma_atR = calloc(TsNumFilterSteps, sizeof(double));
  sigma_Tmin = calloc(TsNumFilterSteps, sizeof(double));
  ST_over_PS = calloc(TsNumFilterSteps, sizeof(double));
  sum_lyn = calloc(TsNumFilterSteps, sizeof(double));
  if (run_globals.params.Flag_IncludeLymanWerner)
    sum_lyn_LW = calloc(TsNumFilterSteps, sizeof(double));

  kappa_10(1.0, 1); // 1 is the flag, allocates memory.
  if (kappa_10_elec(1.0, 1) < 0)
    return -2;
  if (kappa_10_pH(1.0, 1) < 0)
    return -3;
  if (T_RECFAST(100, 1) < 0)
    return -4;
  if (xion_RECFAST(100, 1) < 0)
    return -5;
  if (spectral_emissivity(0, 1) < 0)
    return -6;

  initialize_interp_arrays();

  return 0;
}

void destruct_heat()
{
  spectral_emissivity(0.0, 2); // 2 is the flag, frees memory.
  xion_RECFAST(100.0, 2);
  T_RECFAST(100.0, 2);
  kappa_10_pH(1.0, 2);
  kappa_10_elec(1.0, 2);
  kappa_10(1.0, 2);

  free(sum_lyn);
  free(ST_over_PS);
  free(sigma_Tmin);
  free(sigma_atR);
  free(zpp_edge);
  if (run_globals.params.Flag_IncludeLymanWerner)
    free(sum_lyn_LW);
}

// ******************************************************************** //
//  ************************ RECFAST quantities ************************ //
//  ******************************************************************** //

// * IGM temperature from RECFAST; includes Compton heating and adiabatic expansion only. * //
double T_RECFAST(float z, int flag)
{
  double ans;
  static double zt[RECFAST_NPTS], TK[RECFAST_NPTS];
  static gsl_interp_accel* acc;
  static gsl_spline* spline;
  float currz, currTK, trash;
  int i;
  FILE* F;

  char fname[STRLEN];

  if (flag == 1) {

    if (run_globals.mpi_rank == 0) {

      sprintf(fname, "%s/recfast_LCDM.dat", run_globals.params.TablesForXHeatingDir);

      // Read in the data
      if (!(F = fopen(fname, "r"))) {
        mlog("T_RECFAST: Unable to open file: %s for reading\nAborting\n", MLOG_MESG, fname);
        return -1;
      }

      for (i = (RECFAST_NPTS - 1); i >= 0; i--) {
        fscanf(F, "%f %E %E %E", &currz, &trash, &trash, &currTK);
        zt[i] = currz;
        TK[i] = currTK;
      }
      fclose(F);
    }

    // broadcast the values to all cores
    MPI_Bcast(zt, sizeof(zt), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(TK, sizeof(TK), MPI_BYTE, 0, run_globals.mpi_comm);

    // Set up spline table
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, RECFAST_NPTS);
    gsl_spline_init(spline, zt, TK, RECFAST_NPTS);

    return 0;
  }

  if (flag == 2) {
    // Free memory
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  if (z > zt[RECFAST_NPTS - 1]) { // Called at z>500! Bail out
    mlog("Called T_RECFAST with z=%f, bailing out!\n", MLOG_MESG, z);
    return -1;
  } else { // Do spline
    ans = gsl_spline_eval(spline, z, acc);
  }
  return ans;
}

// * Ionization fraction from RECFAST. * //
double xion_RECFAST(float z, int flag)
{
  static double zt[RECFAST_NPTS], xion[RECFAST_NPTS];
  static gsl_interp_accel* acc;
  static gsl_spline* spline;
  float trash, currz, currxion;
  double ans;
  int i;
  FILE* F;

  char fname[STRLEN];

  if (flag == 1) {

    if (run_globals.mpi_rank == 0) {

      sprintf(fname, "%s/recfast_LCDM.dat", run_globals.params.TablesForXHeatingDir);

      // Read in the data
      if (!(F = fopen(fname, "r"))) {
        mlog("xion_RECFAST: Unable to open file: %s for reading\nAborting\n", MLOG_MESG, fname);
        return -1;
      }

      for (i = (RECFAST_NPTS - 1); i >= 0; i--) {
        fscanf(F, "%f %E %E %E", &currz, &currxion, &trash, &trash);
        zt[i] = currz;
        xion[i] = currxion;
      }
      fclose(F);
    }

    // broadcast the values to all cores
    MPI_Bcast(zt, sizeof(zt), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(xion, sizeof(xion), MPI_BYTE, 0, run_globals.mpi_comm);

    // Set up spline table
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, RECFAST_NPTS);
    gsl_spline_init(spline, zt, xion, RECFAST_NPTS);

    return 0;
  }

  if (flag == 2) {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  if (z > zt[RECFAST_NPTS - 1]) { // Called at z>500! Bail out
    mlog("Called xion_RECFAST with z=%f, bailing out!\n", MLOG_MESG, z);
    return -1;
  } else { // Do spline
    ans = gsl_spline_eval(spline, z, acc);
  }
  return ans;
}

//  Returns the frequency threshold where \tau_X = 1, given parameter values of
//  electron fraction in the IGM outside of HII regions, x_e,
//  recieved redshift, zp, and emitted redshift, zpp.
typedef struct
{
  double x_e, zp, zpp, fcoll, HI_filling_factor_zp;
  int snap_i;
} nu_tau_one_params;
double nu_tau_one_helper(double nu, void* params)
{
  nu_tau_one_params* p = (nu_tau_one_params*)params;
  return tauX(nu, p->x_e, p->zp, p->zpp, p->fcoll, p->HI_filling_factor_zp, p->snap_i) - 1;
}
double nu_tau_one(double zp, double zpp, double x_e, double fcoll, double HI_filling_factor_zp, int snap_i)
{
  int status, iter, max_iter;
  const gsl_root_fsolver_type* T;
  gsl_root_fsolver* s;
  gsl_function F;
  double x_lo, x_hi, r = 0;
  double relative_error = 0.02;
  nu_tau_one_params p;

  // check if too ionized
  if (x_e > 0.9999) {
    mlog("WARNING: x_e value is too close to 1 for convergence in nu_tau_one\n", MLOG_MESG);
    return -1;
  }

  // select solver and allocate memory
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T); // non-derivative based Brent method
  if (!s) {
    mlog("Unable to allocate memory in function nu_tau_one\n", MLOG_MESG);
    return -1;
  }

  // check if lower bound has null
  if (tauX(HeI_NUIONIZATION, x_e, zp, zpp, fcoll, HI_filling_factor_zp, snap_i) < 1)
    return HeI_NUIONIZATION;

  // set frequency boundary values
  x_lo = HeI_NUIONIZATION;
  x_hi = 1e6 * NU_over_EV;

  // select function we wish to solve
  p.x_e = x_e;
  p.zp = zp;
  p.zpp = zpp;
  p.fcoll = fcoll;
  p.HI_filling_factor_zp = HI_filling_factor_zp;
  p.snap_i = snap_i;
  F.function = &nu_tau_one_helper;
  F.params = &p;
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  // iterate until we guess close enough
  iter = 0;
  max_iter = 100;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, relative_error);
  }

  while (status == GSL_CONTINUE && iter < max_iter);

  // deallocate and return
  gsl_root_fsolver_free(s);

  return r;
}

//  Calculates the optical depth for a photon arriving at z = zp with frequency nu, emitted at z = zpp.
//  The filling factor of neutral IGM at zp is HI_filling_factor_zp.
typedef struct
{
  double nu_0, x_e, ion_eff;
  int snap_i;
} tauX_params;
double tauX_integrand(double zhat, void* params)
{
  double n, drpropdz, nuhat, HI_filling_factor_zhat, sigma_tilde, fcoll;
  tauX_params* p = (tauX_params*)params;

  drpropdz = SPEED_OF_LIGHT * dtdz((float)zhat);
  n = N_b0 * pow(1 + zhat, 3);
  nuhat = p->nu_0 * (1 + zhat);
  fcoll = interpolate_fcoll(zhat, p->snap_i);
  if (fcoll < 1e-20)
    HI_filling_factor_zhat = 1;
  else
    HI_filling_factor_zhat =
      1 - p->ion_eff * fcoll /
            (1.0 - x_e_ave); // simplification to use the <x_e> value at zp and not zhat.  should'nt matter much since
                             // the evolution in x_e_ave is slower than fcoll.  in principle should make an ar$
  if (HI_filling_factor_zhat < 1e-4)
    HI_filling_factor_zhat = 1e-4; // set a floor for post-reionization stability

  sigma_tilde = species_weighted_x_ray_cross_section(nuhat, p->x_e);
  return drpropdz * n * HI_filling_factor_zhat * sigma_tilde;
}
double tauX(double nu, double x_e, double zp, double zpp, double fcoll, double HI_filling_factor_zp, int snap_i)
{
  double result, error;
  gsl_function F;
  double rel_tol = 0.005; //<- relative tolerance
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
  tauX_params p;

  F.function = &tauX_integrand;
  p.nu_0 = nu / (1 + zp);
  p.x_e = x_e;
  // effective efficiency for the PS (not ST) mass function; quicker to compute...
  if (HI_filling_factor_zp > FRACT_FLOAT_ERR) {
    p.ion_eff = run_globals.params.physics.ReionEfficiency;
  } else
    p.ion_eff = run_globals.params.physics.ReionEfficiency;

  p.snap_i = snap_i;

  F.params = &p;
  gsl_integration_qag(&F, zpp, zp, 0, rel_tol, 1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

double dtdz(float z)
{
  double x, dxdz, const1, denom, numer, OMl, OMm;

  OMl = run_globals.params.OmegaLambda;
  OMm = run_globals.params.OmegaM;

  x = sqrt(OMl / OMm) * pow(1 + z, -3.0 / 2.0);
  dxdz = sqrt(OMl / OMm) * pow(1 + z, -5.0 / 2.0) * (-3.0 / 2.0);
  const1 = 2 * sqrt(1 + OMm / OMl) / (3.0 * HUBBLE * run_globals.params.Hubble_h);

  numer = dxdz * (1 + x * pow(pow(x, 2) + 1, -0.5));
  denom = x + sqrt(pow(x, 2) + 1);
  return (const1 * numer / denom);
}

// returns the hubble "constant" (in 1/sec) at z //
double hubble(float z)
{

  double OMr, OMl, OMm;

  OMl = run_globals.params.OmegaLambda;
  OMm = run_globals.params.OmegaM;
  OMr = run_globals.params.OmegaR;

  return HUBBLE * run_globals.params.Hubble_h * sqrt(OMm * pow(1 + z, 3) + OMr * pow(1 + z, 4) + OMl);
}

//  The total weighted HI + HeI + HeII  cross-section in pcm^-2
//  technically, the x_e should be local, line of sight (not global) here,
//  but that would be very slow...
double species_weighted_x_ray_cross_section(double nu, double x_e)
{
  double HI_factor, HeI_factor, HeII_factor;

  HI_factor = f_H * (1 - x_e) * HI_ion_crosssec(nu);
  HeI_factor = f_He * (1 - x_e) * HeI_ion_crosssec(nu);
  HeII_factor = f_He * x_e * HeII_ion_crosssec(nu);

  return HI_factor + HeI_factor + HeII_factor;
}

// function HeI_ion_crosssec returns the HI ionization cross section at parameter frequency (taken from Verner et al
// (1996)
double HeI_ion_crosssec(double nu)
{
  double x, y;

  if (nu < HeI_NUIONIZATION)
    return 0;

  x = nu / NU_over_EV / 13.61 - 0.4434;
  y = sqrt(x * x + pow(2.136, 2));
  return 9.492e-16 * ((x - 1) * (x - 1) + 2.039 * 2.039) * pow(y, (0.5 * 3.188 - 5.5)) *
         pow(1.0 + sqrt(y / 1.469), -3.188);
}

// function HeII_ion_crosssec returns the HeII ionization cross section at parameter frequency (taken from Osterbrock,
// pg. 14)
double HeII_ion_crosssec(double nu)
{
  double epsilon, Z = 2;

  if (nu < HeII_NUIONIZATION)
    return 0;

  if (nu == HeII_NUIONIZATION)
    nu += TINY;

  epsilon = sqrt(nu / HeII_NUIONIZATION - 1);
  return (6.3e-18) / Z / Z * pow(HeII_NUIONIZATION / nu, 4) * exp(4 - (4 * atan(epsilon) / epsilon)) /
         (1 - exp(-2 * M_PI / epsilon));
}

// function HI_ion_crosssec returns the HI ionization cross section at parameter frequency (taken from Osterbrock, pg.
// 14)
double HI_ion_crosssec(double nu)
{
  double epsilon, Z = 1;

  if (nu < NUIONIZATION)
    return 0;

  if (nu == NUIONIZATION)
    nu += TINY;

  epsilon = sqrt(nu / NUIONIZATION - 1);
  return SIGMA_HI / Z / Z * pow(NUIONIZATION / nu, 4) * exp(4 - (4 * atan(epsilon) / epsilon)) /
         (1 - exp(-2 * M_PI / epsilon));
}

// comoving distance (in cm) per unit redshift
double drdz(float z)
{
  return (1.0 + z) * SPEED_OF_LIGHT * dtdz(z);
}

/* function INVSINH returns the inverse hyperbolic sine of parameter x */
double invsinh(double x)
{
  return log(x + sqrt(pow(x, 2) + 1));
}

/* function GETTIME returns the age of the universe, at a given redshift parameter, z.
   This assumes zero curvature.  (from Weinberg 1989) */
double gettime(double z)
{
  double term1, const1, OMl, OMm;

  OMl = run_globals.params.OmegaLambda;
  OMm = run_globals.params.OmegaM;

  term1 = invsinh(sqrt(OMl / OMm) * pow(1 + z, -3.0 / 2.0));
  const1 = 2 * sqrt(1 + OMm / OMl) / (3 * HUBBLE * run_globals.params.Hubble_h);
  return (term1 * const1);
}

// Returns the maximum redshift at which a Lyn transition contributes to Lya flux at z
float zmax(float z, int n)
{
  double num, denom;
  num = 1 - pow(n + 1, -2);
  denom = 1 - pow(n, -2);
  return (float)((1 + z) * num / denom - 1);
}

// Returns frequency of Lyman-n, in units of Lyman-alpha
double nu_n(int n)
{
  double ans;

  ans = 1.0 - pow(n, -2.0);
  ans /= 0.75;
  return ans;
}

// Returns recycling fraction (=fraction of photons converted into Lyalpha for Ly-n resonance
double frecycle(int n)
{
  switch (n) {
    case 0:
      return 1;
    case 1:
      return 1;
    case 2:
      return 1;
    case 3:
      return 0;
    case 4:
      return 0.2609;
    case 5:
      return 0.3078;
    case 6:
      return 0.3259;
    case 7:
      return 0.3353;
    case 8:
      return 0.3410;
    case 9:
      return 0.3448;
    case 10:
      return 0.3476;
    case 11:
      return 0.3496;
    case 12:
      return 0.3512;
    case 13:
      return 0.3524;
    case 14:
      return 0.3535;
    case 15:
      return 0.3543;
    case 16:
      return 0.3550;
    case 17:
      return 0.3556;
    case 18:
      return 0.3561;
    case 19:
      return 0.3565;
    case 20:
      return 0.3569;
    case 21:
      return 0.3572;
    case 22:
      return 0.3575;
    case 23:
      return 0.3578;
    case 24:
      return 0.3580;
    case 25:
      return 0.3582;
    case 26:
      return 0.3584;
    case 27:
      return 0.3586;
    case 28:
      return 0.3587;
    case 29:
      return 0.3589;
    case 30:
      return 0.3590;
    default:
      return 0;
  }
}

// Reads in and constructs table of the piecewise power-law fits to Pop 2 and Pop 3 stellar spectra, from Barkana
double spectral_emissivity(double nu_norm, int flag)
{
  static int n[NSPEC_MAX];
  static float nu_n[NSPEC_MAX], alpha_S_2[NSPEC_MAX];
  static float alpha_S_3[NSPEC_MAX], N0_2[NSPEC_MAX], N0_3[NSPEC_MAX];
  double n0_fac;
  double ans = 0;
  int i;
  FILE* F;

  char fname[STRLEN];

  switch (flag) {
    case 2:
      for (i = 1; i < (NSPEC_MAX - 1); i++) {
        if ((nu_norm >= nu_n[i]) && (nu_norm < nu_n[i + 1])) {
          if (Pop == 2) {
            ans = N0_2[i] / (alpha_S_2[i] + 1) * (pow(nu_n[i + 1], alpha_S_2[i] + 1) - pow(nu_norm, alpha_S_2[i] + 1));
          } else if (Pop == 3) {
            ans = N0_3[i] / (alpha_S_3[i] + 1) * (pow(nu_n[i + 1], alpha_S_3[i] + 1) - pow(nu_norm, alpha_S_3[i] + 1));
          } else {
            mlog("Invalid value for Stellar Population", MLOG_MESG);
          }
          return ans > 0 ? ans : 1e-40;
        }
      }

    case 1:
      if (run_globals.mpi_rank == 0) {

        sprintf(fname, "%s/stellar_spectra.dat", run_globals.params.TablesForXHeatingDir);

        // Read in the data
        if (!(F = fopen(fname, "r"))) {
          mlog("spectral_emissivity: Unable to open file: stellar_spectra.dat at %s for reading\nAborting\n",
               MLOG_MESG,
               fname);
          return -1;
        }

        for (i = 1; i < NSPEC_MAX; i++) {
          fscanf(F, "%i %e %e %e %e", &n[i], &N0_2[i], &alpha_S_2[i], &N0_3[i], &alpha_S_3[i]);
        }
        fclose(F);

        for (i = 1; i < NSPEC_MAX; i++) {
          nu_n[i] = (float)(4.0 / 3.0 * (1.0 - 1.0 / pow(n[i], 2.0)));
        }

        for (i = 1; i < (NSPEC_MAX - 1); i++) {
          n0_fac = (pow(nu_n[i + 1], alpha_S_2[i] + 1) - pow(nu_n[i], alpha_S_2[i] + 1));
          N0_2[i] *= (alpha_S_2[i] + 1) / n0_fac * Pop2_ion;
          n0_fac = (pow(nu_n[i + 1], alpha_S_3[i] + 1) - pow(nu_n[i], alpha_S_3[i] + 1));
          N0_3[i] *= (alpha_S_3[i] + 1) / n0_fac * Pop3_ion;
        }
      }

      // broadcast the values to all cores
      MPI_Bcast(nu_n, sizeof(nu_n), MPI_BYTE, 0, run_globals.mpi_comm);
      MPI_Bcast(alpha_S_2, sizeof(alpha_S_2), MPI_BYTE, 0, run_globals.mpi_comm);
      MPI_Bcast(alpha_S_3, sizeof(alpha_S_3), MPI_BYTE, 0, run_globals.mpi_comm);
      MPI_Bcast(N0_2, sizeof(N0_2), MPI_BYTE, 0, run_globals.mpi_comm);
      MPI_Bcast(N0_3, sizeof(N0_3), MPI_BYTE, 0, run_globals.mpi_comm);

      return 0.0;

    default:
      for (i = 1; i < (NSPEC_MAX - 1); i++) {
        if ((nu_norm >= nu_n[i]) && (nu_norm < nu_n[i + 1])) {
          // We are in the correct spectral region
          if (Pop == 2)
            ans = N0_2[i] * pow(nu_norm, alpha_S_2[i]);
          else
            ans = N0_3[i] * pow(nu_norm, alpha_S_3[i]);

          return ans / Ly_alpha_HZ;
        }
      }

      i = NSPEC_MAX - 1;
      if (Pop == 2)
        return N0_2[i] * pow(nu_norm, alpha_S_2[i]) / Ly_alpha_HZ;
      else
        return N0_3[i] * pow(nu_norm, alpha_S_3[i]) / Ly_alpha_HZ;
  }
}

typedef struct
{
  double x_e, NU_X_THRESH, X_RAY_SPEC_INDEX;
} int_over_nu_params;

//  Evaluates the frequency integral in the Tx evolution equation
//  photons starting from zpp arive at zp, with mean IGM electron
//  fraction of x_e (used to compute tau), and local electron
//  fraction local_x_e
//  FLAG = 0 for heat integral
//  FLAG = 1 for ionization integral
//  FLAG = 2 for Lya integral
double integrand_in_nu_heat_integral(double nu, void* params)
{

  double species_sum;

  int_over_nu_params* p = (int_over_nu_params*)params;

  // HI
  species_sum = interp_fheat((float)((nu - NUIONIZATION) / NU_over_EV), (float)p->x_e) * PLANCK * (nu - NUIONIZATION) *
                f_H * (1 - p->x_e) * HI_ion_crosssec(nu);

  // HeI
  species_sum += interp_fheat((float)((nu - HeI_NUIONIZATION) / NU_over_EV), (float)p->x_e) * PLANCK *
                 (nu - HeI_NUIONIZATION) * f_He * (1 - p->x_e) * HeI_ion_crosssec(nu);

  // HeII
  species_sum += interp_fheat((float)((nu - HeII_NUIONIZATION) / NU_over_EV), (float)p->x_e) * PLANCK *
                 (nu - HeII_NUIONIZATION) * f_He * p->x_e * HeII_ion_crosssec(nu);

  return species_sum * pow(nu / (p->NU_X_THRESH * NU_over_EV), -p->X_RAY_SPEC_INDEX - 1);
}

double integrand_in_nu_ion_integral(double nu, void* params)
{
  double species_sum, F_i;

  int_over_nu_params* p = (int_over_nu_params*)params;

  // photoionization of HI, prodicing e- of energy h*(nu - nu_HI)
  F_i = interp_nion_HI((float)((nu - NUIONIZATION) / NU_over_EV), (float)p->x_e) +
        interp_nion_HeI((float)((nu - NUIONIZATION) / NU_over_EV), (float)p->x_e) +
        interp_nion_HeII((float)((nu - NUIONIZATION) / NU_over_EV), (float)p->x_e) + 1;
  species_sum = F_i * f_H * (1 - p->x_e) * HI_ion_crosssec(nu);

  // photoionization of HeI, prodicing e- of energy h*(nu - nu_HeI)
  F_i = interp_nion_HI((float)((nu - HeI_NUIONIZATION) / NU_over_EV), (float)p->x_e) +
        interp_nion_HeI((float)((nu - HeI_NUIONIZATION) / NU_over_EV), (float)p->x_e) +
        interp_nion_HeII((float)((nu - HeI_NUIONIZATION) / NU_over_EV), (float)p->x_e) + 1;
  species_sum += F_i * f_He * (1 - p->x_e) * HeI_ion_crosssec(nu);

  // photoionization of HeII, prodicing e- of energy h*(nu - nu_HeII)
  F_i = interp_nion_HI((float)((nu - HeII_NUIONIZATION) / NU_over_EV), (float)p->x_e) +
        interp_nion_HeI((float)((nu - HeII_NUIONIZATION) / NU_over_EV), (float)p->x_e) +
        interp_nion_HeII((float)((nu - HeII_NUIONIZATION) / NU_over_EV), (float)p->x_e) + 1;
  species_sum += F_i * f_He * p->x_e * HeII_ion_crosssec(nu);

  return species_sum * pow(nu / (p->NU_X_THRESH * NU_over_EV), -p->X_RAY_SPEC_INDEX - 1);
}

double integrand_in_nu_lya_integral(double nu, void* params)
{
  double species_sum;

  int_over_nu_params* p = (int_over_nu_params*)params;

  // HI
  species_sum =
    interp_n_Lya((float)((nu - NUIONIZATION) / NU_over_EV), (float)p->x_e) * f_H * (1 - p->x_e) * HI_ion_crosssec(nu);

  // HeI
  species_sum += interp_n_Lya((float)((nu - HeI_NUIONIZATION) / NU_over_EV), (float)p->x_e) * f_He * (1 - p->x_e) *
                 HeI_ion_crosssec(nu);

  // HeII
  species_sum +=
    interp_n_Lya((float)((nu - HeII_NUIONIZATION) / NU_over_EV), (float)p->x_e) * f_He * p->x_e * HeII_ion_crosssec(nu);

  return species_sum * pow(nu / (p->NU_X_THRESH * NU_over_EV), -p->X_RAY_SPEC_INDEX - 1);
}

double integrate_over_nu(double zp,
                         double local_x_e,
                         double lower_int_limit,
                         double thresh_energy,
                         double spec_index,
                         int FLAG)
{
  double result, error;
  double rel_tol = 0.01; //<- relative tolerance
  gsl_function F;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  int_over_nu_params p;

  p.x_e = local_x_e;
  p.NU_X_THRESH = thresh_energy;
  p.X_RAY_SPEC_INDEX = spec_index;

  F.params = &p;

  if (FLAG == 0)
    F.function = &integrand_in_nu_heat_integral;
  else if (FLAG == 1)
    F.function = &integrand_in_nu_ion_integral;
  else {
    F.function = &integrand_in_nu_lya_integral;
  }

  gsl_integration_qag(&F,
                      lower_int_limit,
                      run_globals.params.physics.NuXrayMax * NU_over_EV,
                      0,
                      rel_tol,
                      1000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &error);
  gsl_integration_workspace_free(w);

  // if it is the Lya integral, add prefactor
  if (FLAG == 2)
    return result * SPEED_OF_LIGHT / (4 * M_PI) / Ly_alpha_HZ / hubble((float)zp);

  return result;
}

// Call once to read in data files and set up arrays for interpolation.
// All data files should be in a local subdirectory "x_int_tables/"; if moved, change input_base
// below to new location.
void initialize_interp_arrays()
{
  FILE* input_file;
  char input_file_name[500];

  char input_base[] = "x_int_tables/";
  char input_tail[100] = ".dat";
  char mode[10] = "r";

  float xHI, xHeI, xHeII, z, T;
  float trash;
  char label[64];

  int i;
  int n_ion;

  // Initialize array of ionized fractions
  x_int_XHII[0] = 1.0e-4;
  x_int_XHII[1] = 2.318e-4;
  x_int_XHII[2] = 4.677e-4;
  x_int_XHII[3] = 1.0e-3;
  x_int_XHII[4] = 2.318e-3;
  x_int_XHII[5] = 4.677e-3;
  x_int_XHII[6] = 1.0e-2;
  x_int_XHII[7] = 2.318e-2;
  x_int_XHII[8] = 4.677e-2;
  x_int_XHII[9] = 1.0e-1;
  x_int_XHII[10] = 0.5;
  x_int_XHII[11] = 0.9;
  x_int_XHII[12] = 0.99;
  x_int_XHII[13] = 0.999;

  if (run_globals.mpi_rank == 0) {

    for (n_ion = 0; n_ion < x_int_NXHII; n_ion++) {

      // Construct filename
      if (x_int_XHII[n_ion] < 0.3) {
        sprintf(input_file_name,
                "%s/%slog_xi_%1.1f%s",
                run_globals.params.TablesForXHeatingDir,
                input_base,
                log10(x_int_XHII[n_ion]),
                input_tail);
      } else {
        sprintf(input_file_name,
                "%s/%sxi_%1.3f%s",
                run_globals.params.TablesForXHeatingDir,
                input_base,
                x_int_XHII[n_ion],
                input_tail);
      }

      input_file = fopen(input_file_name, mode);

      if (input_file == NULL) {
        mlog("Can't open input file %s!\n", MLOG_MESG, input_file_name);
        exit(1);
      }

      // Read in first line
      for (i = 1; i <= 5; i++) {
        fscanf(input_file, "%s", label);
      }

      // Read in second line (ionized fractions info)
      fscanf(input_file, "%g %g %g %g %g", &xHI, &xHeI, &xHeII, &z, &T);

      // Read in column headings
      for (i = 1; i <= 11; i++) {
        fscanf(input_file, "%s", label);
      }

      // Read in data table
      for (i = 0; i < x_int_NENERGY; i++) {
        fscanf(input_file,
               "%g %g %g %g %g %g %g %g %g",
               &x_int_Energy[i],
               &trash,
               &x_int_fheat[n_ion][i],
               &trash,
               &x_int_n_Lya[n_ion][i],
               &x_int_nion_HI[n_ion][i],
               &x_int_nion_HeI[n_ion][i],
               &x_int_nion_HeII[n_ion][i],
               &trash);
      }

      fclose(input_file);
    }
  }

  // broadcast the values to all cores
  MPI_Bcast(&x_int_Energy, sizeof(x_int_Energy), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_fheat, sizeof(x_int_fheat), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_n_Lya, sizeof(x_int_n_Lya), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_nion_HI, sizeof(x_int_nion_HI), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_nion_HeI, sizeof(x_int_nion_HeI), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_nion_HeII, sizeof(x_int_nion_HeII), MPI_BYTE, 0, run_globals.mpi_comm);
}

// Function to compute fheat for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the fraction of the this electron energy that is returned to heat, thus
// it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that highest value.  Similarly, if
// xHII_call is less than or smaller than the limits of the ionized fraction array (10^-4 and
// 0.999), it just uses those values.
float interp_fheat(float En, float xHII_call)
{

  int n_low, n_high;
  int m_xHII_low, m_xHII_high;

  float elow_result, ehigh_result, final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999 * x_int_Energy[x_int_NENERGY - 1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = (float)(x_int_Energy[x_int_NENERGY - 1] * 0.999);
  } else if (En < x_int_Energy[0]) {
    return 1.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII - 1] * 0.999) {
    xHII_call = (float)(x_int_XHII[x_int_NXHII - 1] * 0.999);
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = (float)(1.001 * x_int_XHII[0]);
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;

  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;

  // First linear interpolation in energy
  elow_result =
    ((x_int_fheat[m_xHII_low][n_high] - x_int_fheat[m_xHII_low][n_low]) / (x_int_Energy[n_high] - x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_fheat[m_xHII_low][n_low];

  // Second linear interpolation in energy
  ehigh_result = ((x_int_fheat[m_xHII_high][n_high] - x_int_fheat[m_xHII_high][n_low]) /
                  (x_int_Energy[n_high] - x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_fheat[m_xHII_high][n_low];

  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;

  return final_result;
}

// Function to compute nLya for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of Lyman-alpha photons produced as a function of the secondary
// electron energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_n_Lya(float En, float xHII_call)
{
  int n_low, n_high;
  int m_xHII_low, m_xHII_high;

  float elow_result, ehigh_result, final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999 * x_int_Energy[x_int_NENERGY - 1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = (float)(x_int_Energy[x_int_NENERGY - 1] * 0.999);
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII - 1] * 0.999) {
    xHII_call = (float)(x_int_XHII[x_int_NXHII - 1] * 0.999);
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = (float)(1.001 * x_int_XHII[0]);
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;

  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;

  // First linear interpolation in energy
  elow_result =
    ((x_int_n_Lya[m_xHII_low][n_high] - x_int_n_Lya[m_xHII_low][n_low]) / (x_int_Energy[n_high] - x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_n_Lya[m_xHII_low][n_low];

  // Second linear interpolation in energy
  ehigh_result = ((x_int_n_Lya[m_xHII_high][n_high] - x_int_n_Lya[m_xHII_high][n_low]) /
                  (x_int_Energy[n_high] - x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_n_Lya[m_xHII_high][n_low];

  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;

  return final_result;
}

// Function to compute nHI for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of hydrogen ionizations produced as a function of the secondary
// electron energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_nion_HI(float En, float xHII_call)
{
  int n_low, n_high;
  int m_xHII_low, m_xHII_high;

  float elow_result, ehigh_result, final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999 * x_int_Energy[x_int_NENERGY - 1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = (float)(x_int_Energy[x_int_NENERGY - 1] * 0.999);
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII - 1] * 0.999) {
    xHII_call = (float)(x_int_XHII[x_int_NXHII - 1] * 0.999);
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = (float)(1.001 * x_int_XHII[0]);
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;

  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;

  // First linear interpolation in energy
  elow_result = ((x_int_nion_HI[m_xHII_low][n_high] - x_int_nion_HI[m_xHII_low][n_low]) /
                 (x_int_Energy[n_high] - x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_nion_HI[m_xHII_low][n_low];

  // Second linear interpolation in energy
  ehigh_result = ((x_int_nion_HI[m_xHII_high][n_high] - x_int_nion_HI[m_xHII_high][n_low]) /
                  (x_int_Energy[n_high] - x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_nion_HI[m_xHII_high][n_low];

  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;

  return final_result;
}

// Function to compute nHeI for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of HeI ionizations produced as a function of the secondary electron
// energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_nion_HeI(float En, float xHII_call)
{
  int n_low, n_high;
  int m_xHII_low, m_xHII_high;

  float elow_result, ehigh_result, final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999 * x_int_Energy[x_int_NENERGY - 1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = (float)(x_int_Energy[x_int_NENERGY - 1] * 0.999);
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII - 1] * 0.999) {
    xHII_call = (float)(x_int_XHII[x_int_NXHII - 1] * 0.999);
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = (float)(1.001 * x_int_XHII[0]);
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;

  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;

  // First linear interpolation in energy
  elow_result = ((x_int_nion_HeI[m_xHII_low][n_high] - x_int_nion_HeI[m_xHII_low][n_low]) /
                 (x_int_Energy[n_high] - x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_nion_HeI[m_xHII_low][n_low];

  // Second linear interpolation in energy
  ehigh_result = ((x_int_nion_HeI[m_xHII_high][n_high] - x_int_nion_HeI[m_xHII_high][n_low]) /
                  (x_int_Energy[n_high] - x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_nion_HeI[m_xHII_high][n_low];

  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;

  return final_result;
}

// Function to compute nHeII for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of HeII ionizations produced as a function of the secondary
// electron energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_nion_HeII(float En, float xHII_call)
{
  int n_low, n_high;
  int m_xHII_low, m_xHII_high;

  float elow_result, ehigh_result, final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999 * x_int_Energy[x_int_NENERGY - 1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = (float)(x_int_Energy[x_int_NENERGY - 1] * 0.999);
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII - 1] * 0.999) {
    xHII_call = (float)(x_int_XHII[x_int_NXHII - 1] * 0.999);
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = (float)(1.001 * x_int_XHII[0]);
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;

  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;

  // First linear interpolation in energy
  elow_result = ((x_int_nion_HeII[m_xHII_low][n_high] - x_int_nion_HeII[m_xHII_low][n_low]) /
                 (x_int_Energy[n_high] - x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_nion_HeII[m_xHII_low][n_low];

  // Second linear interpolation in energy
  ehigh_result = ((x_int_nion_HeII[m_xHII_high][n_high] - x_int_nion_HeII[m_xHII_high][n_low]) /
                  (x_int_Energy[n_high] - x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_nion_HeII[m_xHII_high][n_low];

  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result) / (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;

  return final_result;
}

// Function to find bounding indices on the energy array, for an input energy En.
// Note it is done exactly for speed, since all the energy arrays have the same structure.
int locate_energy_index(float En)
{
  int n_low;

  // Find energy table location analytically!
  if (En < 1008.88) {
    n_low = (int)(log(En / 10.0) / 1.98026273e-2);
  } else {
    n_low = 232 + (int)(log(En / 1008.88) / 9.53101798e-2);
  }

  return n_low;
}

// Function to find bounding indices on the ionized fraction array, for an input fraction
// xHII_call.  This is done by comparing to each element, since there are only 14 elements.
int locate_xHII_index(float xHII_call)
{
  int m_xHII_low;

  // Determine relevant ionized fractions to interpolate; do iteratively because few elements
  m_xHII_low = x_int_NXHII - 1;
  while (xHII_call < x_int_XHII[m_xHII_low]) {
    m_xHII_low--;
  }
  return m_xHII_low;
}

// ********************************************************************
// ************************** IGM Evolution ***************************
//  This function creates the d/dz' integrands
// *********************************************************************
void evolveInt(float zp,
               float curr_delNL0,
               const double SFR_GAL[],
               const double SFR_QSO[],
               const double freq_int_heat_GAL[],
               const double freq_int_ion_GAL[],
               const double freq_int_lya_GAL[],
               const double freq_int_heat_QSO[],
               const double freq_int_ion_QSO[],
               const double freq_int_lya_QSO[],
               int COMPUTE_Ts,
               const double y[],
               double deriv[])
{

  double dadia_dzp, dcomp_dzp, dxheat_dt_GAL, dxion_source_dt_GAL, dxion_sink_dt;
  double dxheat_dt_QSO, dxion_source_dt_QSO, dxlya_dt_QSO, dstarlya_dt_QSO;
  double zpp, dzpp;
  int zpp_ct;
  double T, x_e, zpp_integrand_GAL, zpp_integrand_QSO;
  double dxe_dzp, n_b, dspec_dzp, dxheat_dzp, dxlya_dt_GAL, dstarlya_dt_GAL, dstarlyLW_dt_GAL;

  x_e = y[0];
  T = y[1];
  n_b = N_b0 * pow(1 + zp, 3) * (1 + curr_delNL0);

  // First, let's do the trapazoidal integration over zpp
  dxheat_dt_GAL = 0;
  dxion_source_dt_GAL = 0;
  dxlya_dt_GAL = 0;
  dstarlya_dt_GAL = 0;
  dstarlyLW_dt_GAL = 0;

  dxheat_dt_QSO = 0;
  dxion_source_dt_QSO = 0;
  dxlya_dt_QSO = 0;
  dstarlya_dt_QSO = 0;

  if (!NO_LIGHT) {
    for (zpp_ct = 0; zpp_ct < run_globals.params.TsNumFilterSteps;
         zpp_ct++) { // Define last redshift that is effective, zpp_edge is defined in init_heat!
      // set redshift of half annulus; dz'' is negative since we flipped limits of integral
      if (zpp_ct == 0) {
        zpp = (zpp_edge[0] + zp) * 0.5;
        dzpp = zp - zpp_edge[0];
      } else {
        zpp = (zpp_edge[zpp_ct] + zpp_edge[zpp_ct - 1]) * 0.5;
        dzpp = zpp_edge[zpp_ct - 1] - zpp_edge[zpp_ct];
      }

      // Use this when using the SFR provided by Meraxes
      // Units should be M_solar/s. Factor of (dt_dzp * dzpp) converts from per s to per z'
      zpp_integrand_GAL =
        SFR_GAL[zpp_ct] * pow(1 + zpp, -run_globals.params.physics.SpecIndexXrayGal); // where is it in the formulas?
      if (run_globals.params.Flag_SeparateQSOXrays) {
        zpp_integrand_QSO = SFR_QSO[zpp_ct] * pow(1 + zpp, -run_globals.params.physics.SpecIndexXrayQSO);
      }

      if (run_globals.params.Flag_SeparateQSOXrays) {

        dxheat_dt_GAL += dt_dzpp * dzpp * zpp_integrand_GAL * freq_int_heat_GAL[zpp_ct];
        dxheat_dt_QSO += dt_dzpp * dzpp * zpp_integrand_QSO * freq_int_heat_QSO[zpp_ct];

        dxion_source_dt_GAL += dt_dzpp * dzpp * zpp_integrand_GAL * freq_int_ion_GAL[zpp_ct];
        dxion_source_dt_QSO += dt_dzpp * dzpp * zpp_integrand_QSO * freq_int_ion_QSO[zpp_ct];

        dxlya_dt_GAL += dt_dzpp * dzpp * zpp_integrand_GAL * freq_int_lya_GAL[zpp_ct];
        dxlya_dt_QSO += dt_dzpp * dzpp * zpp_integrand_QSO * freq_int_lya_QSO[zpp_ct];

        // Use this when using the SFR provided by Meraxes
        // Units should be M_solar/s. Factor of (dt_dzp * dzpp) converts from per s to per z'
        dstarlya_dt_GAL += SFR_GAL[zpp_ct] * pow(1 + zp, 2) * (1 + zpp) * sum_lyn[zpp_ct] * dt_dzpp * dzpp;
        dstarlya_dt_QSO += SFR_QSO[zpp_ct] * pow(1 + zp, 2) * (1 + zpp) * sum_lyn[zpp_ct] * dt_dzpp * dzpp;

      } else {
        dxheat_dt_GAL += dt_dzpp * dzpp * zpp_integrand_GAL *
                         freq_int_heat_GAL[zpp_ct]; // Integral in frequency must be computed for each TsNumFilterSteps
        dxion_source_dt_GAL += dt_dzpp * dzpp * zpp_integrand_GAL * freq_int_ion_GAL[zpp_ct];
        dxlya_dt_GAL += dt_dzpp * dzpp * zpp_integrand_GAL * freq_int_lya_GAL[zpp_ct];

        // Use this when using the SFR provided by Meraxes
        // Units should be M_solar/s. Factor of (dt_dzp * dzpp) converts from per s to per z'
        dstarlya_dt_GAL += SFR_GAL[zpp_ct] * pow(1 + zp, 2) * (1 + zpp) * sum_lyn[zpp_ct] * dt_dzpp * dzpp;
      }
      if (run_globals.params.Flag_IncludeLymanWerner)
        dstarlyLW_dt_GAL += SFR_GAL[zpp_ct] * pow(1 + zp, 2) * (1 + zpp) * sum_lyn_LW[zpp_ct] * dt_dzpp * dzpp;
    }

    // After you finish the loop for each Radius, you add prefactors which are constants for the redshift (snapshot) and
    // defined in ComputeTs.c
    if (run_globals.params.Flag_SeparateQSOXrays) {
      dxheat_dt_GAL *= const_zp_prefactor_GAL;
      dxion_source_dt_GAL *= const_zp_prefactor_GAL;
      dxlya_dt_GAL *= const_zp_prefactor_GAL * n_b;

      // Use this when using the SFR provided by Meraxes
      // Units should be M_solar/s. Factor of (dt_dzp * dzpp) converts from per s to per z'
      // The division by Omb * RHOcrit arises from the differences between eq. 13 and eq. 22 in Mesinger et al. (2011),
      // accounting for the M_solar factor (SFR -> number)
      dstarlya_dt_GAL *= (SPEED_OF_LIGHT / (4. * M_PI)) / (PROTONMASS / SOLAR_MASS);

      dxheat_dt_QSO *= const_zp_prefactor_QSO;
      dxion_source_dt_QSO *= const_zp_prefactor_QSO;
      dxlya_dt_QSO *= const_zp_prefactor_QSO * n_b;

      dstarlya_dt_QSO *= (SPEED_OF_LIGHT / (4. * M_PI)) / (PROTONMASS / SOLAR_MASS);
    } else {
      dxheat_dt_GAL *= const_zp_prefactor_GAL;
      dxion_source_dt_GAL *= const_zp_prefactor_GAL;
      dxlya_dt_GAL *= const_zp_prefactor_GAL * n_b;

      dstarlya_dt_GAL *= (SPEED_OF_LIGHT / (4. * M_PI)) / (PROTONMASS / SOLAR_MASS);
    }
    if (run_globals.params.Flag_IncludeLymanWerner)
      dstarlyLW_dt_GAL *= (SPEED_OF_LIGHT / (4. * M_PI)) / (PROTONMASS / SOLAR_MASS);

  } // end NO_LIGHT if statement

  // **** Now we can solve the evolution equations  ***** //
  // *** First let's do dxe_dzp *** //

  dxion_sink_dt = alpha_A(T) * CLUMPING_FACTOR * x_e * x_e * f_H * n_b;
  dxe_dzp = dt_dzp * ((dxion_source_dt_GAL + dxion_source_dt_QSO) - dxion_sink_dt);
  deriv[0] = dxe_dzp;

  // *** Next, let's get the temperature components *** //
  // first, adiabatic term (3rd term in 11)
  dadia_dzp = 3 / (1.0 + zp);
  if (fabs(curr_delNL0) > FRACT_FLOAT_ERR) // add adiabatic heating/cooling from structure formation
    dadia_dzp += dgrowth_factor_dzp / (1.0 / curr_delNL0 + growth_factor_zp);
  dadia_dzp *= (2.0 / 3.0) * T;

  // next heating due to the changing species
  dspec_dzp = -dxe_dzp * T / (1 + x_e);

  // next, Compton heating
  dcomp_dzp = dT_comp(zp, T, x_e);

  // lastly, X-ray heating
  dxheat_dzp = (dxheat_dt_GAL + dxheat_dt_QSO) * dt_dzp * 2.0 / 3.0 / BOLTZMANN / (1.0 + x_e);

  // summing them up...
  deriv[1] = dxheat_dzp + dcomp_dzp + dspec_dzp + dadia_dzp;

  // *** Finally, if we are at the last redshift step, Lya *** //
  deriv[2] = (dxlya_dt_GAL + dxlya_dt_QSO) + (dstarlya_dt_GAL + dstarlya_dt_QSO);

  // stuff for marcos
  deriv[3] = dxheat_dzp;
  deriv[4] = dt_dzp * (dxion_source_dt_GAL + dxion_source_dt_QSO);
  if (run_globals.params.Flag_IncludeLymanWerner)
    deriv[5] = dstarlyLW_dt_GAL * (PLANCK * 1e21);
}

// * Compton heating term * //
double dT_comp(double z, double TK, double xe)
{
  double Trad, ans;

  Trad = TCMB * (1.0 + z);
  ans = (-1.51e-4) * (xe / (1.0 + xe)) / (hubble((float)z) / (HUBBLE * run_globals.params.Hubble_h)) /
        run_globals.params.Hubble_h * pow(Trad, 4.0) / (1.0 + z);
  ans *= Trad - TK;
  return ans;
}

// * returns the case A hydrogen recombination coefficient (Abel et al. 1997) in cm^3 s^-1 * //
double alpha_A(double T)
{
  double logT, ans;
  logT = log(T / (double)1.1604505e4);
  ans = exp(-28.6130338 - 0.72411256 * logT - 2.02604473e-2 * pow(logT, 2) - 2.38086188e-3 * pow(logT, 3) -
            3.21260521e-4 * pow(logT, 4) - 1.42150291e-5 * pow(logT, 5) + 4.98910892e-6 * pow(logT, 6) +
            5.75561414e-7 * pow(logT, 7) - 1.85676704e-8 * pow(logT, 8) - 3.07113524e-9 * pow(logT, 9));
  return ans;
}

//  FUNCTION dicke(z)
//  Computes the dicke growth function at redshift z, i.e. the z dependance part of sigma
//
//  References: Peebles, "Large-Scale...", pg.53 (eq. 11.16). Includes omega<=1
//  Nonzero Lambda case from Liddle et al, astro-ph/9512102, eqs. 6-8.
//  and quintessence case from Wang et al, astro-ph/9804015
//
//  Normalized to dicke(z=0)=1

double dicke(double z)
{
  double omegaM_z, dick_z, dick_0, x, x_0;
  double tiny = 1e-4;

  double OMm = run_globals.params.OmegaM;
  double OMl = run_globals.params.OmegaLambda;
  double OMr = run_globals.params.OmegaR;
  double OMk = run_globals.params.OmegaK;

  double wl = run_globals.params.wLambda;
  double OMtot = OMm + OMl + OMr + OMk;

  if (fabs(OMm - 1.0) < tiny) { // OMm = 1 (Einstein de-Sitter)
    return 1.0 / (1.0 + z);
  } else if ((OMl > (-tiny)) && (fabs(OMl + OMm + OMr - 1.0) < 0.01) && (fabs(wl + 1.0) < tiny)) {
    // this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    // it is taken from liddle et al.
    omegaM_z = OMm * pow(1 + z, 3) / (OMl + OMm * pow(1 + z, 3) + OMr * pow(1 + z, 4));
    dick_z = 2.5 * omegaM_z / (1.0 / 70.0 + omegaM_z * (209 - omegaM_z) / 140.0 + pow(omegaM_z, 4.0 / 7.0));
    dick_0 = 2.5 * OMm / (1.0 / 70.0 + OMm * (209 - OMm) / 140.0 + pow(OMm, 4.0 / 7.0));
    return dick_z / (dick_0 * (1.0 + z));
  } else if ((OMtot < (1 + tiny)) && (fabs(OMl) < tiny)) { // open, zero lambda case (peebles, pg. 53)
    x_0 = 1.0 / (OMm + 0.0) - 1.0;
    dick_0 = 1 + 3.0 / x_0 + 3 * log(sqrt(1 + x_0) - sqrt(x_0)) * sqrt(1 + x_0) / pow(x_0, 1.5);
    x = fabs(1.0 / (OMm + 0.0) - 1.0) / (1 + z);
    dick_z = 1 + 3.0 / x + 3 * log(sqrt(1 + x) - sqrt(x)) * sqrt(1 + x) / pow(x, 1.5);
    return dick_z / dick_0;
  } else if ((OMl > (-tiny)) && (fabs(OMtot - 1.0) < tiny) && (fabs(wl + 1) > tiny)) {
    mlog("IN WANG\n", MLOG_MESG);
    return -1;
  }
  mlog("No growth function!!! Output will be garbage.", MLOG_MESG);
  return -1;
}

// * redshift derivative of the growth function at z * //
double ddicke_dz(double z)
{
  float dz = 1e-10;

  return (dicke(z + dz) - dicke(z)) / dz;
}

float get_Ts(float z, float delta, float TK, float xe, float Jalpha, float* curr_xalpha)
{
  double Trad, xc, xa_tilde;
  double TS, TSold, TSinv;
  double Tceff;

  Trad = TCMB * (1.0 + z);
  xc = xcoll(z, TK, delta, xe);

  if (Jalpha > 1.0e-20) { // * Must use WF effect * //
    TS = Trad;
    TSold = 0.0;
    while (fabs(TS - TSold) / TS > 1.0e-3) {
      TSold = TS;
      xa_tilde = xalpha_tilde(z, Jalpha, TK, TS, delta, xe);
      Tceff = Tc_eff(TK, TS);
      TSinv = (1.0 / Trad + xa_tilde / Tceff + xc / TK) / (1.0 + xa_tilde + xc);

      TS = 1.0 / TSinv;
    }
    *curr_xalpha = (float)xa_tilde;

    if (isnan(xa_tilde)) {
      printf("z = %e delta = %e TK = %e xe = %e Jalpha = %e xa_tilde = %e xc = %e\n",
             z,
             delta,
             TK,
             xe,
             Jalpha,
             xa_tilde,
             xc);
    }
  } else { // * Collisions only * //
    TSinv = (1.0 / Trad + xc / TK) / (1.0 + xc);
    TS = 1.0 / TSinv;
    *curr_xalpha = 0;
  }

  // It can very rarely result in a negative spin temperature. If negative, it is a very small number.
  // Take the absolute value, the optical depth can deal with very large numbers, so ok to be small
  if (TS < 0.) {
    TS = fabs(TS);
  }

  return (float)TS;
}

double xcoll(double z, double TK, double delta, double xe)
{
  return xcoll_HI(z, TK, delta, xe) + xcoll_elec(z, TK, delta, xe) + xcoll_prot(z, TK, delta, xe);
}

double xcoll_HI(double z, double TK, double delta, double xe)
{
  double krate, nH, Trad;
  double xcoll;

  Trad = TCMB * (1.0 + z);
  nH = (1.0 - xe) * No * pow(1.0 + z, 3.0) * (1.0 + delta);
  krate = kappa_10(TK, 0);
  xcoll = T21 / Trad * nH * krate / A10_HYPERFINE;
  return xcoll;
}

// * Note that this assumes Helium ionized same as Hydrogen * //
double xcoll_elec(double z, double TK, double delta, double xe)
{
  double krate, ne, Trad;
  double xcoll;

  Trad = TCMB * (1.0 + z);
  ne = xe * N_b0 * pow(1.0 + z, 3.0) * (1.0 + delta);
  krate = kappa_10_elec(TK, 0);
  xcoll = T21 / Trad * ne * krate / A10_HYPERFINE;
  return xcoll;
}

double xcoll_prot(double z, double TK, double delta, double xe)
{
  double krate, np, Trad;
  double xcoll;

  Trad = TCMB * (1.0 + z);
  np = xe * No * pow(1.0 + z, 3.0) * (1.0 + delta);
  krate = kappa_10_pH(TK, 0);
  xcoll = T21 / Trad * np * krate / A10_HYPERFINE;
  return xcoll;
}

double kappa_10(double TK, int flag)
{
  int i;
  static double tkin[KAPPA_10_NPTS], kap[KAPPA_10_NPTS];
  static gsl_interp_accel* acc;
  static gsl_spline* spline;
  double ans;

  if (flag == 1) { // * Set up spline table * //
    // * Initialize kappa from Zygelman (2005), Table 2, column 4 * //
    tkin[0] = 1.0;
    kap[0] = 1.38e-13;
    tkin[1] = 2.0;
    kap[1] = 1.43e-13;
    tkin[2] = 4.0;
    kap[2] = 2.71e-13;
    tkin[3] = 6.0;
    kap[3] = 6.60e-13;
    tkin[4] = 8.0;
    kap[4] = 1.47e-12;
    tkin[5] = 10.0;
    kap[5] = 2.88e-12;
    tkin[6] = 15.0;
    kap[6] = 9.10e-12;
    tkin[7] = 20.0;
    kap[7] = 1.78e-11;
    tkin[8] = 25.0;
    kap[8] = 2.73e-11;
    tkin[9] = 30.0;
    kap[9] = 3.67e-11;
    tkin[10] = 40.0;
    kap[10] = 5.38e-11;
    tkin[11] = 50.0;
    kap[11] = 6.86e-11;
    tkin[12] = 60.0;
    kap[12] = 8.14e-11;
    tkin[13] = 70.0;
    kap[13] = 9.25e-11;
    tkin[14] = 80.0;
    kap[14] = 1.02e-10;
    tkin[15] = 90.0;
    kap[15] = 1.11e-10;
    tkin[16] = 100.0;
    kap[16] = 1.19e-10;
    tkin[17] = 200.0;
    kap[17] = 1.75e-10;
    tkin[18] = 300.0;
    kap[18] = 2.09e-10;
    tkin[19] = 501.0;
    kap[19] = 2.565e-10;
    tkin[20] = 701.0;
    kap[20] = 2.91e-10;
    tkin[21] = 1000.0;
    kap[21] = 3.31e-10;
    tkin[22] = 2000.0;
    kap[22] = 4.27e-10;
    tkin[23] = 3000.0;
    kap[23] = 4.97e-10;
    tkin[24] = 5000.0;
    kap[24] = 6.03e-10;
    tkin[25] = 7000.0;
    kap[25] = 6.87e-10;
    tkin[26] = 10000.0;
    kap[26] = 7.87e-10;

    // * Convert to logs for interpolation * //
    for (i = 0; i < KAPPA_10_NPTS; i++) {
      tkin[i] = log(tkin[i]);
      kap[i] = log(kap[i]);
    }

    // * Set up spline table * //
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, KAPPA_10_NPTS);
    gsl_spline_init(spline, tkin, kap, KAPPA_10_NPTS);
    return 0;
  }

  if (flag == 2) { // * Clear memory * //
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  if (log(TK) < tkin[0]) { // * Below 1 K, just use that value * //
    ans = kap[0];
  } else if (log(TK) > tkin[KAPPA_10_NPTS - 1]) {
    // * Power law extrapolation * //
    ans = log(exp(kap[KAPPA_10_NPTS - 1]) * pow(TK / exp(tkin[KAPPA_10_NPTS - 1]), 0.381));
  } else { // * Do spline * //
    TK = log(TK);
    ans = gsl_spline_eval(spline, TK, acc);
  }
  return exp(ans);
}

// * Interpolate exact results for kappa_10^eH.  The table goes up to
// * 10^5 K, but values are only accurate for T<2x10^4 K.  From
// * Furlanetto & Furlanetto 2006 * //
double kappa_10_elec(double T, int flag)
{
  static double TK[KAPPA_10_elec_NPTS], kappa[KAPPA_10_elec_NPTS];
  static gsl_interp_accel* acc;
  static gsl_spline* spline;
  double ans;
  int i;
  float curr_TK, curr_kappa;
  FILE* F;

  char fname[STRLEN];

  if (flag == 1) {

    if (run_globals.mpi_rank == 0) {

      sprintf(fname, "%s/kappa_eH_table.dat", run_globals.params.TablesForXHeatingDir);

      // Read in the data
      if (!(F = fopen(fname, "r"))) {
        mlog("Unable to open the kappa_10^eH file at %s\nAborting\n", MLOG_MESG, fname);
        return 0;
      }

      for (i = 0; i < KAPPA_10_elec_NPTS; i++) {
        fscanf(F, "%f %e", &curr_TK, &curr_kappa);
        TK[i] = curr_TK;
        kappa[i] = curr_kappa;
      }

      for (i = 0; i < KAPPA_10_elec_NPTS; i++) {
        TK[i] = log(TK[i]);
        kappa[i] = log(kappa[i]);
      }
    }

    // broadcast the values to all cores
    MPI_Bcast(TK, sizeof(TK), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(kappa, sizeof(kappa), MPI_BYTE, 0, run_globals.mpi_comm);

    // * Set up spline table * //
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, KAPPA_10_elec_NPTS);
    gsl_spline_init(spline, TK, kappa, KAPPA_10_elec_NPTS);
    return 0;
  }

  if (flag == 2) {
    // * Free memory * //
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  T = log(T);
  if (T < TK[0]) { // * Use TK=1 K value if called at lower temperature * //
    ans = kappa[0];
  } else if (T > TK[KAPPA_10_elec_NPTS - 1]) {
    // * Power law extrapolation * //
    ans = kappa[KAPPA_10_elec_NPTS - 1] +
          ((kappa[KAPPA_10_elec_NPTS - 1] - kappa[KAPPA_10_elec_NPTS - 2]) /
           (TK[KAPPA_10_elec_NPTS - 1] - TK[KAPPA_10_elec_NPTS - 2]) * (T - TK[KAPPA_10_elec_NPTS - 1]));
  } else { // * Do spline * //
    ans = gsl_spline_eval(spline, T, acc);
  }
  return exp(ans);
}

// * Interpolate exact results for kappa_10^pH.  The table goes up to
// * 10^5 K, but values are only accurate for T<2x10^4 K.  From
// * Furlanetto & Furlanetto 2006 * //
double kappa_10_pH(double T, int flag)
{
  static double TK[KAPPA_10_pH_NPTS], kappa[KAPPA_10_pH_NPTS];
  static gsl_interp_accel* acc;
  static gsl_spline* spline;
  double ans;
  int i;
  float curr_TK, curr_kappa;

  FILE* F;

  char fname[STRLEN];

  if (flag == 1) {

    if (run_globals.mpi_rank == 0) {

      sprintf(fname, "%s/kappa_pH_table.dat", run_globals.params.TablesForXHeatingDir);

      if (!(F = fopen(fname, "r"))) {
        mlog("Unable to open the kappa_10^pH file at %s\nAborting\n", MLOG_MESG, fname);
        return 0;
      }

      for (i = 0; i < KAPPA_10_pH_NPTS; i++) {
        fscanf(F, "%f %e", &curr_TK, &curr_kappa);
        TK[i] = curr_TK;
        kappa[i] = curr_kappa;
      }
      fclose(F);

      for (i = 0; i < KAPPA_10_pH_NPTS; i++) {
        TK[i] = log(TK[i]);
        kappa[i] = log(kappa[i]);
      }
    }

    // broadcast the values to all cores
    MPI_Bcast(TK, sizeof(TK), MPI_BYTE, 0, run_globals.mpi_comm);
    MPI_Bcast(kappa, sizeof(kappa), MPI_BYTE, 0, run_globals.mpi_comm);

    // * Set up spline table * //
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, KAPPA_10_pH_NPTS);
    gsl_spline_init(spline, TK, kappa, KAPPA_10_pH_NPTS);
    return 0;
  }

  if (flag == 2) {
    // * Free memory * //
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return 0;
  }

  T = log(T);
  if (T < TK[0]) { // * Use TK=1 K value if called at lower temperature * //
    ans = kappa[0];
  } else if (T > TK[KAPPA_10_pH_NPTS - 1]) {
    // * Power law extrapolation * //
    ans = kappa[KAPPA_10_pH_NPTS - 1] +
          ((kappa[KAPPA_10_pH_NPTS - 1] - kappa[KAPPA_10_pH_NPTS - 2]) /
           (TK[KAPPA_10_pH_NPTS - 1] - TK[KAPPA_10_pH_NPTS - 2]) * (T - TK[KAPPA_10_pH_NPTS - 1]));
  } else { // * Do spline * //
    ans = gsl_spline_eval(spline, T, acc);
  }
  ans = exp(ans);
  return ans;
}

// ********************************************************************
// ********************* Wouthuysen-Field Coupling ********************
// ******************************************************************** //

// * NOTE Jalpha is by number * //
double xalpha_tilde(double z, double Jalpha, double TK, double TS, double delta, double xe)
{
  double tgp, Stilde, x;

  tgp = taugp(z, delta, xe);
  Stilde = Salpha_tilde(TK, TS, tgp);
  x = 1.66e11 / (1.0 + z) * Stilde * Jalpha;
  return x;
}

// Compute the Gunn-Peterson optical depth.
double taugp(double z, double delta, double xe)
{
  return 1.342881e-7 / hubble((float)z) * No * pow(1 + z, 3) * (1.0 + delta) * (1.0 - xe);
}

double Salpha_tilde(double TK, double TS, double tauGP)
{
  double xi;
  double ans;

  xi = pow(1.0e-7 * tauGP / TK / TK, 1.0 / 3.0);
  ans = 1.0 - 0.0631789 / TK + 0.115995 / TK / TK - 0.401403 / TS / TK;
  ans += 0.336463 / TS / TK / TK;
  ans /= 1.0 + 2.98394 * xi + 1.53583 * xi * xi + 3.85289 * xi * xi * xi;
  return ans;
}

double Tc_eff(double TK, double TS)
{
  double ans;

  ans = 1.0 / TK + 0.405535 / TK * (1.0 / TS - 1.0 / TK);
  ans = 1.0 / ans;
  return ans;
}

double interpolate_fcoll(double redshift, int snap_i)
{
  double interp_fcoll;
  if (snap_i == 0) {
    // This should never occur...
    interp_fcoll = stored_fcoll[snap_i];
  } else {
    interp_fcoll = stored_fcoll[snap_i - 1] + (redshift - run_globals.ZZ[snap_i - 1]) *
                                                (stored_fcoll[snap_i] - stored_fcoll[snap_i - 1]) /
                                                (run_globals.ZZ[snap_i] - run_globals.ZZ[snap_i - 1]);
  }

  if (interp_fcoll < 0.0) {
    interp_fcoll = 0.0;
  }

  return interp_fcoll;
}
