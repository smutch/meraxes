/*
 * This code is an amalgamation of the requisite functions for X-ray heating taken from
 * 21cmFAST. Specifically, from heating_helper_progs.c and elec_interp.c.
 * The original code was written by Andrei Mesinger.
 * Written to work within Meraxes by Bradley Greig.
 */

#ifndef XRAY_HEATING_FUNCTIONS_H
#define XRAY_HEATING_FUNCTIONS_H

#include "meraxes.h"

// Below gives grid sizes for the interpolation arrays
#ifndef _x_int_VARIABLES_DEFINED
#define x_int_NXHII 14
#define x_int_NENERGY 258
#define _x_int_VARIABLES_DEFINED
#endif

#define Pop (int)(2)
#define Pop2_ion run_globals.params.physics.ReionNionPhotPerBary
#define Pop3_ion (float)(44021)

#define NSPEC_MAX (int)23
#define RECFAST_NPTS (int)501
#define KAPPA_10_NPTS (int)27
#define KAPPA_10_elec_NPTS (int)20
#define KAPPA_10_pH_NPTS (int)17

#define CLUMPING_FACTOR                                                                                                \
  (double)(2) // sub grid scale.  note that if you want to run-down from a very high redshift (>50), you should set this
              // to one..
#define T21 (double)(0.0628)                /* temperature corresponding to the 21cm photon */
#define A10_HYPERFINE (double)(2.85e-15)    /* spontaneous emission coefficient in s^-1 */
#define Ly_alpha_HZ (double)(2.46606727e15) /* frequency of Lyalpha */
#define NU_over_EV (double)(1.60217646e-12 / PLANCK)
#define NUIONIZATION (double)(13.60 * NU_over_EV)     /* ionization frequency of H */
#define HeI_NUIONIZATION (double)(24.59 * NU_over_EV) /* ionization frequency of HeI */
#define HeII_NUIONIZATION (double)(NUIONIZATION * 4)  /* ionization frequency of HeII */

// Define some global variables; yeah i know it isn't "good practice" but doesn't matter
// NB. Not written by smutch!!! ;)
#ifdef _XRAY_HEATING_FUNCTIONS_C
double x_e_ave;
double dt_dzpp;
double dt_dzp;
double* zpp_edge;

// Have this arbitrarily large for now. Will do this properly later
double stored_fcoll[1000];
double* sum_lyn;
double* sum_lyn_LW;
double growth_factor_zp;
double dgrowth_factor_dzp;
double const_zp_prefactor_GAL;
double const_zp_prefactor_QSO;
float x_int_XHII[x_int_NXHII];
#else
extern double x_e_ave;
extern double dt_dzpp;
extern double dt_dzp;
extern double* zpp_edge;
extern double stored_fcoll[1000];
extern double* sum_lyn;
extern double* sum_lyn_LW;
extern double growth_factor_zp;
extern double dgrowth_factor_dzp;
extern double const_zp_prefactor_GAL;
extern double const_zp_prefactor_QSO;
extern float x_int_XHII[x_int_NXHII];
#endif

#ifdef __cplusplus
extern "C"
{
#endif

  // Initialization; must be called once to
  void initialize_interp_arrays();

  // Primary functions to compute heating fractions and number of Lya photons or ionization produced,
  // Note that En is the energy of the *primary* photon, so the energy in the initial ionization is
  // included in all these.
  // All energies are in eV.
  // xHII_call is the desired ionized fraction.
  float interp_fheat(float En, float xHII_call);
  float interp_n_Lya(float En, float xHII_call);
  float interp_nion_HI(float En, float xHII_call);
  float interp_nion_HeI(float En, float xHII_call);
  float interp_nion_HeII(float En, float xHII_call);

  int locate_energy_index(float En);
  int locate_xHII_index(float xHII_call);

  double dT_comp(double z, double TK, double xe);

  double dicke(double z);
  double alpha_A(double T);

  /* initialization routine */
  int init_heat();

  /* destruction/deallocation routine */
  void destruct_heat();

  /* returns the spectral emissity */
  double spectral_emissivity(double nu_norm, int flag);

  /* Ionization fraction from RECFAST. */
  double xion_RECFAST(float z, int flag);

  /* IGM temperature from RECFAST; includes Compton heating and adiabatic expansion only. */
  double T_RECFAST(float z, int flag);

  double HI_ion_crosssec(double nu);
  double HeII_ion_crosssec(double nu);
  double HeI_ion_crosssec(double nu);

  /* Calculates the optical depth for a photon arriving at z = zp with frequency nu, emitted at z = zpp */
  double tauX(double nu, double x_e, double zp, double zpp, double fcoll, double HI_filling_factor_zp, int snap_i);

  /* The total weighted HI + HeI + HeII  cross-section in pcm^-2 */
  double species_weighted_x_ray_cross_section(double nu, double x_e);

  /* Returns the frequency threshold where \tau = 1 between zp and zpp,
     in the IGM with mean electron fraction x_e */
  double nu_tau_one(double zp, double zpp, double x_e, double fcoll, double HI_filling_factor_zp, int snap_i);

  /* Main integral driver for the frequency integral in the evolution equations */
  double integrate_over_nu(double zp,
                           double local_x_e,
                           double lower_int_limit,
                           double thresh_energy,
                           double spec_index,
                           int FLAG);

  /* Returns the maximum redshift at which a Lyn transition contributes to Lya
     flux at z */
  float zmax(float z, int n);

  /* Returns frequency of Lyman-n, in units of Lyman-alpha */
  double nu_n(int n);

  /* Returns recycling fraction (=fraction of photons converted into Lyalpha for Ly-n resonance */
  double frecycle(int n);

  /* returns the spin temperature */
  float get_Ts(float z, float delta, float TK, float xe, float Jalpha, float* curr_xalpha);

  /* Spin Temperature helper functions */
  double xcoll(double z, double TK, double delta, double xe);
  double xcoll_HI(double z, double TK, double delta, double xe);
  double xcoll_elec(double z, double TK, double delta, double xe);
  double xcoll_prot(double z, double TK, double delta, double xe);
  double kappa_10_pH(double T, int flag);
  double kappa_10_elec(double T, int flag);

  double kappa_10(double TK, int flag);
  double xalpha_tilde(double z, double Jalpha, double TK, double TS, double delta, double xe);
  double taugp(double z, double delta, double xe);
  double Salpha_tilde(double TK, double TS, double tauGP);
  double Tc_eff(double TK, double TS);

  double ddicke_dz(double z);
  double dtdz(float z);
  double drdz(float z);

  double gettime(double z);
  double hubble(float z);

  double interpolate_fcoll(double redshift, int snap_i);

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
                 double deriv[]);

#ifdef __cplusplus
}
#endif

#endif
