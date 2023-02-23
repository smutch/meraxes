#include <math.h>

#include "meraxes.h"
#include "virial_properties.h"
#include <gsl/gsl_integration.h>

static float x_int_zvals[x_int_NCFVALS];
static float x_int_radvals[x_int_NCFVALS];
static float x_int_CFvals[x_int_NCFVALS];

static inline double E_z(double z, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double one_plus_z;
  double one_plus_z_sq;
  double one_plus_z_cu;
  double result;

  one_plus_z = 1. + z;
  one_plus_z_sq = one_plus_z * one_plus_z;
  one_plus_z_cu = one_plus_z_sq * one_plus_z;
  result = sqrt(OmegaM * one_plus_z_cu + OmegaK * one_plus_z_sq + OmegaLambda);

  return result;
}

static inline double Omega_z(double redshift, double OmegaM, double OmegaK, double OmegaLambda)
{
  // Function stolen and adapted from gbpCosmo
  double Ez;
  double one_plus_z_cube;

  Ez = E_z(redshift, OmegaM, OmegaK, OmegaLambda);
  one_plus_z_cube = (1. + redshift) * (1. + redshift) * (1. + redshift);

  return OmegaM * one_plus_z_cube / (Ez * Ez);
}

static inline double Delta_vir(double redshift)
{
  // Function stolen and adapted from gbpCosmo
  double x;
  double Omega;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;

  Omega = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);
  x = Omega - 1.;

  return (18. * M_PI * M_PI + 82 * x - 39 * x * x) / Omega;
}

//! Calculates Mvir in internal units (1.e10 h^{-1}Msol), given Tvir (in K) and a redshift (z)
double Tvir_to_Mvir(double T, double z)
{
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;

  double mu; //!< Mean molecular weight (ionized gas)

  if (T < 9.99999e3) // Neutral IGM
    mu = 1.22;
  else // Ionised IGM
    mu = 0.59;

  double z_term = pow((1. + z) / 10., -1.5);
  double T_term = pow(T / 1.98e4, 1.5);
  double cosmo_term = pow(OmegaM / Omega_z(z, OmegaM, OmegaK, OmegaLambda) * Delta_vir(z) / 18. / (M_PI * M_PI), -0.5);
  double mol_term = pow(mu / 0.6, -1.5);

  return 0.01 * mol_term * cosmo_term * T_term * z_term;
}

double calculate_Mvir(double Mvir, int len)
{
  if ((len < 0) && (Mvir > 0))
    return Mvir;
  else
    return (double)len * run_globals.params.PartMass;
}

double hubble_at_snapshot(int snapshot)
{
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = run_globals.ZZ[snapshot] + 1;

  return Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda);
}

double hubble_time(int snapshot)
{
  return 1.0 / hubble_at_snapshot(snapshot);
}

double calculate_Rvir(double Mvir, int snapshot)
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;

  hubble_of_z_sq = pow(hubble_at_snapshot(snapshot), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);

  Delta = Delta_vir(run_globals.ZZ[snapshot]);

  fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);

  return cbrt(Mvir * fac);
}

double calculate_Rvir_2(double Mvir, double redshift) //Mvir in 10^10 Msol/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;

  hubble_of_z_sq = pow(Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);

  Delta = Delta_vir(redshift);

  fac = 1 / (Delta * 4 * M_PI / 3.0 * rhocrit);

  return cbrt(Mvir * fac);
}

double calculate_Mvir_2(double Rvir, double redshift) //Rvir in Mpc/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;
  
  hubble_of_z_sq = pow(Hubble * sqrt(OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);

  fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM;

  return pow(Rvir , 3) * fac;
}

double calculate_gasMass(int snapshot, double length) //length in comoving units
{
  double hubble_of_z_sq;
  double rhocrit;
  double rhob;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaB = OmegaM * run_globals.params.BaryonFrac;
  
  hubble_of_z_sq = pow(hubble_at_snapshot(snapshot), 2);

  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * run_globals.G);
  rhob = rhocrit * OmegaB;

  return rhob * pow((length / (1.0 + run_globals.ZZ[snapshot])), 3.0);
}

double calculate_Vvir(double Mvir, double Rvir)
{
  return sqrt((run_globals.G) * Mvir / Rvir);
}

double calculate_spin_param(halo_t* halo)
{
  double angmom_mag =
    sqrt(halo->AngMom[0] * halo->AngMom[0] + halo->AngMom[1] * halo->AngMom[1] + halo->AngMom[2] * halo->AngMom[2]);
  return angmom_mag / (1.414213562 * halo->Vvir * halo->Rvir);
}

// Here you are adding functions that you need to compute the correlation function (needed for boost the probability of getting metal enriched, motivated by clustering
// As a first test these routines are copied from a previous work of Manu when he was a dumb Master student (the originals were written in Python).
// In the future it might be work to see if these can be improved. Parameters from Eisenstein & Hu 1998 or 1999 (EH98, EH99)

double calculate_zeq(double OmegaM)
{
  double Theta = 2.728 / 2.7;
  double little_h = run_globals.params.Hubble_h;
  
  return 2.5e4 * OmegaM * pow(little_h, 2) * pow(Theta, -4); //EH99
}

double Transfer_function(double k) //EH99
{
  double OmegaM = run_globals.params.OmegaM;
  double OmegaB = OmegaM * run_globals.params.BaryonFrac;
  double OmegaC = OmegaM - OmegaB;
  double little_h = run_globals.params.Hubble_h;
  double Theta = 2.728 / 2.7;
  
  double fc = OmegaC / OmegaM;
  double fb = OmegaB / OmegaM;
  double fcb = fc + fb;
  double alpha = fc / fcb; // Eq. 15
  
  double s_hor = 44.5 * log(9.83 / (OmegaM * pow(little_h, 2))) / pow(1.0 + 10.0 * pow((OmegaB * pow(little_h, 2)), 0.75), 0.5); // Eq. 4
  double Gamma = OmegaM * pow(little_h, 2) * (pow(alpha, 0.5) + (1 - pow(alpha, 0.5)) / (1 + pow(0.43 * k * s_hor, 4))); //Eq. 16
  
  double q = k * Theta * Theta / Gamma;
  double Beta = 1. / (1 - 0.949 * fb); //Eq. 21
  double L = log(exp(1) + 1.84 * Beta * pow(alpha, 0.5) * q); //Eq. 19
  double C = 14.4 + (325. / (1 + 60.5 * pow(q, 1.11))); // Eq. 20
  
  return L / (L + C * q * q); // Eq. 18 and 24
} 

/*double Transfer_function(double k) //EH98
{
  double OmegaM = run_globals.params.OmegaM;
  double OmegaB = OmegaM * run_globals.params.BaryonFrac;
  double OmegaC = OmegaM - OmegaB;
  double fc = OmegaC / OmegaM;
  double fb = OmegaB / OmegaM;
  double fcb = fc + fb;
  double little_h = run_globals.params.Hubble_h;
  double Theta = 2.728 / 2.7;
  
  double zequiv = calculate_zeq(OmegaM);
  double keq = 7.46e-2 * OmegaM * little_h * little_h / (Theta * Theta);
  
  double b1d = 0.313 * pow(OmegaM * little_h * little_h, -0.419) * (1.0 + 0.607 * pow(OmegaM * little_h * little_h, 0.674));
  double b2d = 0.238 * pow(OmegaM * little_h * little_h, 0.223);
  double zd = 1291.0 * pow(OmegaM * little_h * little_h, 0.251) / (1.0 + 0.659 * pow(OmegaM * little_h *little_h, 0.828)) * (1.0 + b1d * pow(OmegaM * little_h * little_h, b2d));
  
  double Rd = 31.5 * OmegaB * little_h * little_h / pow(Theta, 4) / (zd / 1e3);
  double Req = 31.5 * OmegaB * little_h * little_h / pow(Theta, 4) / (zequiv / 1e3);
  
  double s = 2.0 / 3.0 / keq * sqrt(6.0 / Req) * log((sqrt(1.0 + Rd) + sqrt(Rd + Req)) / (1.0 + sqrt(Req)));
  double ksilk = 1.6 * pow(OmegaB * little_h * little_h, 0.52) * pow(OmegaM * little_h * little_h, 0.73) * (1.0 + pow(10.4 * OmegaM * little_h * little_h, -0.95));
  double q = k / 13.41 / keq;
  
  double a1 = pow(46.9 * OmegaM * little_h * little_h, 0.670) * (1.0 + pow(32.1 * OmegaM * little_h * little_h, -0.532));
  double a2 = pow(12.0 * OmegaM * little_h * little_h, 0.424) * (1.0 + pow(45.0 * OmegaM * little_h * little_h, -0.582));
  double ac = pow(a1, -fb) * pow(a2, -fb * fb * fb);
  
  double b1 = 0.944 / (1.0 + pow(458.0 * OmegaM * little_h * little_h, -0.708));
  double b2 = pow(0.395 * OmegaM * little_h * little_h, -0.0266);
  double bc = 1.0 / (1.0 + b1 * (pow(fc / OmegaM, b2) - 1.0));
  
  double y = (1.0 + zequiv) / (1.0 + zd);
  double Gy = y * (-6.0 * sqrt(1.0 + y) + (2.0 + 3.0 * y) * log((sqrt(1.0 + y) + 1.0) / (sqrt(1.0 + y) - 1.0)));
  
  double ab = 2.07 * keq * s * pow(1.0 + Rd, -3.0 / 4.0) * Gy;
  
  double f = 1.0 / (1.0 + pow(k * s / 5.4, 4));

  double C = 14.2 / ac + 386.0 / (1.0 + 69.9 * pow(q, 1.08));

  double T0t = log(exp(1) + 1.8 * bc * q) / (log(exp(1) + 1.8 * bc * q) + C * q * q);
  
  double C1bc = 14.2 + 386.0 / (1.0 + 69.9 * pow(q, 1.08));
  double T0t1bc = log(exp(1) + 1.8 * bc * q) / (log(exp(1) + 1.8 * bc * q) + C1bc * q * q);
  double Tc = f * T0t1bc + (1.0 - f) * T0t;

  double bb = 0.5 + fb + (3.0 - 2.0 * fb) * sqrt((17.2 * OmegaM * little_h * little_h) * (17.2 * OmegaM * little_h * little_h) + 1.0);

  double bnode = 8.41 * pow(OmegaM * little_h * little_h, 0.435);

  double st = s / pow(1.0 + (bnode / k / s) * (bnode / k / s) * (bnode / k / s), (1.0 / 3.0));

  double C11 = 14.2 + 386.0 / (1.0 + 69.9 * pow(q, 1.08));
  double T0t11 = log(exp(1) + 1.8 * q) / (log(exp(1) + 1.8 * q) + C11 * q * q);
  double Tb = (T0t11 / (1.0 + pow(k * s / 5.2, 2)) + ab / (1.0 + pow(bb / k / s, 3)) * exp(-pow(k / ksilk, 1.4))) * sin(k * st) / (k * st);

  double Tk = fb * Tb + fc / OmegaM * Tc;
  
  return Tk; // Eq. 18 and 24
}  */

double integrand_GF(double redshift) //EH99
{
  #define WORKSIZE 1000
  //double zplus1 = run_globals.ZZ[snapshot] + 1;
  double zplus1 = redshift + 1;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaLambda = run_globals.params.OmegaLambda;
  
  return zplus1 / pow(OmegaM * pow(zplus1, 3) + (1 - OmegaM - OmegaLambda) * pow(zplus1, 2) + OmegaLambda, 1.5);
}

double Growth_Factor(double redshift) //Now it works!
{
  double OmegaM = run_globals.params.OmegaM;
  double OmegaLambda = run_globals.params.OmegaLambda;

  double zplus1 = redshift + 1; 
  double zequiv = calculate_zeq(OmegaM);
  double normalization = GF_norm();
  double Pref = pow(OmegaM * pow(zplus1, 3) + (1 - OmegaM - OmegaLambda) * pow(zplus1, 2) + OmegaLambda, 0.5); 
  
  gsl_function F;
  gsl_integration_workspace* workspace;
  double result; 
  double abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_GF;
  F.params = &(run_globals.params);

  gsl_integration_qag(
    &F, redshift, zequiv, 1.0 / run_globals.Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);
  
  return Pref * result / normalization;  
}
double GF_norm() //For Normalization
{
  double OmegaM = run_globals.params.OmegaM;
  double zequiv = calculate_zeq(OmegaM);
  
  gsl_function F;
  gsl_integration_workspace* workspace;
  double norm; 
  double normerr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_GF;
  F.params = &(run_globals.params);

  gsl_integration_qag(
    &F, 0, zequiv, 1.0 / run_globals.Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &norm, &normerr);

  gsl_integration_workspace_free(workspace);
  
  return norm;
}


double PowerSpectrum(double redshift, double scale)
{
  double spectral_index = run_globals.params.SpectralIndex;
  double N = spectral_index - 1;
  double OmegaM = run_globals.params.OmegaM;
  double little_h = run_globals.params.Hubble_h;
  //double scale_h = scale * little_h; // convert k from Mpc to Mpc / h
  
  double deltah = 1.94 * 1.0e-5 * pow(OmegaM, (-0.785 - 0.05 * log(OmegaM))) * exp(-0.95 * N - 0.169 * pow(N,2));
  double TF = Transfer_function(scale); 
  double Dz = Growth_Factor(redshift); 
  double D0 = Growth_Factor(0); 
  double Pk;
  
  Pk = 2 * M_PI * M_PI + deltah * deltah * pow((SPEED_OF_LIGHT * 1e-5 * scale / (little_h * 100)), 3 + spectral_index) * TF * TF * Dz * Dz / (D0 * D0 * pow(scale, 3));
  
  return Pk;
}

typedef struct
{
  double redshift, HaloMass;
} int_S2_params;

double integrand_S2(double k, void* params)
{
  int_S2_params* p = (int_S2_params*)params;
  
  double Radius = calculate_Rvir_2(p->HaloMass, p->redshift);
  //mlog("Radius is %f", MLOG_MESG, Radius);
  
  double PS = PowerSpectrum(p->redshift, k);
  
  double j1 = (sin(k * Radius) - (k * Radius * cos(k * Radius))) / (k * Radius);
  
  return k * k * PS / (2 * M_PI * M_PI) * pow(3 * j1 / (k * Radius), 2);
}

double Sigma(double redshift, double Halo_Mass) // Still a tiny difference
{
  double Hubble = run_globals.Hubble;
  double Normalization = SigmaNorm(redshift);
  double Sigma8 = run_globals.params.Sigma8; //Need this to check normalization
  double little_h = run_globals.params.Hubble_h;
  
  int_S2_params p;

  p.redshift = redshift;
  p.HaloMass = Halo_Mass;
  
  gsl_function F;
  gsl_integration_workspace* workspace;
  
  double result; 
  double abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_S2;
  F.params = &p;

  gsl_integration_qag(
    &F, 0, 500, 1.0 / Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr); //500 should be infinite

  gsl_integration_workspace_free(workspace);
  
  return Sigma8 * sqrt(result / Normalization);   
}

double SigmaNorm(double redshift) //Need this for normalization 
{

  double Hubble = run_globals.Hubble;
  double little_h = run_globals.params.Hubble_h;
  
  double M8 = calculate_Mvir_2(8.0, 0); //Mvir correspondent to a halo of (8Mpc/h virial radius)
  
  int_S2_params p;

  p.redshift = redshift;
  //p.HaloMass = 2.75173293e14; //Halo mass correspondent to Rvir = 8h^-1, this is written extremely badly, use it now just to check that the function is working.
  p.HaloMass = M8;
  
  gsl_function F;
  gsl_integration_workspace* workspace;
  
  double norma; 
  double normaerr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_S2;
  F.params = &p;

  gsl_integration_qag(
    &F, 0, 500, 1.0 / Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &norma, &normaerr); //500 should be infinite

  gsl_integration_workspace_free(workspace);
  
  return norma;   
}

typedef struct
{
  double Radius, redshift;
} int_2CF_params;

double integrand_2pointCF(double k, void* params)
{

  int_2CF_params* p = (int_2CF_params*)params;
  
  //double Radius = calculate_Rvir_2(p->HaloMass, p->redshift);
  
  double PS = PowerSpectrum(p->redshift, k);
  
  //double j1 = (sin(k * Radius) - (k * Radius * cos(k * Radius))) / (k * Radius);
  
  //return k * k * PS / (2 * M_PI * M_PI) * pow(3 * j1 / (k * Radius), 2);
  return k * k * PS / (2 * M_PI * M_PI) * sin(k * p->Radius) / (k * p->Radius);
}

void initialize_interpCF_arrays()
{
  FILE* input_fileCF;
  char input_file_nameCF[500];
  char input_baseCF[] = "SpatialCF.dat";
  char modeCF[10] = "r";
  
  int i;
  
  if (run_globals.mpi_rank == 0) {
    sprintf(input_file_nameCF,"%s/%s", run_globals.params.TablesForXHeatingDir, input_baseCF); // ATM is in the same location, you might want change it later!
    input_fileCF = fopen(input_file_nameCF, modeCF);
    
    if (input_fileCF == NULL) {
        mlog("Can't open input file %s!\n", MLOG_MESG, input_file_nameCF);
        exit(1);
      }
    
    // Read in data table
      for (i = 0; i < x_int_NCFVALS; i++) {
        fscanf(input_fileCF,
               "%g %g %g",
               &x_int_zvals[i],
               &x_int_radvals[i],
               &x_int_CFvals[i]);
      }

      fclose(input_fileCF);
    }
    
  // broadcast the values to all cores
  MPI_Bcast(&x_int_zvals, sizeof(x_int_zvals), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_radvals, sizeof(x_int_radvals), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&x_int_CFvals, sizeof(x_int_CFvals), MPI_BYTE, 0, run_globals.mpi_comm);
}

double read_SpatialCF(double redshift, double Radius) //Radius in cMpc/h
{
  int z_index = 0;
  int R_index = 0;
  int i = 0;
  int ii = 0;
  
  for (i = 0; i < x_int_NCFVALS; i++) {
    if (abs(x_int_zvals[i] - redshift) <= 0.05) {
      z_index = i;
      for (ii = z_index; ii < x_int_NCFVALS; ii++) {
        if (abs(x_int_zvals[ii] - redshift) > 0.05 && Radius < MAX_RAD) {
          mlog("Error, you didn't find the radius value!\n", MLOG_MESG);
          exit(1);
          }
        if (abs((Radius - x_int_radvals[ii]) / Radius) < 0.1) {
          R_index = ii;
          break;
          }
        }
      if (Radius >= MAX_RAD) // If that's the case take the largest value
        R_index = ii;
      break;
      }
    }
    mlog("Index value %d %d:", MLOG_MESG, z_index, R_index);
    mlog("Red value %f :", MLOG_MESG, x_int_zvals[6078]);
    mlog("Radius value %f :", MLOG_MESG, x_int_radvals[R_index]);         
  return x_int_CFvals[R_index];
}

double TwoPointCF_2(double redshift, double Halo_Mass, double Radius) // 2nd attempt, reading from tables
{

  //int_2CF_params p;

  //p.redshift = redshift; 
  //p.Radius = Radius;
  double Hubble = run_globals.Hubble;
  double nuu = nuc(redshift, Halo_Mass);
  //double nuu = nuc_2(redshift, Halo_Mass);
  double DeltaCrit = 1.686 / Growth_Factor(redshift); // Double check this later, in Mo & White they just do 1.686 * (1 + redshift_2)
  //double DeltaCrit = 1.686 * (1 + redshift);
  
  double SpatialCFval = read_SpatialCF(redshift, Radius);
  double LinearBias = 1 + ((nuu * nuu - 1) / DeltaCrit);
  
  /*gsl_function F;
  gsl_integration_workspace* workspace;
  
  double result; 
  double abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_2pointCF;
  F.params = &p;

  gsl_integration_qag(
    &F, 1.0e-7, 1000, 1.0 / Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS61, workspace, &result, &abserr); //500 should be infinite, qao or qag?

  gsl_integration_workspace_free(workspace);*/
  
  mlog("AutoCF %f", MLOG_MESG, SpatialCFval);
  
  return SpatialCFval * LinearBias * LinearBias;   
}

double nuc(double redshift, double Halo_Mass)
{
  double DeltaCrit = 1.686 / Growth_Factor(redshift);
  double ss = Sigma(redshift, Halo_Mass);
  
  return DeltaCrit / ss;
}

double nuc_2(double redshift, double Halo_Mass)
{
  double DeltaCrit = 1.686 * (1 + redshift);
  double ss = Sigma(redshift, Halo_Mass);
  
  return DeltaCrit / ss;
}

double R0(double redshift, double Halo_Mass) // 7,8 from Barone-Nugent and 12 from Sheth&Tormen, result in Mpc/h, Probably there is a problem with this!
{
  double little_h = run_globals.params.Hubble_h;
  double Sigma8 = run_globals.params.Sigma8;
  
  double DeltaCrit = 1.686 / Growth_Factor(redshift);
  double nuu = nuc(redshift, Halo_Mass);
  
  double gamma = 1.6;
  double a = 0.707;
  double p = 0.3;
  
  return ( 8.0 * pow(Sigma8 * Sigma8 / 72.0 * (3 - gamma) * (4 - gamma) * (6 - gamma) * pow(2.0, gamma) * pow(1 + (a * nuu - 1) / DeltaCrit + 2 * p / nuu / (1 + a * pow(nuu, p)), 2), (1.0 / gamma)));
  
}

double TwoPointCF(double Radius, double Corr_length) //both in Mpc/h
{
  double gamma = 1.6;
  
  return (pow(Radius / Corr_length, -gamma));

}
