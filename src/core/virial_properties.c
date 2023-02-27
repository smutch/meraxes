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

/*double calculate_Rvir_2(double Mvir, double redshift) //from Mvir in 10^10 Msol/h to Rvir in comoving Mpc/h
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
  //fac = 1 / (4 * M_PI / 3.0 * OmegaM * rhocrit);

  return cbrt(Mvir * fac) * zplus1; 
  //return cbrt(Mvir * fac);
}*/

double calculate_Rvir_2(double Mvir, double redshift) //from Mvir in 10^10 Msol/h to Rvir in comoving Mpc/h
{
  double hubble_of_z_sq;
  double rhocrit;
  double fac;
  double Delta;
  double Hubble = run_globals.Hubble;
  double little_h = run_globals.params.Hubble_h;
  double OmegaM = run_globals.params.OmegaM;
  double OmegaK = run_globals.params.OmegaK;
  double OmegaLambda = run_globals.params.OmegaLambda;
  double zplus1 = redshift + 1;
  double GG = 4.3009*1e-9;

  //rhocrit = 3 * (little_h * 100 * little_h * 100 * OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda) / (8 * M_PI * GG);
  rhocrit = 3 * little_h * little_h * 100 * 100 / (8 * M_PI * GG);

  //Delta = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);
  //fac = 1 / (OmegaM*zplus1*zplus1*zplus1 * 4 * M_PI / 3.0 * rhocrit);
  fac = 1 / (4 * M_PI / 3.0 * OmegaM * rhocrit);

  return cbrt(Mvir * 1e10 * fac); 
  //return cbrt(Mvir * fac);
}

/*double calculate_Mvir_2(double Rvir, double redshift) //from Rvir in comoving Mpc/h to 10^10Msol/h
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

  fac = 4.0 / 3.0 * M_PI * rhocrit * Delta;
  //fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM * Delta;

  return pow(Rvir / zplus1 , 3) * fac;
  //return pow(Rvir , 3) * fac;
}*/

double calculate_Mvir_2(double Rvir, double redshift) //from Rvir in comoving Mpc/h to 10^10Msol/h
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
  double GG = 4.3009*1e-9;
  double little_h = run_globals.params.Hubble_h;
  
  //Delta = Delta_vir(redshift);
  //Delta = Omega_z(redshift, OmegaM, OmegaK, OmegaLambda);

  //rhocrit = 3 * (little_h * 100 * little_h * 100 * OmegaM * zplus1 * zplus1 * zplus1 + OmegaK * zplus1 * zplus1 + OmegaLambda) / (8 * M_PI * GG);
  rhocrit = 3 * little_h * little_h * 100 * 100 / (8 * M_PI * GG);

  //fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM * zplus1 * zplus1 * zplus1;
  fac = 4.0 / 3.0 * M_PI * rhocrit * OmegaM;

  return (pow(Rvir , 3) * fac) / 1e10;
  //return pow(Rvir , 3) * fac;
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
  
  return 2.5 * 1e4 * OmegaM * pow(little_h, 2) * pow(Theta, -4); //EH99
}

double Transfer_function(double k) //EH99 You want k in h/Mpc
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
  
  Pk = 2 * M_PI * M_PI * deltah * deltah * pow((SPEED_OF_LIGHT * 1e-5 * scale / (little_h * 100)), 3 + spectral_index) * TF * TF * Dz * Dz / (D0 * D0 * pow(scale, 3));
  
  return Pk;
}

typedef struct
{
  double redshift, HaloMass;
} int_S2_params;

double integrand_S2(double k, void* params)
{
  int_S2_params* p = (int_S2_params*)params;
  double little_h = run_globals.params.Hubble_h;
  
  //double Radius = calculate_Rvir_2(p->HaloMass, p->redshift); //Compute Rvir in cMpc/h
  double Radius = calculate_Rvir_2(p->HaloMass, 0);
  
  double PS = PowerSpectrum(p->redshift, k);
  
  double j1 = (sin(k * Radius) - (k * Radius * cos(k * Radius))) / pow((k * Radius),2);
  
  return k * k * PS / (2 * M_PI * M_PI) * pow(3 * j1 / (k * Radius) , 2);
}

double Sigma(double redshift, double Halo_Mass) // Still a tiny difference
{
  double Hubble = run_globals.Hubble;
  double Normalization = SigmaNorm(0);
  double Sigma8 = run_globals.params.Sigma8; //Need this to check normalization
  double little_h = run_globals.params.Hubble_h;
  //double RedFactor = Growth_Factor(redshift) / Growth_Factor(0);
  
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
    &F, 0, 5000, 1.0 / Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr); //500 should be infinite

  gsl_integration_workspace_free(workspace);
  
  //return Sigma8 * sqrt(result / Normalization) * RedFactor;   
  return Sigma8 * sqrt(result / Normalization); 
}

double SigmaNorm(double redshift) //Need this for normalization 
{

  double Hubble = run_globals.Hubble;
  double little_h = run_globals.params.Hubble_h;
  
  double M8 = calculate_Mvir_2(8.0, 0) ; //Mvir correspondent to a halo of (8Mpc/h virial radius)
  
  int_S2_params p;

  p.redshift = redshift;
  p.HaloMass = M8;
  
  gsl_function F;
  gsl_integration_workspace* workspace;
  
  double norma; 
  double normaerr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_S2;
  F.params = &p;

  gsl_integration_qag(
    &F, 0, 5000, 1.0 / Hubble, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &norma, &normaerr); //500 should be infinite

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
  //mlog("Tot values %d:", MLOG_MESG, x_int_NCFVALS);
  for (i = 0; i < x_int_NCFVALS; i++) {
    if (fabs(x_int_zvals[i] - redshift) <= 0.07) {
      z_index = i;
      for (ii = z_index; ii < x_int_NCFVALS; ii++) {
        if (fabs(x_int_zvals[ii] - redshift) > 0.07 && Radius < MAX_RAD) {
          mlog("Error, you didn't find the radius value!\n", MLOG_MESG);
          exit(1);
          }
        if (fabs((Radius - x_int_radvals[ii]) / Radius) < 0.09) {
          R_index = ii;
          break;
          }
        }
      if (Radius >= MAX_RAD) // If that's the case take the largest value
        R_index = ii;
      break;
      }
    }
    //mlog("Index value %d %d:", MLOG_MESG, z_index, R_index);
    //mlog("Red value %f", MLOG_MESG, x_int_zvals[z_index]);
    //mlog("Radius value %f", MLOG_MESG, x_int_radvals[R_index]);         
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
 
  //mlog("LinearBias %f", MLOG_MESG, LinearBias);
  
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
