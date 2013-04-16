#include "meraxes.h"
#include <math.h>

int init_galaxies(galaxy_struct *Gal, int n_halos_max)
{

  int n_total_galaxies = n_halos_max*ALLOCFACTOR;
  
  // malloc the galaxies
  Gal = SID_malloc(sizeof(galaxy_struct)*n_total_galaxies);

  // Initialise the properties
  for(int i_gal=0; i_gal<n_total_galaxies; i_gal++)
  {
    Gal[i_gal].Type          = -1;
    Gal[i_gal].CentralGal    = NULL;
    Gal[i_gal].NextGal       = NULL;
    Gal[i_gal].HaloDesc      = -1;
    Gal[i_gal].HaloNGal      = 0;
    Gal[i_gal].CentralGal    = -1;
    Gal[i_gal].Len           = 0;
    Gal[i_gal].CentralMvir   = 0.0;
    Gal[i_gal].Mvir          = 0.0;
    Gal[i_gal].dM            = 0.0;
    Gal[i_gal].dMdt          = 0.0;
    Gal[i_gal].Rvir          = 0.0;
    Gal[i_gal].Vvir          = 0.0;
    Gal[i_gal].Vmax          = 0.0;
    Gal[i_gal].StellarMass   = 0.0;
    Gal[i_gal].BlackHoleMass = 0.0;
    Gal[i_gal].Cos_Inc       = 0.0;
    Gal[i_gal].MergTime      = 99999.9;
    
    for(int ii=0; ii<3; ii++)
    {
      Gal[i_gal].Pos[ii] = -99999.9;
      Gal[i_gal].Vel[ii] = -99999.9;
    }
    for(int ii=0; ii<NOUT; ii++)
      Gal[i_gal].Sfr[ii] = -99999.9;
  }

  return n_total_galaxies;
  
}

static double calculate_Vvir(run_globals_struct *run_globals, halo_struct *halo)
{
  return sqrt(run_globals->G * halo->Mvir/1.0e10 / halo->Rvir);
}

void copy_halo_to_galaxy(run_globals_struct *run_globals, halo_struct *halo, galaxy_struct *gal)
{
  gal->Type            = halo->Type;
  gal->Len             = halo->Len;
  gal->HaloDesc        = halo->DescIndex;
  gal->Mvir            = halo->Mvir/1.0e10;
  gal->Rvir            = halo->Rvir;
  gal->Vvir            = calculate_Vvir(run_globals, halo);
  gal->Vmax            = halo->Vmax;
}
