#include "meraxes.h"
#include <math.h>

double calculate_merging_time(run_globals_struct *run_globals, galaxy_struct *gal, int snapshot)
{
  galaxy_struct *central;
  double coulomb, mergtime, sat_mass, sat_rad, central_rvir;

  central = gal->FOFGroup->Centralgal;

  if((central == gal) || (central->FOFGroup->Centralgal==-1))
  {
    SID_log_error("Invalid merger...!");
    ABORT(EXIT_FAILURE);
  }

  coulomb = log((double)(central->Len) / (double)(gal->Len) + 1);

	sat_mass = gal->Mvir;

  sat_rad = sqrt(
    pow(central->Pos[0] - gal->Pos[0], 2.0) +
    pow(central->Pos[1] - gal->Pos[1], 2.0) +
    pow(central->Pos[2] - gal->Pos[2], 2.0) );

  // convert to physical length 
  sat_rad /= (1 + run_globals->ZZ[snapshot]);

  central_rvir = central->Rvir;

  if(sat_rad > central_rvir)
    sat_rad = central_rvir;

  if(sat_mass > 0.0)
    mergtime =
    run_globals->params.MergerTimeFactor *
    1.17 * sat_rad * sat_rad * central->Vvir / (coulomb * run_globals->G * sat_mass);
  else
    mergtime = -99.9;

  return mergtime;

}
