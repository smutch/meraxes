#include "meraxes.h"
#include <math.h>

double calculate_merging_time(run_globals_struct *run_globals, galaxy_struct *sat, int snapshot)
{
  galaxy_struct *parent;
  double coulomb, mergtime, sat_mass, sat_rad, parent_rvir;

  parent = sat->Halo->Galaxy;
  gal->MergerTarget = parent;

  if(parent == gal)
  {
    SID_log_error("Invalid merger...!");
    ABORT(EXIT_FAILURE);
  }

  coulomb = log((double)(parent->Len) / (double)(sat->Len) + 1);

	sat_mass = sat->Mvir;

  sat_rad = sqrt(
    pow(parent->Pos[0] - sat->Pos[0], 2.0) +
    pow(parent->Pos[1] - sat->Pos[1], 2.0) +
    pow(parent->Pos[2] - sat->Pos[2], 2.0) );

  // convert to physical length 
  sat_rad /= (1 + run_globals->ZZ[snapshot]);

  parent_rvir = parent->Rvir;

  if(sat_rad > parent_rvir)
    sat_rad = parent_rvir;

  if(sat_mass > 0.0)
    mergtime =
    run_globals->params.MergerTimeFactor *
    1.17 * sat_rad * sat_rad * parent->Vvir / (coulomb * run_globals->G * sat_mass);
  else
    mergtime = -99.9;

  return mergtime;

}
