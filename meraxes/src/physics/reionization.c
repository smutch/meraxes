#include "meraxes.h"

bool check_reionization_cooling(galaxy_struct *gal)
{

  double Tvir;
  double cell_ionization = gal->CellIonization;
  bool   flag;

  if(cell_ionization>0.995)
  {
    Tvir = 35.9 * (gal->Vvir*gal->Vvir); 
    flag = (Tvir < 1e5) ? false : true;
  } else
    flag = true;

  return flag;

}
