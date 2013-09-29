#include "meraxes.h"

bool check_reionization_cooling(float cell_ionization, float Vvir)
{

  float Tvir;
  bool   flag;

  if(cell_ionization>0.995)
  {
    Tvir = 35.9 * (Vvir*Vvir); 
    flag = (Tvir < 1e5) ? false : true;
  } else
    flag = true;

  return flag;

}
