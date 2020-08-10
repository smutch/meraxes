#ifndef FIND_HII_BUBBLES_H
#define FIND_HII_BUBBLES_H

#include "utils.h"

#ifdef __cplusplus
extern "C"
{
#endif

  double RtoM(double R);
  void find_HII_bubbles(int snapshot, timer_info* timer_total);

#ifdef __cplusplus
}
#endif

#endif
