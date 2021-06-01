#ifndef COMPUTE_TS_H
#define COMPUTE_TS_H

#include "utils.h"

#define R_XLy_MAX (float)(500)

#ifdef __cplusplus
extern "C"
{
#endif

  void ComputeTs(int snapshot, timer_info* timer_total);

#ifdef __cplusplus
}
#endif

#endif
