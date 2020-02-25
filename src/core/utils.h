#ifndef UTILS_H
#define UTILS_H

#include <sys/time.h>

//! Utility timer struct for GPU runs
typedef struct timer_info {
    struct timeval start;
    struct timeval stop;
} timer_info;

#ifdef _MAIN
float timer_gpu = 0.f;
#else
extern float timer_gpu;
#endif

#ifdef __cplusplus
extern "C" {
#endif

void timer_start(timer_info* timer);
void timer_stop(timer_info* timer);
float timer_delta(timer_info timer);

#ifdef __cplusplus
}
#endif

#endif
