#ifndef _UTILS_H
#define _UTILS_H

#include <time.h>

struct timer_info { time_t start; time_t stop;};
typedef struct timer_info timer_info;

void timer_start(timer_info *timer);
void timer_stop (timer_info *timer);
long timer_delta(timer_info timer);
#endif
