#include "meraxes.h"
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif
void timer_start(timer_info* timer)
{
    gettimeofday(&(timer->start), NULL);
}

void timer_stop(timer_info* timer)
{
    gettimeofday(&(timer->stop), NULL);
}

float timer_delta(timer_info timer)
{
    struct timeval diff;
    timersub(&(timer.stop), &(timer.start), &diff);
    return (float)((float)diff.tv_sec + (1e-6 * (float)diff.tv_usec));
}

#ifdef __cplusplus
}
#endif
