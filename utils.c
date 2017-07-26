#include <time.h>
#include "utils.h"

void timer_start(timer_info *timer){
    time(&(timer->start));
}

void timer_stop (timer_info *timer){
    time(&(timer->stop));
}

long timer_delta(timer_info timer){
    return((long)difftime(timer.stop,timer.start));
}

