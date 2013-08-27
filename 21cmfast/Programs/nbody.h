#ifndef _NBODY
#define _NBODY

#define ABORT(sigterm)                                                                   \
do {                                                                                     \
  fprintf(stderr, "in file: %s\tfunc: %s\tline: %i", __FILE__, __FUNCTION__, __LINE__);  \
  exit(sigterm);                                                                         \
} while(0)

#define MAXSNAPS 59
#define FIRSTSNAP 0
#define LASTSNAP 58

#define ROOT_PATH "/home/gpoole/GiggleZ/"
#define SIM_NAME "GiggleZ_MR"

#endif
