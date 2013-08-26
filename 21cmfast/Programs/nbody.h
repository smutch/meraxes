#ifndef _NBODY
#define _NBODY

#define ABORT(sigterm)                                                                   \
do {                                                                                     \
  fprintf(stderr, "in file: %s\tfunc: %s\tline: %i", __FILE__, __FUNCTION__, __LINE__);  \
  exit(sigterm);                                                                         \
} while(0)

#define MAXSNAPS 100
#define FIRSTSNAP 0
#define LASTSNAP 99

#define ROOT_PATH "/home/smutch/Tiamat/"
#define SIM_NAME "Tiamat_full"

#endif
