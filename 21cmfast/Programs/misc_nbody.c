#include <stdio.h>
#include <stdlib.h>
#include "nbody.h"

void read_zlist(float zlist[MAXSNAPS])
{
  FILE *fin;
  char fname[576];
  char buffer[400];
  char val_buffer[100];
  int i_snap = 0;
  float a_val;
  
  // Parse the expansion factor list and create a redshift list
  sprintf(fname, ROOT_PATH "/" SIM_NAME "/run/" SIM_NAME ".a_list");
  if (!(fin = fopen(fname, "r")))
  {
    fprintf(stderr, "Can't open expansion factor list...\n");
    fflush(stderr);
    fprintf(stderr, "fname = %s\n", fname);
    fflush(stderr);
    ABORT(EXIT_FAILURE);
  }

  i_snap = 0;
  while(!feof(fin))
  {
      *buffer = 0;
      fgets(buffer, 200, fin);
      if(sscanf(buffer, "%s", val_buffer) < 1)
        continue;
      a_val = strtof(val_buffer, NULL);
      zlist[i_snap++] = (1./a_val)-1.;
  }
  fclose(fin);

  if (i_snap!=MAXSNAPS)
  {
    fprintf(stderr, "Didn't read expected number (%d) of expansion factor values...\n", MAXSNAPS);
    fprintf(stderr, "Read %d values...\n", i_snap);
    fprintf(stderr, "Last read value was a=%.3f (z=%.3f)\n", a_val, zlist[i_snap-1]);
    fprintf(stderr, "fname = %s\n", fname);
    ABORT(EXIT_FAILURE);
  } 
}

