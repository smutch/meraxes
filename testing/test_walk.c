#define _MAIN
#include "trees.h"


int main(int argc, char* argv[])
{
  // Quick argument parsing test
  if(argc!=7){
    printf("Usage:\n\t%s sim total_sim_snaps n_every_snaps n_scan_snaps snapshot_first snapshot_last\n", argv[0]);
    return EXIT_FAILURE;
  }
  char *sim             = argv[1];
  const int   total_sim_snaps = atoi(argv[2]);
  const int   n_every_snaps   = atoi(argv[3]);
  const int   n_scan_snaps    = atoi(argv[4]);
  const int   snapshot_first  = atoi(argv[5]);
  const int   snapshot_last   = atoi(argv[6]);

  // Initialise gbpSID
  SID_init(&argc, &argv, NULL);
  SID.fp_log = stdout;
  
  // We must hold n_scan_snaps worth of snapshots in order to be able to
  // gaurantee successful connection of halos between snapshots
  Halo        **halos;
  TreesHeader  *headers;
  halos     = malloc(sizeof(Halo *)      *(n_scan_snaps+1));
  headers   = malloc(sizeof(TreesHeader) *(n_scan_snaps+1));

  // Initialise the halo pointers
  for(int i=0; i<n_scan_snaps+1; i++)
    halos[i] = NULL;

  // Read the first n_scan_snaps+1 snapshots in...
  for(int i_snap=0, snapshot=snapshot_first; snapshot<=snapshot_last; snapshot++, i_snap=(i_snap+1)%(n_scan_snaps+1)){
    headers[i_snap] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, snapshot, &(halos[i_snap]));
  }
  
  // Now loop through each snapshot and try connecting the halos
  // TODO: Fix this for when we have multiple large file_offset values.
  FILE *temp_fout = fopen("test_dump_walk.txt", "w");
  int test_halo_id = 0;
  Halo chalo;
  for(int snapshot=snapshot_first, i_snap=0, i_roll=0; snapshot<=snapshot_last; snapshot++, i_roll=(i_roll+1)%(n_scan_snaps+1)){
    chalo = halos[i_snap+i_roll][test_halo_id];
    fprintf(temp_fout, "%d    %d    %d    %d    %.3e\n", snapshot, chalo.id, chalo.desc_id, chalo.file_offset, chalo.M_vir);
    test_halo_id = chalo.desc_id;
    if(halos[i_snap+i_roll]!=NULL)
      free_trees(&(halos[i_snap+i_roll]));
    headers[i_snap+i_roll] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, snapshot, &(halos[i_snap+i_roll]));
    i_snap = chalo.file_offset-1;
  }
  fclose(temp_fout);

  // // Take the first central halo as a test and connect the history
  // int idx = 0;
  // for(int i_snap=0, snapshot=snapshot_first; snapshot<=snapshot_last; snapshot++, i_snap++){
  //   printf("snap %d => id=%d, Mvir=%.3e, n_subgroups=%d, desc_id=%d\n", 
  //       i_snap, groups[i_snap][idx].id, groups[i_snap][idx].halo->M_vir, groups[i_snap][idx].n_subgroups, groups[i_snap][idx].desc_id); 
  //   idx = groups[i_snap][idx].desc_id;
  //   if(idx<0)
  //     break;
  // }

  // Free arrays
  for(int i_snap=0; i_snap<n_scan_snaps; i_snap++)
    free_trees(&(halos[i_snap]));
  free(headers);
  free(halos);

  SID_exit(ERROR_NONE);

}
