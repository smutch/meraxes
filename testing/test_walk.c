#define _MAIN
#include "trees.h"


int main(int argc, char* argv[])
{
  // Quick argument parsing test
  if(argc!=7){
    printf("Usage:\n\t%s sim total_sim_snaps n_every_snaps n_scan_snaps snapshot_first snapshot_last\n", argv[0]);
    return EXIT_FAILURE;
  }
  char *sim                   = argv[1];
  const int   total_sim_snaps = atoi(argv[2]);
  const int   n_every_snaps   = atoi(argv[3]);
  const int   n_scan_snaps    = atoi(argv[4]);
  const int   snapshot_first  = atoi(argv[5]);
  const int   snapshot_last   = atoi(argv[6]);

  // Initialise gbpSID
  SID_init(&argc, &argv, NULL);
  // SID.fp_log = stdout;
  
  // We must hold n_scan_snaps+1 worth of snapshots in order to be able to
  // gaurantee successful connection of halos between snapshots
  Halo        **halos;
  TreesHeader  *headers;
  halos     = SID_malloc(sizeof(Halo *)      *(n_scan_snaps+1));
  headers   = SID_malloc(sizeof(TreesHeader) *(n_scan_snaps+1));

  // Initialise the halo pointers
  for(int i=0; i<n_scan_snaps+1; i++)
    halos[i] = NULL;

  // Read the first n_scan_snaps+1 snapshots in...
  int last_read_snap = snapshot_first;
  for(int i_snap=0, snapshot=snapshot_first; (snapshot<=snapshot_last) && (i_snap<(n_scan_snaps+1)); snapshot++, i_snap++){
    headers[snapshot%(n_scan_snaps+1)] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, snapshot, &(halos[snapshot%(n_scan_snaps+1)]));
    last_read_snap = snapshot;
  }

  // DEBUG - Quick check
  for(int i=0; i<(n_scan_snaps+1); i++)
  {
    for(int j=0; j<(headers[i].n_groups+headers[i].n_subgroups); j++)
    {
      if(halos[i][j].id!=j+1)
      {
        SID_log("Epic fail! -> halos[%d][%d].id=%d !!!", SID_LOG_COMMENT, i, j, halos[i][j].id);
        ABORT(34494);
      }
    }
  }

  
  // Now loop through each snapshot and try connecting the halos
  // TODO: Fix this for when we have multiple large file_offset values.
  /*
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
     */

  FILE *temp_fout = fopen("test_dump_walk.txt", "w");
  fprintf(temp_fout, "# snapshot  ID  desc_ID  file_offset  type  M_vir\n");
  int halo_id = 0;
  Halo chalo;
  for(int snapshot=snapshot_first; snapshot<=snapshot_last; snapshot++){
    chalo = halos[snapshot%(n_scan_snaps+1)][halo_id];
    fprintf(temp_fout, "%d    %d    %d    %d    %d    %.3e\n", snapshot, chalo.id, chalo.desc_id, chalo.file_offset, chalo.type, chalo.M_vir);
    halo_id = chalo.desc_id;

    if (halo_id<0)
      break;

    if(halos[snapshot%(n_scan_snaps+1)]!=NULL)
      free_trees(&(halos[snapshot%(n_scan_snaps+1)]));
    last_read_snap++;
    if (last_read_snap<=snapshot_last)
      headers[snapshot%(n_scan_snaps+1)] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, last_read_snap, &(halos[snapshot%(n_scan_snaps+1)]));

    // If the halo skips snapshots then interpolate...
    if (chalo.file_offset > 1)
    {
      double M_vir = chalo.M_vir;
      double M_vir_step = (halos[(snapshot+chalo.file_offset)%(n_scan_snaps+1)][halo_id].M_vir - M_vir)/(double)chalo.file_offset; 
      for(int i_skip=1; i_skip<(chalo.file_offset+1); i_skip++)
      {
        snapshot++;
        M_vir+=(M_vir+M_vir_step); 
        fprintf(temp_fout, "%d    %d    %d    %d    %d    %.3e\n", snapshot+1, -1, -1, chalo.file_offset-i_skip, -999, M_vir);
        if(halos[snapshot%(n_scan_snaps+1)]!=NULL)
          free_trees(&(halos[snapshot%(n_scan_snaps+1)]));
        last_read_snap++;
        if (last_read_snap<=snapshot_last)
          headers[snapshot%(n_scan_snaps+1)] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, last_read_snap, &(halos[last_read_snap%(n_scan_snaps+1)]));
      }
    }
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
  {
    if (&(halos[i_snap])!=NULL)
      free_trees(&(halos[i_snap]));
  }
  free(headers);
  free(halos);

  SID_exit(ERROR_NONE);

}
