#define _MAIN
#include "trees.h"

int check_for_merger(int flag)
{
  return (flag&TREE_CASE_MERGER)==TREE_CASE_MERGER;
}

int main(int argc, char* argv[])
{
  // Quick argument parsing test
  if(argc!=8){
    printf("Usage:\n\t%s sim total_sim_snaps n_every_snaps n_scan_snaps snapshot_first snapshot_last init_index\n", argv[0]);
    return EXIT_FAILURE;
  }
  char *sim                   = argv[1];
  const int   total_sim_snaps = atoi(argv[2]);
  const int   n_every_snaps   = atoi(argv[3]);
  const int   n_scan_snaps    = atoi(argv[4]);
  const int   snapshot_first  = atoi(argv[5]);
  const int   snapshot_last   = atoi(argv[6]);
  const int   init_index      = atoi(argv[6]);

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
  for(int i_snap=0, snapshot=snapshot_first; (snapshot<=snapshot_last) && (i_snap<(n_scan_snaps+1)); snapshot++, i_snap++)
  {
    headers[snapshot%(n_scan_snaps+1)] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, snapshot, &(halos[snapshot%(n_scan_snaps+1)]));
    SID_log("DEBUG: snapshot%(n_scan_snaps+1)=%d", SID_LOG_COMMENT, snapshot%(n_scan_snaps+1));
    last_read_snap = snapshot;
  }

  
  // Now loop through each snapshot and try connecting the halos
  FILE *temp_fout = fopen("test_dump_walk.txt", "w");
  fprintf(temp_fout, "# snapshot  ID  index  desc_ID  file_index  file_offset  type  M_vir  next_M_vir  interpolated?\n");
  int halo_index = init_index;
  int i_roll = 0;
  Halo chalo;
  for(int snapshot=snapshot_first; snapshot<=snapshot_last; snapshot++){
    i_roll = snapshot%(n_scan_snaps+1);

    // Make sure we have some halos in this snapshot!
    if ((headers[i_roll].n_subgroups<1) || (headers[i_roll].n_subgroups<=halo_index))
    {
      if(halos[i_roll]!=NULL)
        free_trees(&(halos[i_roll]));
      last_read_snap++;
      if (last_read_snap<=snapshot_last)
        headers[i_roll] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, last_read_snap, &(halos[snapshot%(n_scan_snaps+1)]));
      continue;
    }

    chalo = halos[i_roll][halo_index];
    halo_index = chalo.file_index;
    fprintf(temp_fout, "%d    %d    %d    %d    %d    %d    %d    %.3e    %.3e    0\n", snapshot, chalo.id, halo_index, chalo.desc_id, chalo.file_index, 
        chalo.file_offset, chalo.type, 
        chalo.M_vir, halos[(snapshot+chalo.file_offset)%(n_scan_snaps+1)][halo_index].M_vir);

    if (chalo.id<0)
      break;

    if (check_for_merger(chalo.tree_flags))
      break;

    if(halos[i_roll]!=NULL)
      free_trees(&(halos[i_roll]));
    last_read_snap++;
    if (last_read_snap<=snapshot_last)
      headers[i_roll] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, last_read_snap, &(halos[snapshot%(n_scan_snaps+1)]));

    // If the halo skips snapshots then interpolate...
    // TODO - This needs debugged...
    if (chalo.file_offset > 1)
    {
      double M_vir = chalo.M_vir;
      double M_vir_step = (halos[(snapshot+chalo.file_offset)%(n_scan_snaps+1)][halo_index].M_vir - M_vir)/(double)chalo.file_offset; 

      // DEBUG
      printf("INTERPOLATION:\n");
      printf("snapshot = %d\n", snapshot);
      printf("cur Mvir = %.2e\n", M_vir);
      printf("future Mvir = %.2e\n", halos[(snapshot+chalo.file_offset)%(n_scan_snaps+1)][halo_index].M_vir);
      printf("skipped snaps = %d\n", chalo.file_offset);
      printf("Mvir step = %.2e\n\n", M_vir_step);

      for(int i_skip=1; i_skip<chalo.file_offset; i_skip++)
      {
        snapshot++;
        i_roll = snapshot%(n_scan_snaps+1);
        M_vir+=M_vir_step; 
        fprintf(temp_fout, "%d    %d    %d    %d    %d    %d    %d    %.3e    %d   1\n", snapshot, -1, -1, -1, -1, chalo.file_offset-i_skip, -999, M_vir, -1);
        if(halos[i_roll]!=NULL)
          free_trees(&(halos[i_roll]));
        last_read_snap++;
        if (last_read_snap<=snapshot_last)
          headers[i_roll] = read_trees(sim, total_sim_snaps, n_every_snaps, n_scan_snaps, last_read_snap, &(halos[last_read_snap%(n_scan_snaps+1)]));
      }
    }
  }
  fclose(temp_fout);

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
