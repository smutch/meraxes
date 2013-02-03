#define _MAIN
#include "trees.h" 

int main(int argc, char *argv[])
{
 
  // Initialise gbpSID
  SID_init(&argc, &argv, NULL);
  SID.fp_log = stdout;
  
  // Quick argument parsing test
  if(argc!=7){
    printf("Usage:\n\t%s sim total_sim_snaps n_every_snaps n_scan_snaps snapshot_first snapshot_last\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  // int snapshot_first = atoi(argv[3]);
  int snapshot_last = atoi(argv[6]);

  Halo *halos;
  TreesHeader header = read_trees(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), snapshot_last, &halos);

#ifdef DEBUG
  // Dump the results as an ascii table
  FILE *fout;
  fout = fopen("test_dump_trees.txt", "w");

  fprintf(fout, "n_groups : %d\n"        , header.n_groups);
  fprintf(fout, "n_subgroups : %d\n"     , header.n_subgroups);
  fprintf(fout, "n_halos_max : %d\n"     , header.n_halos_max);
  fprintf(fout, "n_trees_subgroup : %d\n", header.n_trees_subgroup);
  fprintf(fout, "n_trees_group : %d\n"   , header.n_trees_group);

  fprintf(fout, "# type    id    desc_id    tree_id    file_offset    n_subgroups    M_vir    n_particles    V_max\n");

  Halo *this_halo;
  int type;
  for (int i=0; i<100; i++){
    this_halo = &(halos[i]);
    if(this_halo->n_subgroups>-1)
      type = 0;
    else
      type = 1;
    fprintf(fout, "%d    %d    %d    %d    %d    %d    %.3e    %d    %.3e\n", type, this_halo->id, this_halo->desc_id,
        this_halo->tree_id, this_halo->file_offset, this_halo->n_subgroups, this_halo->M_vir, this_halo->n_particles, this_halo->V_max);
  }

  fclose(fout);

#endif

  // Free allocated arrays
  free_trees(&halos);

  SID_exit(ERROR_NONE);

}
