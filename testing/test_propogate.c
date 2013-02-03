#define _MAIN
#include "trees.h"

#define BARYON_FRAC 0.17
#define ALLOC_FAC 1.5

static Galaxy* init_galaxies(int n_gals)
{

  Galaxy *galaxies = malloc(sizeof(Galaxy)*n_gals);

  for(int i_gal=0; i_gal<n_gals; i_gal++){
    // Pointers
    galaxies[i_gal].group = NULL;
    galaxies[i_gal].subgroup = NULL;
    galaxies[i_gal].halo = NULL;

    // Properties
    galaxies[i_gal].type = CENTRALT;
    galaxies[i_gal].cold_gas = 0.0;
  }

  return galaxies;
}

static void free_galaxies(Galaxy *galaxies)
{
  free(galaxies);
}

static void evolve_galaxies(Galaxy *galaxies, Subgroup *subgroups, int n_gals)
{
  Halo *old_halo;
  Group *old_group;
  Subgroup *old_subgroup;

  for(int i_gal=0; i_gal<n_gals; i_gal++){
    // Store old pointers for comparison
    old_halo = galaxies[i_gal].halo;
    old_group = galaxies[i_gal].group;
    old_subgroup = galaxies[i_gal].subgroup;

    // Update pointers
    galaxies[i_gal].group = NULL;
    galaxies[i_gal].subgroup = &(subgroups[i_gal]);
    galaxies[i_gal].halo = subgroups[i_gal].halo;
   
    // Do physics
    galaxies[i_gal].cold_gas += (galaxies[i_gal].halo->M_vir)-(old_halo->M_vir) * BARYON_FRAC; 
  }
}

int main(int argc, char* argv[])
{
 
  // Initialise gbpSID
  SID_init(&argc, &argv, NULL);
  SID.fp_log = stdout;
  
  // Quick argument parsing test
  if(argc!=7){
    printf("Usage:\n\t%s sim total_sim_snaps n_every_snaps n_scan_snaps snapshot_first snapshot_last\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  int n_scan_snaps = atoi(argv[4]);
  int snapshot_first = atoi(argv[5]);
  int snapshot_last = atoi(argv[6]);

  // We must hold n_scan_snaps worth of snapshots in order to be able to
  // gaurantee successful connection of halos between snapshots
  Group **groups;
  groups = malloc(sizeof(Group*) * n_scan_snaps);
  Subgroup **subgroups;
  subgroups = malloc(sizeof(Subgroup*) * n_scan_snaps);
  Halo **halos;
  halos = malloc(sizeof(Halo*) * n_scan_snaps);
  TreesHeader *headers;
  headers = malloc(sizeof(TreesHeader) * n_scan_snaps);

  // Read the trees for the first snapshot and allocate/initialise the galaxies
  headers[0] = read_trees(argv[1], atoi(argv[2]), atoi(argv[3]), n_scan_snaps, snapshot_first, &(halos[0]));
  int n_gals = headers[0].n_subgroups;
  Galaxy *galaxies = init_galaxies((int)(n_gals*ALLOC_FAC));
  
  // Evolve the galaxies to the first snapshot
  evolve_galaxies(galaxies, subgroups[0], n_gals);
  free_trees(&(halos[0]));

  // for(int i_snap=1, snapshot=snapshot_first+1; snapshot<=snapshot_last; snapshot++, i_snap++){
  //   headers[i_snap] = read_trees(argv[1], atoi(argv[2]), atoi(argv[3]), n_scan_snaps, snapshot, &(groups[i_snap]), &(subgroups[i_snap]), &(halos[i_snap]));
  //   evolve_galaxies(galaxies, subgroups[0], n_gals);
  //   fprintf(stderr, "galaxies[0].cold_gas = %.3f\n", galaxies[0].cold_gas);
  //   free_trees(&(groups[i_snap]), &(subgroups[i_snap]), &(halos[i_snap]));
  // }
  

  // Free arrays
  free_galaxies(galaxies);
  free(headers);
  free(halos);
  free(subgroups);
  free(groups);

  SID_exit(ERROR_NONE);

  return 0;
}
