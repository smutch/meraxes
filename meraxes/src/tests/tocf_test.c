#define _MAIN
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <meraxes.h>
#include <cmocka.h>

typedef struct state_t {
  galaxy_t *gals;
  int n_gals;
  int global_n_gals;
} state_t;
#define mystate ((state_t *)*state)


static int setup_map_galaxies_to_slabs(void **state) {
  *state = SID_malloc(sizeof(state_t));

  tocf_params.HII_dim = 64;
  tocf_params.uvb_feedback = 1;
  run_globals.params.TOCF_Flag = 1;
  run_globals.params.BoxSize = 10.;  // Not the size of Tiamat but easy for checking

  malloc_reionization_grids();

  mystate->global_n_gals = 8;

  switch(SID.My_rank)
  {
    case 0:
      {
        mystate->n_gals = 4;
        float xpos[] = {43, 21, 67, 88};
        mystate->gals = SID_malloc(mystate->n_gals * sizeof(galaxy_t));
        galaxy_t *gals = mystate->gals;
        for(int ii=0; ii < mystate->n_gals; ii++)
        {
          if (ii < mystate->n_gals-1)
            gals[ii].Next = &(gals[ii+1]);
          else
            gals[ii].Next = NULL;

          gals[ii].Type = 0;
          gals[ii].Pos[0] = xpos[ii];
          gals[ii].Pos[1] = gals[ii].Pos[0];
          gals[ii].Pos[2] = gals[ii].Pos[0];

          for(int jj=0; jj<3; jj++)
            SID_log("R%d: [%d][%d] -> %.2f", SID_LOG_COMMENT|SID_LOG_ALLRANKS, SID.My_rank, ii, jj, gals[ii].Pos[jj]);
        }
      }
      break;

    case 1:
      {
        mystate->n_gals = 3;
        float xpos[] = {1, 2, 99, 4};
        mystate->gals = SID_malloc(mystate->n_gals * sizeof(galaxy_t));
        galaxy_t *gals = mystate->gals;
        for(int ii=0; ii < mystate->n_gals; ii++)
        {
          if (ii < mystate->n_gals-1)
            gals[ii].Next = &(gals[ii+1]);
          else
            gals[ii].Next = NULL;

          gals[ii].Type = 0;
          gals[ii].Pos[0] = xpos[ii];
          gals[ii].Pos[1] = gals[ii].Pos[0];
          gals[ii].Pos[2] = gals[ii].Pos[0];

          for(int jj=0; jj<3; jj++)
            SID_log("R%d: [%d][%d] -> %.2f", SID_LOG_COMMENT|SID_LOG_ALLRANKS, SID.My_rank, ii, jj, gals[ii].Pos[jj]);
        }
      }
      break;

    case 2:
      {
        mystate->n_gals = 1;
        float xpos[] = {50};
        mystate->gals = SID_malloc(mystate->n_gals * sizeof(galaxy_t));
        galaxy_t *gals = mystate->gals;
        for(int ii=0; ii < mystate->n_gals; ii++)
        {
          if (ii < mystate->n_gals-1)
            gals[ii].Next = &(gals[ii+1]);
          else
            gals[ii].Next = NULL;

          gals[ii].Type = 0;
          gals[ii].Pos[0] = xpos[ii];
          gals[ii].Pos[1] = gals[ii].Pos[0];
          gals[ii].Pos[2] = gals[ii].Pos[0];

          for(int jj=0; jj<3; jj++)
            SID_log("R%d: [%d][%d] -> %.2f", SID_LOG_COMMENT|SID_LOG_ALLRANKS, SID.My_rank, ii, jj, gals[ii].Pos[jj]);
        }
      }
      break;

    case 3:
      {
        mystate->n_gals = 0;
        mystate->gals = NULL;
      }
      break;

    default:
      SID_log_error("Fail!");
      break;
  }
  
  run_globals.FirstGal = mystate->gals;

  return 0;
}


static int teardown_map_galaxies_to_slabs(void **state) {
  SID_free(SID_FARG ((state_t *)*state)->gals);
  free_reionization_grids();
  SID_free(SID_FARG *state);

  return 0;
}


static void test_map_galaxies_to_slabs(void **state) {

  int correct_nix[] = {16, 16, 16, 16};
  for(int ii=0; ii<SID.n_proc; ii++)
    assert_int_equal((int)tocf_params.slab_nix[ii], correct_nix[ii]);

  int correct_ix_start[] = {0, 16, 32, 48};
  for(int ii=0; ii<SID.n_proc; ii++)
    assert_int_equal((int)tocf_params.slab_ix_start[ii], correct_ix_start[ii]);

  int n_mapped = map_galaxies_to_slabs(mystate->n_gals);
  assert_int_equal(n_mapped, mystate->n_gals);

  /* gal_to_slab_t *galaxy_to_slab_map = run_globals.tocf_grids.galaxy_to_slab_map; */
  
}


int main(int argc, char *argv[])
{
  SID_init(&argc, &argv, NULL, NULL);

  // Ensure we are running with the expected number of processors
  assert_int_equal(SID.n_proc, 4);

  const struct CMUnitTest tests[] = {
    cmocka_unit_test_setup_teardown(test_map_galaxies_to_slabs, setup_map_galaxies_to_slabs, teardown_map_galaxies_to_slabs),
  };

  int result = cmocka_run_group_tests(tests, NULL, NULL);

  SID_exit(result);
  return result;
}
