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


static int setup_tocf_tests(void **state)
{
  *state = SID_malloc(sizeof(state_t));

  tocf_params.HII_dim = 64;
  tocf_params.uvb_feedback = 1;
  run_globals.params.TOCF_Flag = 1;
  run_globals.params.BoxSize = 100.;  // Not the size of Tiamat but easy for checking

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

          // for(int jj=0; jj<3; jj++)
          //   SID_log("R%d: [%d][%d] -> %.2f", SID_LOG_COMMENT|SID_LOG_ALLRANKS, SID.My_rank, ii, jj, gals[ii].Pos[jj]);
        }
      }
      break;

    case 1:
      {
        mystate->n_gals = 3;
        float xpos[] = {1, 2, 99};
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

          // for(int jj=0; jj<3; jj++)
          //   SID_log("R%d: [%d][%d] -> %.2f", SID_LOG_COMMENT|SID_LOG_ALLRANKS, SID.My_rank, ii, jj, gals[ii].Pos[jj]);
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

          // for(int jj=0; jj<3; jj++)
          //   SID_log("R%d: [%d][%d] -> %.2f", SID_LOG_COMMENT|SID_LOG_ALLRANKS, SID.My_rank, ii, jj, gals[ii].Pos[jj]);
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

  // grids
  tocf_grids_t *grids = &(run_globals.tocf_grids);
  int HII_dim = tocf_params.HII_dim;
  ptrdiff_t *slab_nix = tocf_params.slab_nix;
  ptrdiff_t slab_n_real = slab_nix[SID.My_rank] * HII_dim * HII_dim; // TODO: NOT WORKING?

  for (int ii = 0; ii < slab_n_real; ii++)
    grids->Mvir_crit[ii] = (float)(1000*SID.My_rank + ii);


  return 0;
}


static int teardown_tocf_tests(void **state)
{

  SID_free(SID_FARG run_globals.tocf_grids.galaxy_to_slab_map);
  SID_free(SID_FARG ((state_t *)*state)->gals);
  free_reionization_grids();
  SID_free(SID_FARG *state);

  return 0;
}


static void test_map_galaxies_to_slabs(void **state)
{

  int correct_nix[] = {16, 16, 16, 16};
  for(int ii=0; ii<SID.n_proc; ii++)
    assert_int_equal((int)tocf_params.slab_nix[ii], correct_nix[ii]);

  int correct_ix_start[] = {0, 16, 32, 48};
  for(int ii=0; ii<SID.n_proc; ii++)
    assert_int_equal((int)tocf_params.slab_ix_start[ii], correct_ix_start[ii]);

  int n_mapped = map_galaxies_to_slabs(mystate->n_gals);
  assert_int_equal(n_mapped, mystate->n_gals);

  gal_to_slab_t *galaxy_to_slab_map = run_globals.tocf_grids.galaxy_to_slab_map; 
  int *correct_slab;

  switch (SID.My_rank)
  {
    case 0:
      {
        // N.B. The galaxy_to_slab_map array is ordered and hence the
        // correct_slab values do not correspond to the order in which the
        // galaxies were created!
        correct_slab = SID_malloc(sizeof(int) * mystate->n_gals);
        correct_slab[0] = 0;
        correct_slab[1] = 1;
        correct_slab[2] = 2;
        correct_slab[3] = 3;
      }
      break;

    case 1:
      {
        correct_slab = SID_malloc(sizeof(int) * mystate->n_gals);
        correct_slab[0] = 0;
        correct_slab[1] = 0;
        correct_slab[2] = 3;
      }
      break;

    case 2:
      {
        correct_slab = SID_malloc(sizeof(int) * mystate->n_gals);
        correct_slab[0] = 2;
      }
      break;

    case 3:
      return;
      break;

    default:
      SID_log_error("Fail!");
  }

  for(int ii=0; ii<mystate->n_gals; ii++)
    assert_int_equal(galaxy_to_slab_map[ii].slab_ind, correct_slab[ii]);

  SID_free(SID_FARG correct_slab);
  
}


// static void test_assign_Mvir_crit_to_galaxies(void **state)
// {
  
//   switch (SID.My_rank)
//   {
//     case 0:
//       {
//         float correct_Mvir_crit[] = {};
//       }
//   }

// }


int main(int argc, char *argv[])
{
  SID_init(&argc, &argv, NULL, NULL);

  // Ensure we are running with the expected number of processors
  assert_int_equal(SID.n_proc, 4);

  const struct CMUnitTest tests[] = {
    cmocka_unit_test_setup_teardown(test_map_galaxies_to_slabs, setup_tocf_tests, teardown_tocf_tests),
  };

  int result = cmocka_run_group_tests(tests, NULL, NULL);

  SID_exit(result);
  return result;
}
