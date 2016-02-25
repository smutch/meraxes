#define _MAIN
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <meraxes.h>
#include <cmocka.h>

static int setup_map_galaxies_to_slabs(void **state) {
  (void) state; /* unused */

  tocf_params.HII_dim = 64;
  tocf_params.uvb_feedback = 1;
  run_globals.params.TOCF_Flag = 1;

  malloc_reionization_grids();

  return 0;
}

static int teardown_map_galaxies_to_slabs(void **state) {
  (void) state; /* unused */

  free_reionization_grids();

  return 0;
}


static void test_map_galaxies_to_slabs(void **state) {
  (void) state; /* unused */
}


int main(int argc, char *argv[])
{
  SID_init(&argc, &argv, NULL, NULL);

  const struct CMUnitTest tests[] = {
    cmocka_unit_test_setup_teardown(test_map_galaxies_to_slabs, setup_map_galaxies_to_slabs, teardown_map_galaxies_to_slabs),
  };

  int result = cmocka_run_group_tests(tests, NULL, NULL);

  SID_exit(result);
  return result;
}
