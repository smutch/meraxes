#define _MAIN
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <meraxes.h>
#include <cmocka.h>
#include <hdf5.h>
#include <hdf5_hl.h>

typedef struct state_t {
  float redshift;
  int snapshot;
} state_t;
#define mystate ((state_t *)*state)

#define MERAXES_TEST_OUTPUT_DIR "/mnt/home/smutch/models/21cm_sam/meraxes/src/tests/output/"

static int setup(void **state)
{
  *state = SID_malloc(sizeof(state_t));

  mystate->redshift = 5.0;
  mystate->snapshot = 100;

  char param_file[] = "/home/smutch/models/21cm_sam/meraxes/src/tests/input/dummy.par";
  read_parameter_file(param_file, 0);

  run_globals.params.BoxSize     = 100.;
  run_globals.params.TOCF_Flag   = 1;
  tocf_params.HII_dim            = 32;
  tocf_params.HII_eff_factor     = 1.0;
  tocf_params.uvb_feedback       = 1;
  tocf_params.r_bubble_min       = 0.5;
  tocf_params.r_bubble_max       = 5.;
  tocf_params.delta_r_HII_factor = 1.1;
  tocf_params.RtoM_filter        = 0;
  tocf_params.HII_filter         = 0;

  set_units();
  malloc_reionization_grids();

  int HII_dim   = tocf_params.HII_dim;
  float *deltax = run_globals.tocf_grids.deltax;
  float *stars  = run_globals.tocf_grids.stars;
  float *sfr    = run_globals.tocf_grids.sfr;
  int n_ix = tocf_params.slab_nix[SID.My_rank];
  for(int ii=0; ii<n_ix; ii++)
    for(int jj=0; jj<HII_dim; jj++)
      for(int kk=0; kk<HII_dim; kk++)
      {
        deltax[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)] = 0.0;
        stars[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)] = 0.0;
        sfr[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)] = 0.0;
      }

  // This should be the central cell of the box
  if(SID.My_rank == SID.n_proc/2)
  {
    stars[grid_index(0, HII_dim/2, HII_dim/2, HII_dim, INDEX_PADDED)] = 1.0e6;
    sfr[grid_index(0, HII_dim/2, HII_dim/2, HII_dim, INDEX_PADDED)] = 100.;
  }

  return 0;

}


static void test_setup(void **state)
{
  int HII_dim   = tocf_params.HII_dim;
  int n_ix    = tocf_params.slab_nix[SID.My_rank];

  char fname[128];
  sprintf(fname, MERAXES_TEST_OUTPUT_DIR"/test_setup-%d.h5", SID.My_rank);
  hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  float *temp = SID_malloc(n_ix * HII_dim * HII_dim * sizeof(float));

  for(int ii=0; ii<n_ix; ii++)
    for(int jj=0; jj<HII_dim; jj++)
      for(int kk=0; kk<HII_dim; kk++)
        temp[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)] = 
          run_globals.tocf_grids.stars[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)];
  H5LTmake_dataset_float(fid, "stars", 3, (hsize_t[3]){n_ix, HII_dim, HII_dim}, run_globals.tocf_grids.stars);

  for(int ii=0; ii<n_ix; ii++)
    for(int jj=0; jj<HII_dim; jj++)
      for(int kk=0; kk<HII_dim; kk++)
        temp[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)] = 
          run_globals.tocf_grids.sfr[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)];
  H5LTmake_dataset_float(fid, "sfr", 3, (hsize_t[3]){n_ix, HII_dim, HII_dim}, run_globals.tocf_grids.sfr);

  H5Fclose(fid);
  SID_free(SID_FARG temp);
}


static void test_find_HII_bubbles(void **state)
{
  find_HII_bubbles(mystate->redshift);

  int HII_dim = tocf_params.HII_dim;
  float *xH = run_globals.tocf_grids.xH;
  int n_ix    = tocf_params.slab_nix[SID.My_rank];

  // if (SID.My_rank == SID.n_proc/2)
  // {
  //   assert_int_equal(xH[grid_index(0, HII_dim/2, HII_dim/2, HII_dim, INDEX_REAL)], 0);
  // }

  char fname[512];
  sprintf(fname, MERAXES_TEST_OUTPUT_DIR"/test_find_bubbles-%d.h5", SID.My_rank);
  hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  H5LTmake_dataset_float(fid, "xH", 3, (hsize_t[3]){n_ix, HII_dim, HII_dim}, xH);
  H5Fclose(fid);
}


static int teardown(void **state)
{
  free_reionization_grids();
  SID_free(SID_FARG *state);

  return 0;
}


int main(int argc, char *argv[])
{
  SID_init(&argc, &argv, NULL, NULL);

  // Turn of SID logging
  // FILE *fp_log = fopen("/dev/null", "w");
  // SID.fp_log = fp_log;


  // Ensure we are running with the expected number of processors
  assert_int_equal(SID.n_proc, 4);

  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_setup),
    cmocka_unit_test(test_find_HII_bubbles),
  };

  int result = cmocka_run_group_tests(tests, setup, teardown);

  SID_exit(result);
  // fclose(fp_log);
  return result;
}
