#define _MAIN
#include <criterion/criterion.h>
#include <criterion/parameterized.h>
#include <meraxes.h>
#include <mpi.h>
#include "../core/parse_paramfile.h"

void setup(void)
{
  int argc = 0;
  char** argv = NULL;

  MPI_Init(&argc, &argv);
  run_globals.mpi_comm = MPI_COMM_WORLD;
  run_globals.mpi_rank = 0;
  run_globals.mpi_size = 1;

  run_globals.params.SnaplistLength = 11;
}

void teardown(void)
{
  MPI_Finalize();
}

TestSuite(parse_snaplist, .init = setup, .fini = teardown);

Test(parse_snaplist, spaces)
{
  int expected[11] = { 0 };

  parse_output_snaps("1, 2, 3");
  memcpy(expected, (int[]){ 1, 2, 3 }, sizeof(int) * 3);
  cr_expect_eq(run_globals.NOutputSnaps, 3);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 3);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("1,2, 3");
  memcpy(expected, (int[]){ 1, 2, 3 }, sizeof(int) * 3);
  cr_expect_eq(run_globals.NOutputSnaps, 3);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 3);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("1,2, 3");
  memcpy(expected, (int[]){ 1, 2, 3 }, sizeof(int) * 3);
  cr_expect_eq(run_globals.NOutputSnaps, 3);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 3);
  free(run_globals.ListOutputSnaps);
}

Test(parse_snaplist, negative)
{
  int expected[11] = { 0 };

  parse_output_snaps("-1");
  memcpy(expected, (int[]){ 10 }, sizeof(int) * 1);
  cr_expect_eq(run_globals.NOutputSnaps, 1);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 1);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("-1, -3:");
  memcpy(expected, (int[]){ 8, 9, 10 }, sizeof(int) * 3);
  cr_expect_eq(run_globals.NOutputSnaps, 3);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 3);
  free(run_globals.ListOutputSnaps);
}

Test(parse_snaplist, slices)
{
  int expected[11] = { 0 };

  parse_output_snaps("1, 2:6, 6");
  memcpy(expected, (int[]){ 1, 2, 3, 4, 5, 6 }, sizeof(int) * 6);
  cr_expect_eq(run_globals.NOutputSnaps, 6);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 6);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("1, 2:6, 5");
  memcpy(expected, (int[]){ 1, 2, 3, 4, 5 }, sizeof(int) * 5);
  cr_expect_eq(run_globals.NOutputSnaps, 5);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 5);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("1, 2:-1, 6");
  memcpy(expected, (int[]){ 1, 2, 3, 4, 5, 6, 7, 8, 9 }, sizeof(int) * 9);
  cr_expect_eq(run_globals.NOutputSnaps, 9);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 9);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("1, 2:");
  memcpy(expected, (int[]){ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, sizeof(int) * 10);
  cr_expect_eq(run_globals.NOutputSnaps, 10);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 10);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps(":");
  memcpy(expected, (int[]){ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, sizeof(int) * 11);
  cr_expect_eq(run_globals.NOutputSnaps, 11);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 11);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps(":8");
  memcpy(expected, (int[]){ 0, 1, 2, 3, 4, 5, 6, 7 }, sizeof(int) * 8);
  cr_expect_eq(run_globals.NOutputSnaps, 8);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 8);
  free(run_globals.ListOutputSnaps);

  parse_output_snaps("0:9");
  memcpy(expected, (int[]){ 0, 1, 2, 3, 4, 5, 6, 7, 8 }, sizeof(int) * 9);
  cr_expect_eq(run_globals.NOutputSnaps, 9);
  cr_expect_arr_eq(run_globals.ListOutputSnaps, expected, sizeof(int) * 9);
  free(run_globals.ListOutputSnaps);
}
