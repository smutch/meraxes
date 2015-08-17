#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>

void read_Mcrit_table(run_globals_t *run_globals)
{

  // TODO: Make sure Flag_ReionizationModifier gets set differently if using uvb fb and 21cmFAST
  if (run_globals->params.physics.Flag_ReionizationModifier != 3)
  {
    // Set this to NULL if we don't need it...
    run_globals->params.MvirCrit = NULL;
    return;
  }

  if(SID.My_rank == 0)
  {
    hid_t fd;
    hsize_t dims;

    // open the file
    fd = H5Fopen(run_globals->params.MvirCritFile, H5P_DEFAULT, H5P_DEFAULT);
    assert(fd >= 0);

    // read the dataset size and ensure it is of the correct size
    H5LTget_dataset_info(fd, "mean_Mvir_crit", &dims, NULL, NULL);
    assert((int)dims == run_globals->params.SnaplistLength);

    // read the dataset
    run_globals->params.MvirCrit = SID_malloc(sizeof(double) * (int)dims);
    H5LTread_dataset_double(fd, "mean_Mvir_crit", run_globals->params.MvirCrit);

    // close the file
    H5Fclose(fd);
  }
  else
    run_globals->params.MvirCrit = SID_malloc(sizeof(double) * run_globals->params.SnaplistLength);

  // broadcast the result to the other ranks
  SID_Bcast(run_globals->params.MvirCrit, run_globals->params.SnaplistLength, 0, SID.COMM_WORLD);

}
