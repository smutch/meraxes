#include <assert.h>
#include <complex.h>
#include <fenv.h>
#include <fftw3-mpi.h>
#include <hdf5_hl.h>
#include <math.h>
#include <sys/stat.h>

#include "meraxes.h"
#include "misc_tools.h"
#include "virial_properties.h"
#include "metal_evo.h"

void assign_slabs_metals() // Not sure if I have to duplicate this function from reionization.c
{
  mlog("Assigning slabs to MPI cores for Metals...", MLOG_OPEN);

  // Assign the slab size
  int n_rank = run_globals.mpi_size;
  int dim = run_globals.params.MetalGridDim;

  // Use fftw to find out what slab each rank should get
  ptrdiff_t local_nix_metals, local_ix_start_metals;
  ptrdiff_t local_n_complex_metals =
    fftwf_mpi_local_size_3d(dim, dim, dim / 2 + 1, run_globals.mpi_comm, &local_nix_metals, &local_ix_start_metals);

  // let every core know...
  ptrdiff_t** slab_nix_metals = &run_globals.metal_grids.slab_nix_metals;
  *slab_nix_metals = malloc(sizeof(ptrdiff_t) * n_rank); ///< array of number of x cells of every rank
  MPI_Allgather(&local_nix_metals, sizeof(ptrdiff_t), MPI_BYTE, *slab_nix_metals, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm);

  ptrdiff_t** slab_ix_start_metals = &run_globals.metal_grids.slab_ix_start_metals;
  *slab_ix_start_metals = malloc(sizeof(ptrdiff_t) * n_rank); ///< array first x cell of every rank
  (*slab_ix_start_metals)[0] = 0;
  for (int ii = 1; ii < n_rank; ii++)
    (*slab_ix_start_metals)[ii] = (*slab_ix_start_metals)[ii - 1] + (*slab_nix_metals)[ii - 1];

  ptrdiff_t** slab_n_complex_metals = &run_globals.metal_grids.slab_n_complex_metals; ///< array of allocation counts for every rank
  *slab_n_complex_metals = malloc(sizeof(ptrdiff_t) * n_rank);                 ///< array of allocation counts for every rank
  MPI_Allgather(
    &local_n_complex_metals, sizeof(ptrdiff_t), MPI_BYTE, *slab_n_complex_metals, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm);

  mlog("...done", MLOG_CLOSE);
}

void construct_metal_grids(int snapshot, int local_ngals) //work in progress
{
  double box_size = run_globals.params.BoxSize; // Sei arrivato qua!
  
}
