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
#include "reionization.c" //For write_grid_float and other stuff

void assign_slabs_metals() // Not sure if I have to duplicate this function from reionization.c
{
  mlog("Assigning slabs to MPI cores for Metals...", MLOG_OPEN);

  // Assign the slab size
  int n_rank = run_globals.mpi_size;
  int dim = run_globals.params.MetalGridDim;

  // Use fftw to find out what slab each rank should get
  ptrdiff_t local_nix_metals, local_ix_start_metals;
  fftwf_mpi_local_size_3d(dim, dim, dim / 2 + 1, run_globals.mpi_comm, &local_nix_metals, &local_ix_start_metals); // Important to define local_nix and local_ix_start

  // let every core know...
  ptrdiff_t** slab_nix_metals = &run_globals.metal_grids.slab_nix_metals;
  *slab_nix_metals = malloc(sizeof(ptrdiff_t) * n_rank); ///< array of number of x cells of every rank
  MPI_Allgather(&local_nix_metals, sizeof(ptrdiff_t), MPI_BYTE, *slab_nix_metals, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm); 

  ptrdiff_t** slab_ix_start_metals = &run_globals.metal_grids.slab_ix_start_metals;
  *slab_ix_start_metals = malloc(sizeof(ptrdiff_t) * n_rank); ///< array first x cell of every rank
  (*slab_ix_start_metals)[0] = 0;
  for (int ii = 1; ii < n_rank; ii++)
    (*slab_ix_start_metals)[ii] = (*slab_ix_start_metals)[ii - 1] + (*slab_nix_metals)[ii - 1];

  mlog("...done", MLOG_CLOSE);
}

void construct_metal_grids(int snapshot, int local_ngals) // You can put here the computation of probability
{
  double box_size = run_globals.params.BoxSize; 
  float* stellar_grid_metals = run_globals.metal_grids.stars_metals;
  float* sfr_grid_metals = run_globals.metal_grids.sfr_metals;
  float* mass_metals_grid_metals = run_globals.metal_grids.mass_metals;
  float* mass_gas_grid_metals = run_globals.metal_grids.mass_gas;
  int MetalGridDim = run_globals.params.MetalGridDim;
  double sfr_timescale_metals = run_globals.params.ReionSfrTimescale * hubble_time(snapshot); 
  
  gal_to_slab_t* galaxy_to_slab_map_metals = run_globals.metal_grids.galaxy_to_slab_map_metals;
  ptrdiff_t* slab_ix_start_metals = run_globals.metal_grids.slab_ix_start_metals;
  ptrdiff_t* slab_nix_metals = run_globals.metal_grids.slab_nix_metals;
  ptrdiff_t slab_n_real_metals = slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim;
  //int local_nix_metals = (int)(run_globals.metal_grids.slab_nix_metals[run_globals.mpi_rank]);

  mlog("Constructing stellar mass and sfr grids for metals...", MLOG_OPEN | MLOG_TIMERSTART);

  // init the grid
  //for (int ii = 0; ii < local_n_complex_metals * 2; ii++) {
  for (int ii = 0; ii < slab_n_real_metals; ii++) { 
    stellar_grid_metals[ii] = 0.0;
    sfr_grid_metals[ii] = 0.0;
    mass_metals_grid_metals[ii] = 0.0;
    mass_gas_grid_metals[ii] = 0.0;
  }

  // loop through each slab
  //
  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  //ptrdiff_t* slab_nix_metals = run_globals.metal_grids.slab_nix_metals;
  ptrdiff_t buffer_size_metals = run_globals.metal_grids.buffer_size_metals;
  float* buffer_metals = run_globals.metal_grids.buffer_metals;

  enum property
  {
    prop_stellar_metals,
    prop_sfr_metals,
    prop_mass_ej_metals,
    prop_mass_ej_gas
  };
  
  for (int prop = prop_stellar_metals; prop <= prop_mas_ej_gas; prop++) {

    int i_gal = 0;
    int skipped_gals = 0;

    for (int i_r = 0; i_r < run_globals.mpi_size; i_r++) {
      // init the buffer
      for (int ii = 0; ii < buffer_size_metals; ii++)
        buffer_metals[ii] = (float)0.;

      // if this core holds no galaxies then we don't need to fill the buffer
      if (local_ngals != 0)
        // fill the local buffer for this slab
        while (((i_gal - skipped_gals) < local_ngals) && (galaxy_to_slab_map_metals[i_gal].slab_ind == i_r)) {
          galaxy_t* gal = galaxy_to_slab_map_metals[i_gal].galaxy;

          // Dead galaxies should not be included here and are not in the
          // local_ngals count.  They will, however, have been assigned to a
          // slab so we will need to ignore them here...
          if (gal->Type > 2) {
            i_gal++;
            skipped_gals++;
            continue;
          }

          assert(galaxy_to_slab_map_metals[i_gal].index >= 0);
          assert((galaxy_to_slab_map_metals[i_gal].slab_ind >= 0) &&
                 (galaxy_to_slab_map_metals[i_gal].slab_ind < run_globals.mpi_size));

          int ix = (int)(pos_to_ngp(gal->Pos[0], box_size, MetalGridDim) - slab_ix_start_metals[i_r]);
          int iy = pos_to_ngp(gal->Pos[1], box_size, MetalGridDim);
          int iz = pos_to_ngp(gal->Pos[2], box_size, MetalGridDim);

          assert((ix < slab_nix_metals[i_r]) && (ix >= 0));
          assert((iy < MetalGridDim) && (iy >= 0));
          assert((iz < MetalGridDim) && (iz >= 0));

          int ind = grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL); // Check what this function does

          assert((ind >= 0) && (ind < slab_nix_metals[i_r] * MetalGridDim * MetalGridDim));
          
          switch (prop) {
            case prop_stellar_metals: // Not sure if you really need this for metals

              buffer_metals[ind] += gal->GrossStellarMass;
              break;


            case prop_sfr_metals:
            
              buffer_metals[ind] += gal->GrossStellarMass;
			  // this sfr grid is used for X-ray and Lyman, we are ignoring BH here, in principle I should use this for metals.
              break;
              
            case prop_mass_ej_metals:
            
              buffer_metals[ind] += gal->MetalsEjectedGas;
              break;
              
            case prop_mass_ej_gas:
            
              buffer_metals[ind] += gal->EjectedGas;
              break; 

            default:
              mlog_error("Unrecognised property in slab creation.");
              ABORT(EXIT_FAILURE);
              break;
          }

          i_gal++;
        }

      // reduce on to the correct rank
      if (run_globals.mpi_rank == i_r)
        MPI_Reduce(MPI_IN_PLACE, buffer_metals, (int)buffer_size_metals, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);
      else
        MPI_Reduce(buffer_metals, buffer_metals, (int)buffer_size_metals, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);

      if (run_globals.mpi_rank == i_r)

        // Do one final pass and divide the sfr_grid by the sfr timescale
        // in order to convert the stellar masses recorded into SFRs before
        // finally copying the values into the appropriate slab.
        switch (prop) {
         
          case prop_sfr_metals:
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  val = (val > 0) ? val / sfr_timescale_metals : 0; 
                  sfr_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = (float)val;
                }
            break;
            
          case prop_mass_ej_metals: // Need this to compute average metallicity of the cell
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  if (val < 0)
                    val = 0;
                  mass_metals_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = (float)val;
                }
            break;
            
          case prop_mass_ej_gas: // Need this to compute average metallicity of the cell
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  if (val < 0)
                    val = 0;
                  mass_gas_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = (float)val;
                }
            break;

          case prop_stellar_metals:
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  float val = buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  if (val < 0)
                    val = 0;
                  stellar_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = val;
                }
            break;

          default:
            mlog_error("Eh!?!");
            ABORT(EXIT_FAILURE);
        }
    }
  }

  mlog("done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

void gen_metal_grids_fname(const int snapshot, char* name, const bool relative) //In the future we can merge this function
{
  if (!relative)
    sprintf(name, "%s/%s_metal_grids_%d.hdf5", run_globals.params.OutputDir, run_globals.params.FileNameGalaxies, snapshot);
  else
    sprintf(name, "%s_metal_grids_%d.hdf5", run_globals.params.FileNameGalaxies, snapshot);
}

void save_metal_input_grids(int snapshot)
{
  metal_grids_t* grids = &(run_globals.metal_grids);
  int MetalGridDim = run_globals.params.MetalGridDim;
  int local_nix_metals = (int)(run_globals.metal_grids.slab_nix_metals[run_globals.mpi_rank]);
  double UnitTime_in_s = run_globals.units.UnitTime_in_s;
  double UnitMass_in_g = run_globals.units.UnitMass_in_g;

  mlog("Saving tocf input metal grids...", MLOG_OPEN);

  char name[STRLEN];
  gen_metal_grids_fname(snapshot, name, false);

  // create the file (in parallel)
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // create the filespace
  hsize_t dims[3] = { (hsize_t)MetalGridDim, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim };
  hid_t fspace_id = H5Screate_simple(3, dims, NULL);

  // create the memspace
  hsize_t mem_dims[3] = { (hsize_t)local_nix_metals, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim };
  hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

  // select a hyperslab in the filespace
  hsize_t start[3] = { (hsize_t)run_globals.metal_grids.slab_ix_start_metals[run_globals.mpi_rank], 0, 0 };
  hsize_t count[3] = { (hsize_t)local_nix_metals, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim };
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  // set the dataset creation property list to use chunking along x-axis
  hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl_id, 3, (hsize_t[3]){ 1, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim });

  // fftw padded grids
  float* grid = (float*)calloc((size_t)local_nix_metals * (size_t)MetalGridDim * (size_t)MetalGridDim, sizeof(float));


  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (grids->stars_metals)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)];
  write_grid_float("stars_metals", grid, file_id, fspace_id, memspace_id, dcpl_id);

  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (float)((grids->sfr_metals)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  write_grid_float("sfr_metals", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (float)((grids->mass_metals)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] * UnitMass_in_g / SOLAR_MASS);
  write_grid_float("mass_metals", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (float)((grids->mass_gas)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] * UnitMass_in_g / SOLAR_MASS);
  write_grid_float("mass_gas", grid, file_id, fspace_id, memspace_id, dcpl_id);


  // tidy up
  free(grid);
  H5Pclose(dcpl_id);
  H5Sclose(memspace_id);
  H5Sclose(fspace_id);
  H5Fclose(file_id);

  mlog("...done", MLOG_CLOSE);
}

void init_metal_grids()
{
  metal_grids_t* grids = &(run_globals.metal_grids);
  int MetalGridDim = run_globals.params.MetalGridDim;
  ptrdiff_t* slab_nix_metals = run_globals.metal_grids.slab_nix_metals;
  ptrdiff_t slab_n_real_metals = slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim;

  mlog("Initialising metal grids...", MLOG_MESG);

  grids->volume_ave_ZIGM = 0.0;
  grids->volume_ave_mass_metals = 0.0;
  grids->Probability_metals = 0.0;

  for (int ii = 0; ii < slab_n_real_metals; ii++) {
    grids->stars_metals[ii] = 0;
    grids->sfr_metals[ii] = 0;
    grids->mass_metals[ii] = (float)0.;
    grids->mass_gas[ii] = (float)0.;
    grids->Zigm_box[ii] = (float)0.;
  }
}

void malloc_metal_grids()
{
  mlog("Allocating reionization grids...", MLOG_OPEN);

  metal_grids_t* grids = &(run_globals.metal_grids);

  // run_globals.NStoreSnapshots is set in `initialize_halo_storage`
  //run_globals.SnapshotDeltax = (float**)calloc((size_t)run_globals.NStoreSnapshots, sizeof(float*)); //?
  //run_globals.SnapshotVel = (float**)calloc((size_t)run_globals.NStoreSnapshots, sizeof(float*));  //?

  grids->galaxy_to_slab_map_metals = NULL;

  grids->stars_metals = NULL;
  grids->sfr_metals = NULL;
  grids->mass_metals = NULL;
  grids->mass_gas = NULL;
  grids->Zigm_box = NULL; 
  
  assign_slabs_metals();

    int MetalGridDim = run_globals.params.MetalGridDim;
    ptrdiff_t* slab_nix_metals = run_globals.metal_grids.slab_nix_metals;
    ptrdiff_t slab_n_real_metals = slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim;
    // create a buffer on each rank which is as large as the largest LOGICAL allocation on any single rank
    int max_cells = 0;

    for (int ii = 0; ii < run_globals.mpi_size; ii++)
      if (slab_nix_metals[ii] > max_cells)
        max_cells = (int)slab_nix_metals[ii];

    max_cells *= MetalGridDim * MetalGridDim;
    grids->buffer_size_metals = max_cells;

    grids->buffer_metals = fftwf_alloc_real((size_t)max_cells);

    grids->stars_metals = fftwf_alloc_real((size_t)slab_n_real_metals);



    grids->mass_metals = fftwf_alloc_real((size_t)slab_n_real_metals);
    grids->mass_gas = fftwf_alloc_real((size_t)slab_n_real_metals);
    grids->Zigm_box = fftwf_alloc_real((size_t)slab_n_real_metals);

    grids->sfr_metals = fftwf_alloc_real((size_t)slab_n_real_metals);

    init_metal_grids();

  mlog("...done", MLOG_CLOSE);
}

void free_metal_grids()
{
  mlog("Freeing metal grids...", MLOG_OPEN);

  metal_grids_t* grids = &(run_globals.metal_grids);

  free(run_globals.metal_grids.slab_ix_start_metals);
  free(run_globals.metal_grids.slab_nix_metals);

  //fftwf_free(grids->volume_ave_ZIGM);
  //fftwf_free(grids->volume_ave_mass_metals);
  fftwf_free(grids->mass_gas);
  fftwf_free(grids->mass_metals);
  fftwf_free(grids->Zigm_box);

  fftwf_free(grids->sfr_metals);

  fftwf_free(grids->stars_metals);

  fftwf_free(grids->buffer_metals);

  mlog(" ...done", MLOG_CLOSE);
}

int map_galaxies_to_slabs_metals(int ngals)
{
  double box_size = run_globals.params.BoxSize;
  int MetalGridDim = run_globals.params.MetalGridDim;

  mlog("Mapping galaxies to slabs for metals...", MLOG_OPEN);

  // Loop through each valid galaxy and find what slab it sits in
  if (ngals > 0)
    run_globals.metal_grids.galaxy_to_slab_map_metals = malloc(sizeof(gal_to_slab_metals_t) * ngals);
  else
    run_globals.metal_grids.galaxy_to_slab_map_metals = NULL;

  gal_to_slab_metals_t* galaxy_to_slab_map_metals = run_globals.metal_grids.galaxy_to_slab_map_metals;
  ptrdiff_t* slab_ix_start_metals = run_globals.metal_grids.slab_ix_start_metals;

  galaxy_t* gal = run_globals.FirstGal;
  int gal_counter = 0;
  while (gal != NULL) {
    // TODO: Note that I am including ghosts here.  We will need to check the
    // validity of this.  By definition, if they are ghosts then their host
    // halo hasn't been identified at this time step and hence they haven't
    // been evolved.  Their properties (Sfr, StellarMass, etc.) will all have
    // been set when they were last identified.
    if (gal->Type < 3) {
      // TODO: for type 2 galaxies these positions will be set from the last
      // time they were identified.  If their host halo has moved significantly
      // since then, these positions won't reflect that and the satellites will
      // be spatially disconnected from their hosts.  We will need to fix this
      // at some point.

      ptrdiff_t ix = pos_to_ngp(gal->Pos[0], box_size, MetalGridDim);

      assert((ix >= 0) && (ix < MetalGridDim));

      galaxy_to_slab_map_metals[gal_counter].index = gal_counter;
      galaxy_to_slab_map_metals[gal_counter].slab_ind =
        searchsorted(&ix, slab_ix_start_metals, run_globals.mpi_size, sizeof(ptrdiff_t), compare_ptrdiff, -1, -1);
      galaxy_to_slab_map_metals[gal_counter++].galaxy = gal;
    }

    gal = gal->Next;
  }

  // sort the slab indices IN PLACE (n.b. compare_slab_assign is a stable comparison)
  if (galaxy_to_slab_map_metals != NULL)
    qsort(galaxy_to_slab_map_metals, (size_t)gal_counter, sizeof(gal_to_slab_t), compare_slab_assign);

  assert(gal_counter == ngals);

  mlog("...done.", MLOG_CLOSE);

  return gal_counter;
}

void save_metal_output_grids(int snapshot)
{

  metal_grids_t* metal_grids = &(run_globals.metal_grids);
  int MetalGridDim = run_globals.params.MetalGridDim;
  int local_nix_metals = (int)(run_globals.metal_grids.slab_nix_metals[run_globals.mpi_rank]);


  mlog("Saving tocf output grids for metals...", MLOG_OPEN);

  char name[STRLEN];
  gen_metal_grids_fname(snapshot, name, false);

  // open the file (in parallel)
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(name, H5F_ACC_RDWR, plist_id);
  H5Pclose(plist_id);

  // create the filespace
  hsize_t dims[3] = { (hsize_t)MetalGridDim, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim };
  hid_t fspace_id = H5Screate_simple(3, dims, NULL);

  // create the memspace
  hsize_t mem_dims[3] = { (hsize_t)local_nix_metals, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim };
  hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

  // select a hyperslab in the filespace
  hsize_t start[3] = { (hsize_t)run_globals.metal_grids.slab_ix_start_metals[run_globals.mpi_rank], 0, 0 };
  hsize_t count[3] = { (hsize_t)local_nix_metals, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim };
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

  // set the dataset creation property list to use chunking along x-axis
  hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl_id, 3, (hsize_t[3]){ 1, (hsize_t)MetalGridDim, (hsize_t)MetalGridDim });

  // create and write the datasets
  write_grid_float("mass_metals", metal_grids->mass_metals, file_id, fspace_id, memspace_id, dcpl_id);
  write_grid_float("Zigm_box", metal_grids->Zigm_box, file_id, fspace_id, memspace_id, dcpl_id);
  //write_grid_float("Probability_metals", grids->Probability_metals, file_id, fspace_id, memspace_id, dcpl_id);

  // fftw padded grids
  float* grid = (float*)calloc((size_t)(local_nix_metals * MetalGridDim * MetalGridDim), sizeof(float));
  
  H5LTset_attribute_double(file_id, "Zigm_box", "volume_ave_ZIGM", &(metal_grids->volume_ave_ZIGM), 1);
  H5LTset_attribute_double(file_id, "mass_metals", "volume_ave_mass_metals", &(metal_grids->volume_ave_mass_metals), 1);
    

  // tidy up
  free(grid);
  H5Pclose(dcpl_id);
  H5Sclose(memspace_id);
  H5Sclose(fspace_id);
  H5Fclose(file_id);

  mlog("...done", MLOG_CLOSE); // Saving tocf grids
}
