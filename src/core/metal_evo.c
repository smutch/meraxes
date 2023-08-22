#include <assert.h>
#include <complex.h>
#include <fenv.h>
#include <fftw3-mpi.h>
#include <hdf5_hl.h>
#include <math.h>
#include <sys/stat.h>

#if USE_MINI_HALOS
#include "meraxes.h"
#include "misc_tools.h"
#include "virial_properties.h"
#include "metal_evo.h"
#include "reionization.c" 

void assign_slabs_metals() 
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

void construct_metal_grids(int snapshot, int local_ngals) 
{
  double box_size = run_globals.params.BoxSize; 
  
  float* mass_metals_grid_metals = run_globals.metal_grids.mass_metals;
  float* mass_gas_grid_metals = run_globals.metal_grids.mass_gas;
  float* prob_grid_metals = run_globals.metal_grids.Probability_metals;
  float* Rave_grid_metals = run_globals.metal_grids.R_ave; 
  float* Rmax_grid_metals = run_globals.metal_grids.R_max;
  int* count_bubble_metals = run_globals.metal_grids.N_bubbles;
  
  int MetalGridDim = run_globals.params.MetalGridDim;
  double redshift = run_globals.ZZ[snapshot]; 
  double pixel_volume_metals = pow(box_size / (double)MetalGridDim, 3); // (Mpc/h)^3
  
  gal_to_slab_t* galaxy_to_slab_map_metals = run_globals.metal_grids.galaxy_to_slab_map_metals;
  ptrdiff_t* slab_ix_start_metals = run_globals.metal_grids.slab_ix_start_metals;
  ptrdiff_t* slab_nix_metals = run_globals.metal_grids.slab_nix_metals;
  ptrdiff_t slab_n_real_metals = slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim;

  mlog("Constructing stellar mass and sfr grids for metals...", MLOG_OPEN | MLOG_TIMERSTART);

  // init the grid
  for (int ii = 0; ii < slab_n_real_metals; ii++) { 
    mass_metals_grid_metals[ii] = 0.0;
    mass_gas_grid_metals[ii] = 0.0;
    prob_grid_metals[ii] = 0.0;
    count_bubble_metals[ii] = 0;
    Rave_grid_metals[ii] = 0.0;
    Rmax_grid_metals[ii] = 0.0;
  }

  // loop through each slab
  ptrdiff_t buffer_size_metals = run_globals.metal_grids.buffer_size_metals;
  float* buffer_metals = run_globals.metal_grids.buffer_metals;

  enum property
  {
    prop_prob,
    prop_count,
    prop_Rave,
    prop_Rmax,
    prop_mass_ej_metals,
    prop_mass_ej_gas
  };
  
  for (int prop = prop_prob; prop <= prop_mass_ej_gas; prop++) {

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

          int ind = grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL); 

          assert((ind >= 0) && (ind < slab_nix_metals[i_r] * MetalGridDim * MetalGridDim));
          
          switch (prop) {
            case prop_prob:
            
              if (gal->RmetalBubble >= 3 * gal->Rvir) {  // A bubble can actually pollute the IGM only if it's bigger than its virial radius (Try with 3 atm)
                buffer_metals[ind] += (4.0 / 3.0 * M_PI * pow((gal->RmetalBubble) * (1 + redshift), 3.0)); // cMpc/h (same units of cell volume)
                  
                  if (gal->RmetalBubble > 0.62 * (box_size / MetalGridDim)) { 
                  double Excess_volume = (4.0 / 3.0 * M_PI * pow((gal->RmetalBubble) * (1 + redshift), 3.0)) - pow((box_size / MetalGridDim), 3.0);
                  
                  //Adiacent cells in the same axis 
                  buffer_metals[ind + 1] += 0.072 * Excess_volume;
                  buffer_metals[ind - 1] += 0.072 * Excess_volume;
                  buffer_metals[ind + MetalGridDim] += 0.072 * Excess_volume;
                  buffer_metals[ind - MetalGridDim] += 0.072 * Excess_volume;
                  buffer_metals[ind + (MetalGridDim * MetalGridDim)] += 0.072 * Excess_volume;
                  buffer_metals[ind - (MetalGridDim * MetalGridDim)] += 0.072 * Excess_volume;
                  
                  //Adiacent cells obliquos
                  buffer_metals[ind + 1 + MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind - 1 + MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind - 1 - MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind + 1 - MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind + (MetalGridDim * MetalGridDim) + 1] += 0.0473 * Excess_volume;
                  buffer_metals[ind - (MetalGridDim * MetalGridDim) + 1] += 0.0473 * Excess_volume;
                  buffer_metals[ind + (MetalGridDim * MetalGridDim) - 1] += 0.0473 * Excess_volume;
                  buffer_metals[ind - (MetalGridDim * MetalGridDim) - 1] += 0.0473 * Excess_volume;
                  buffer_metals[ind + (MetalGridDim * MetalGridDim) + MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind - (MetalGridDim * MetalGridDim) + MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind + (MetalGridDim * MetalGridDim) - MetalGridDim] += 0.0473 * Excess_volume;
                  buffer_metals[ind - (MetalGridDim * MetalGridDim) - MetalGridDim] += 0.0473 * Excess_volume;
                }
              }

              break;
              
            case prop_Rave:
            
              if (gal->RmetalBubble > 0.)
                buffer_metals[ind] += gal->RmetalBubble * (1 + redshift); // cMpc/h
              
              break;
              
            case prop_Rmax:
              
              if (gal->RmetalBubble * (1 + redshift) >= buffer_metals[ind])
                buffer_metals[ind] = gal->RmetalBubble * (1 + redshift); //cMpc/h  
              
              break;
                            
            case prop_count:
            
              if (gal->RmetalBubble > 0.)
                buffer_metals[ind] += 1;

              break;
              
            case prop_mass_ej_metals:
            
              buffer_metals[ind] += gal->MetalsEjectedGas; //Internal units (same of gas_cell)
              break;
              
            case prop_mass_ej_gas:
            
              buffer_metals[ind] += (gal->EjectedGas - gal->HotGas - gal->ColdGas);
              break; 

            default:
              mlog_error("Unrecognised property in slab creation.");
              ABORT(EXIT_FAILURE);
              break;
          }

          i_gal++;
        }

      // reduce on to the correct rank
      if (prop == prop_Rmax) {
        if (run_globals.mpi_rank == i_r)
          MPI_Reduce(MPI_IN_PLACE, buffer_metals, (int)buffer_size_metals, MPI_FLOAT, MPI_MAX, i_r, run_globals.mpi_comm);
        else
          MPI_Reduce(buffer_metals, buffer_metals, (int)buffer_size_metals, MPI_FLOAT, MPI_MAX, i_r, run_globals.mpi_comm);
      }
      
      else {
        if (run_globals.mpi_rank == i_r)
          MPI_Reduce(MPI_IN_PLACE, buffer_metals, (int)buffer_size_metals, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);
        else
          MPI_Reduce(buffer_metals, buffer_metals, (int)buffer_size_metals, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);
      }

      if (run_globals.mpi_rank == i_r)

        // Do one final pass and divide the Rave_grid by the number of bubbles
        // finally copying the values into the appropriate slab.
        switch (prop) {
         
          case prop_prob:
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] / pixel_volume_metals; // You want this comoving
                  if (val < 0)
                    val = 0;
                  if (val > 1) {
                  // If R > cbrt(3/4*pi)*(L^3) the bubble is large enough to contribute to the filling factor of the nearby cells
                  // These factors assume that the large bubble originates at the centre of the cell and that Rbubble is smaller than 1.5 * (box_size/MetalGridDim).
                    /*double ExcessVal = val - 1.0;
                    
                    //Adiacent cells in the same axis 
                    prob_grid_metals[grid_index(ix + 1, iy, iz, MetalGridDim, INDEX_REAL)] += (float)(0.072 * ExcessVal);
                    prob_grid_metals[grid_index(ix - 1, iy, iz, MetalGridDim, INDEX_REAL)] += (float)(0.072 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy + 1, iz, MetalGridDim, INDEX_REAL)] += (float)(0.072 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy - 1, iz, MetalGridDim, INDEX_REAL)] += (float)(0.072 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy, iz + 1, MetalGridDim, INDEX_REAL)] += (float)(0.072 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy, iz - 1, MetalGridDim, INDEX_REAL)] += (float)(0.072 * ExcessVal);
                    
                    //Adiacent cells obliquos
                    prob_grid_metals[grid_index(ix + 1, iy + 1, iz, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix - 1, iy + 1, iz, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix - 1, iy - 1, iz, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix + 1, iy - 1, iz, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy + 1, iz + 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy + 1, iz - 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix + 1, iy, iz + 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix - 1, iy, iz - 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy - 1, iz + 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix, iy - 1, iz - 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix - 1, iy, iz + 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);
                    prob_grid_metals[grid_index(ix + 1, iy, iz - 1, MetalGridDim, INDEX_REAL)] += (float)(0.0473 * ExcessVal);*/
                  
                    val = 1; // It's a probability!
                  }
                  prob_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] += (float)val; // Added += after the last modification!
                  
                  if (prob_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] > 1)
                    prob_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = 1.0;  
                }
            break;
            
          case prop_count:
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  if (val < 0)
                    val = 0;
                  count_bubble_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = (int)val; 
                }
            break;
            
          case prop_Rave:
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  if (val < 0)
                    val = 0;
                  if (val > 0)
                    Rave_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = (float)val / count_bubble_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                }
            break;
            
          case prop_Rmax:
            for (int ix = 0; ix < slab_nix_metals[i_r]; ix++)
              for (int iy = 0; iy < MetalGridDim; iy++)
                for (int iz = 0; iz < MetalGridDim; iz++) {
                  double val = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
                  if (val >= Rmax_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)])
                    Rmax_grid_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] = (float)val; 
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
          (float)((grids->Probability_metals)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)]);
  write_grid_float("Probability_metals", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (float)((grids->R_ave)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)]);
  write_grid_float("Average Radius", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (float)((grids->R_max)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)]);
  write_grid_float("Max Radius", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (float)((grids->mass_IGM)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] * UnitMass_in_g / SOLAR_MASS);
  write_grid_float("mass_IGM", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
  for (int ii = 0; ii < local_nix_metals; ii++)
    for (int jj = 0; jj < MetalGridDim; jj++)
      for (int kk = 0; kk < MetalGridDim; kk++)
        grid[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)] =
          (int)((grids->N_bubbles)[grid_index(ii, jj, kk, MetalGridDim, INDEX_REAL)]);
  write_grid_float("N_bubbles", grid, file_id, fspace_id, memspace_id, dcpl_id);
  
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

  for (int ii = 0; ii < slab_n_real_metals; ii++) {
    grids->mass_metals[ii] = (float)0.;
    grids->mass_gas[ii] = (float)0.;
    grids->Zigm_box[ii] = (float)0.;
    grids->Probability_metals[ii] = (float)0.;
    grids->N_bubbles[ii] = (int)0;
    grids->R_ave[ii] = (float)0.;
    grids->R_max[ii] = (float)0.;
    grids->mass_IGM[ii] = (float)0.;
  }
}

void malloc_metal_grids()
{
  mlog("Allocating metal grids...", MLOG_OPEN);

  metal_grids_t* grids = &(run_globals.metal_grids);

  // run_globals.NStoreSnapshots is set in `initialize_halo_storage` Leave it commented because you might use that in the future when looking for scaling relations
  //run_globals.SnapshotDeltax = (float**)calloc((size_t)run_globals.NStoreSnapshots, sizeof(float*)); //?
  //run_globals.SnapshotVel = (float**)calloc((size_t)run_globals.NStoreSnapshots, sizeof(float*));  //?

  grids->galaxy_to_slab_map_metals = NULL;

  grids->mass_metals = NULL;
  grids->mass_gas = NULL;
  grids->Zigm_box = NULL; 
  grids->Probability_metals = NULL;
  grids->N_bubbles = NULL;
  grids->R_ave = NULL;
  grids->R_max = NULL;
  grids->mass_IGM = NULL;
  
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

  grids->mass_metals = fftwf_alloc_real((size_t)slab_n_real_metals);
  grids->mass_gas = fftwf_alloc_real((size_t)slab_n_real_metals);
  grids->Zigm_box = fftwf_alloc_real((size_t)slab_n_real_metals);

  grids->Probability_metals = fftwf_alloc_real((size_t)slab_n_real_metals);
  grids->N_bubbles = fftwf_alloc_real((size_t)slab_n_real_metals);
  grids->R_ave = fftwf_alloc_real((size_t)slab_n_real_metals);
  grids->R_max = fftwf_alloc_real((size_t)slab_n_real_metals);
  grids->mass_IGM = fftwf_alloc_real((size_t)slab_n_real_metals);
 
  init_metal_grids();

  mlog("...done", MLOG_CLOSE);
}

void free_metal_grids()
{
  mlog("Freeing metal grids...", MLOG_OPEN);

  metal_grids_t* grids = &(run_globals.metal_grids);

  free(run_globals.metal_grids.slab_ix_start_metals);
  free(run_globals.metal_grids.slab_nix_metals);

  fftwf_free(grids->mass_gas);
  fftwf_free(grids->mass_metals);
  fftwf_free(grids->Zigm_box);
  fftwf_free(grids->Probability_metals);
  fftwf_free(grids->N_bubbles);
  fftwf_free(grids->R_ave);
  fftwf_free(grids->R_max);
  fftwf_free(grids->mass_IGM);

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
    if (gal->Type < 3) {

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

void assign_probability_to_galaxies(int ngals_in_metal_slabs, int snapshot, int flag_property) //Right now you are assigning a probability to all galaxies! Actually you need only to the newly formed
{
  // Same way in which we assing Mcrit due to Reio and LW feedback in reionization.c

  gal_to_slab_t* galaxy_to_slab_map_metals = run_globals.metal_grids.galaxy_to_slab_map_metals;
  float* buffer_metals = run_globals.metal_grids.buffer_metals;
  ptrdiff_t* slab_nix_metals = run_globals.metal_grids.slab_nix_metals;
  ptrdiff_t* slab_ix_start_metals = run_globals.metal_grids.slab_ix_start_metals;
  int MetalGridDim = run_globals.params.MetalGridDim;
  double box_size = run_globals.params.BoxSize;
  float* Probability_metals = run_globals.metal_grids.Probability_metals;
  float* mass_metals = run_globals.metal_grids.mass_metals;
  float* mass_IGM = run_globals.metal_grids.mass_IGM;
  float* Rave_metals = run_globals.metal_grids.R_ave;
  float* Rmax_metals = run_globals.metal_grids.R_max;
  int total_assigned = 0;
  
  if (flag_property == 0) 
    mlog("Assigning probability for metals...", MLOG_OPEN);
    
  if (flag_property == 1)
    mlog("Assigning metals IGM for metals...", MLOG_OPEN);
    
  if (flag_property == 2)
    mlog("Assigning gas IGM for metals...", MLOG_OPEN);
    
  if (flag_property == 3)
    mlog("Assigning Rave for metals...", MLOG_OPEN);
    
  if (flag_property == 4)
    mlog("Assigning Rmax for metals...", MLOG_OPEN);

  int slab_map_offsets[run_globals.mpi_size];
  for (int ii = 0, i_gal = 0; ii < run_globals.mpi_size; ii++) {
    if (galaxy_to_slab_map_metals != NULL) {
      while ((i_gal < (ngals_in_metal_slabs - 1)) && (galaxy_to_slab_map_metals[i_gal].slab_ind < ii))
        i_gal++;

      if (galaxy_to_slab_map_metals[i_gal].slab_ind == ii)
        slab_map_offsets[ii] = i_gal;
      else
        slab_map_offsets[ii] = -1;
    } else
      slab_map_offsets[ii] = -1;
  }
  
    for (int i_skip = 0; i_skip < run_globals.mpi_size; i_skip++) {
      int recv_from_rank = (run_globals.mpi_rank + i_skip) % run_globals.mpi_size;
      int send_to_rank = (run_globals.mpi_rank - i_skip + run_globals.mpi_size) % run_globals.mpi_size;

      bool send_flag = false;
      bool recv_flag = (slab_map_offsets[recv_from_rank] > -1);
      
      if (flag_property == 0) {

        if (i_skip > 0) {
          MPI_Sendrecv(&recv_flag,
                       sizeof(bool),
                       MPI_BYTE,
                       recv_from_rank,
                       6393762,
                       &send_flag,
                       sizeof(bool),
                       MPI_BYTE,
                       send_to_rank,
                       6393762,
                       run_globals.mpi_comm,
                       MPI_STATUS_IGNORE);

          // need to ensure sends and receives do not clash!
            if (send_to_rank > run_globals.mpi_rank) {
              if (send_flag) {
                int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                MPI_Send(Probability_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
              }
              if (recv_flag) {
                int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
              }
            } 
            else {
              if (recv_flag) {
                int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
              } 
              if (send_flag) {
                int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                MPI_Send(Probability_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
              }
            }
          }
          else {
            int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
            memcpy(buffer_metals, Probability_metals, sizeof(float) * n_cells);
          }
        } 
        
      if (flag_property == 1) {
      
        if (i_skip > 0) {
            MPI_Sendrecv(&recv_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         recv_from_rank,
                         6393762,
                         &send_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         send_to_rank,
                         6393762,
                         run_globals.mpi_comm,
                         MPI_STATUS_IGNORE);

            // need to ensure sends and receives do not clash!
              if (send_to_rank > run_globals.mpi_rank) {
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(mass_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                }
              } 
              else {
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                } 
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(mass_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
              }
            }
            else {
              int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
              memcpy(buffer_metals, mass_metals, sizeof(float) * n_cells);
            }
          }
          
      if (flag_property == 2) {
      
        if (i_skip > 0) {
            MPI_Sendrecv(&recv_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         recv_from_rank,
                         6393762,
                         &send_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         send_to_rank,
                         6393762,
                         run_globals.mpi_comm,
                         MPI_STATUS_IGNORE);

            // need to ensure sends and receives do not clash!
              if (send_to_rank > run_globals.mpi_rank) {
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(mass_IGM, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                }
              } 
              else {
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                } 
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(mass_IGM, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
              }
            }
            else {
              int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
              memcpy(buffer_metals, mass_IGM, sizeof(float) * n_cells);
            }
          }    
          
      if (flag_property == 3) {
      
        if (i_skip > 0) {
            MPI_Sendrecv(&recv_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         recv_from_rank,
                         6393762,
                         &send_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         send_to_rank,
                         6393762,
                         run_globals.mpi_comm,
                         MPI_STATUS_IGNORE);

            // need to ensure sends and receives do not clash!
              if (send_to_rank > run_globals.mpi_rank) {
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(Rave_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                }
              } 
              else {
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                } 
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(Rave_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
              }
            }
            else {
              int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
              memcpy(buffer_metals, Rave_metals, sizeof(float) * n_cells);
            }
          }
          
      if (flag_property == 4) {
      
        if (i_skip > 0) {
            MPI_Sendrecv(&recv_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         recv_from_rank,
                         6393762,
                         &send_flag,
                         sizeof(bool),
                         MPI_BYTE,
                         send_to_rank,
                         6393762,
                         run_globals.mpi_comm,
                         MPI_STATUS_IGNORE);

            // need to ensure sends and receives do not clash!
              if (send_to_rank > run_globals.mpi_rank) {
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(Rmax_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                }
              } 
              else {
                if (recv_flag) {
                  int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
                  MPI_Recv(buffer_metals, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                } 
                if (send_flag) {
                  int n_cells = (int)(slab_nix_metals[run_globals.mpi_rank] * MetalGridDim * MetalGridDim);
                  MPI_Send(Rmax_metals, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
              }
            }
            else {
              int n_cells = (int)(slab_nix_metals[recv_from_rank] * MetalGridDim * MetalGridDim);
              memcpy(buffer_metals, Rmax_metals, sizeof(float) * n_cells);
            }
          }
      
      if (recv_flag) {
        int i_gal = slab_map_offsets[recv_from_rank];
        int ix_start_metals = (int)slab_ix_start_metals[recv_from_rank];
        while ((i_gal < ngals_in_metal_slabs) && (galaxy_to_slab_map_metals[i_gal].slab_ind == recv_from_rank)) {
          // TODO: We should use the position of the FOF group here...
          galaxy_t* gal = galaxy_to_slab_map_metals[i_gal].galaxy;
          int ix = pos_to_ngp(gal->Pos[0], box_size, MetalGridDim) - ix_start_metals;
          int iy = pos_to_ngp(gal->Pos[1], box_size, MetalGridDim);
          int iz = pos_to_ngp(gal->Pos[2], box_size, MetalGridDim);

          assert(ix >= 0);
          assert(ix < slab_nix_metals[recv_from_rank]);
          
          if (flag_property == 0)//{
            gal->Metal_Probability = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
            
          if (flag_property == 1)
              gal->Metals_IGM = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
              
          if (flag_property == 2) {
              //gal->Gas_IGM = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)] + cell_gas; 
              gal->Gas_IGM = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)]; 
    	      gal->Metallicity_IGM = calc_metallicity(gal->Gas_IGM, gal->Metals_IGM); //Maybe this can be put somewhere else
          }
          
          if (flag_property == 3)
            gal->AveBubble = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
            
          if (flag_property == 4)
            gal->MaxBubble = (double)buffer_metals[grid_index(ix, iy, iz, MetalGridDim, INDEX_REAL)];
           
          // increment counters
          i_gal++;
          total_assigned++;
        }
      }
    }

    if (total_assigned != ngals_in_metal_slabs)
      ABORT(EXIT_FAILURE);

  mlog("...done.", MLOG_CLOSE);
}


#endif
