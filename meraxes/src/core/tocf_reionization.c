#ifdef USE_TOCF

#include "meraxes.h"
#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>

void set_HII_eff_factor()
{
  // Use the params passed to Meraxes via the input file to set the HII ionising efficiency factor
  physics_params_t *params = &(run_globals.params.physics);

  // If we are using a redshift dependent escape fraction then reset
  // ReionEscapeFrac to one as we don't want to inlcude it in the
  // HII_eff_factor (it will be included in the stellar mass and SFR grids sent
  // to 21cmFAST instead).
  if (params->Flag_RedshiftDepEscFrac)
  {
    SID_log("Flag_RedshiftDepEscFrac is on => setting ReionEscapeFrac = 1.", SID_LOG_COMMENT);
    params->ReionEscapeFrac = 1.0;
  }

  // The following is based on Sobacchi & Messinger (2013) eqn 7
  // with f_* removed and f_b added since we define f_coll as M_*/M_tot rather than M_vir/M_tot,
  // and also with the inclusion of the effects of the Helium fraction.
  tocf_params.HII_eff_factor = 1.0 / run_globals.params.BaryonFrac
    * params->ReionNionPhotPerBary * params->ReionEscapeFrac / (1.0 - 0.75*tocf_params.Y_He);

  // Account for instantaneous recycling factor so that stellar mass is cumulative
  if (params->Flag_IRA)
    tocf_params.HII_eff_factor /= params->SfRecycleFraction;

  SID_log("Set value of tocf_params.HII_eff_factor = %g", SID_LOG_COMMENT, tocf_params.HII_eff_factor);
}



void assign_slabs()
{
  // Allocations made in this function are free'd in `free_reionization_grids`.

  // Assign the slab size
  int n_rank = SID.n_proc;
  int dim = tocf_params.HII_dim;

  // Use fftw to find out what slab each rank should get
  ptrdiff_t local_nix, local_ix_start;
  ptrdiff_t local_n_complex = fftwf_mpi_local_size_3d(dim, dim, dim/2 + 1, MPI_COMM_WORLD, &local_nix, &local_ix_start);

  // let every core know...
  ptrdiff_t *slab_nix = tocf_params.slab_nix;
  slab_nix = SID_malloc(sizeof(ptrdiff_t) * n_rank);  ///< array of number of x cells of every rank
  MPI_Allgather(&local_nix, sizeof(ptrdiff_t), MPI_BYTE, slab_nix, sizeof(ptrdiff_t), MPI_BYTE, MPI_COMM_WORLD);

  ptrdiff_t *slab_ix_start = tocf_params.slab_ix_start;
  slab_ix_start = SID_malloc(sizeof(ptrdiff_t) * n_rank); ///< array first x cell of every rank
  slab_ix_start[0] = 0;
  for(int ii=1; ii<n_rank; ii++)
    slab_ix_start[ii] = slab_ix_start[ii-1] + slab_nix[ii-1];

  ptrdiff_t *slab_n_complex = tocf_params.slab_n_complex;  ///< array of allocation counts for every rank
  slab_n_complex = SID_malloc(sizeof(ptrdiff_t) * n_rank);  ///< array of allocation counts for every rank
  MPI_Allgather(&local_n_complex, sizeof(ptrdiff_t), MPI_BYTE, slab_n_complex, sizeof(ptrdiff_t), MPI_BYTE, MPI_COMM_WORLD);
}



void call_find_HII_bubbles(int snapshot, int unsampled_snapshot, int nout_gals)
{
  // Thin wrapper round find_HII_bubbles

  int total_n_out_gals = 0;

  tocf_grids_t *grids = &(run_globals.tocf_grids);

  SID_log("Getting ready to call find_HII_bubbles...", SID_LOG_OPEN);

  // Check to see if there are actually any galaxies at this snapshot
  SID_Allreduce(&nout_gals, &total_n_out_gals, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  if (total_n_out_gals == 0)
  {
    SID_log("No galaxies in the simulation - skipping...", SID_LOG_CLOSE);
    return;
  }

  // Construct the stellar mass grid
  construct_stellar_grids(snapshot, nout_gals);

  SID_log("...done", SID_LOG_CLOSE);

  // Read in the dark matter density grid
  read_dm_grid(unsampled_snapshot, 0, (float*)(grids->deltax));

  // Call find_HII_bubbles
  SID_log("Calling find_HII_bubbles...", SID_LOG_OPEN | SID_LOG_TIMER);
  grids->global_xH = find_HII_bubbles(run_globals.ZZ[snapshot]);

  SID_log("grids->global_xH = %g", SID_LOG_COMMENT, grids->global_xH);
  SID_log("...done", SID_LOG_CLOSE);
}


void malloc_reionization_grids()
{

  tocf_grids_t *grids = &(run_globals.tocf_grids);

  grids->galaxy_to_slab_map = NULL;

  grids->xH                 = NULL;
  grids->stars              = NULL;
  grids->stars_unfiltered   = NULL;
  grids->stars_filtered     = NULL;
  grids->deltax             = NULL;
  grids->deltax_unfiltered  = NULL;
  grids->deltax_filtered    = NULL;
  grids->sfr                = NULL;
  grids->sfr_unfiltered     = NULL;
  grids->sfr_filtered       = NULL;
  grids->z_at_ionization    = NULL;
  grids->J_21_at_ionization = NULL;
  grids->J_21               = NULL;

  grids->global_xH = 1.0;
  grids->reion_complete = false;

  if (run_globals.params.TOCF_Flag)
  {

    assign_slabs();

    int HII_dim = tocf_params.HII_dim;
    ptrdiff_t *slab_nix = tocf_params.slab_nix;
    ptrdiff_t slab_n_real = slab_nix[SID.My_rank] * HII_dim * HII_dim;
    ptrdiff_t slab_n_complex = tocf_params.slab_n_complex[SID.My_rank];

    // create a buffer on each rank which is as large as the largest LOGICAL allocation on any single rank
    int max_cells = 0;

    for(int ii=0; ii < SID.n_proc; ii++)
      if(slab_nix[ii] > max_cells)
        max_cells = slab_nix[ii];

    max_cells *= HII_dim * HII_dim;
    grids->buffer_size = max_cells;

    grids->buffer          = fftwf_alloc_real(max_cells);
    grids->stars           = fftwf_alloc_real(slab_n_real * 2);  // padded for in-place FFT
    grids->stars_filtered  = fftwf_alloc_complex(slab_n_complex);
    grids->deltax          = fftwf_alloc_real(slab_n_real * 2);  // padded for in-place FFT
    grids->deltax_filtered = fftwf_alloc_complex(slab_n_complex);
    grids->sfr             = fftwf_alloc_real(slab_n_real * 2);  // padded for in-place FFT
    grids->sfr_filtered    = fftwf_alloc_complex(slab_n_complex);
    grids->xH              = fftwf_alloc_real(slab_n_real);
    grids->z_at_ionization = fftwf_alloc_real(slab_n_real);

    if (tocf_params.uvb_feedback)
    {
      grids->J_21_at_ionization = fftwf_alloc_real(slab_n_real);
      grids->J_21               = fftwf_alloc_real(slab_n_real);
      grids->Mvir_crit          = fftwf_alloc_real(slab_n_real);
    }

    SID_log("Initialising grids...", SID_LOG_COMMENT);

    for (int ii = 0; ii < slab_n_real; ii++)
    {
      grids->xH[ii] = 1.0;
      grids->z_at_ionization[ii] = -1;
    }
    memset(grids->xH, 1.0, sizeof(float) * slab_n_real);
    memset(grids->z_at_ionization, -1.0, sizeof(float) * slab_n_real);

    if (tocf_params.uvb_feedback)
    {
      memset(grids->J_21_at_ionization, 0., sizeof(float) * slab_n_real);
      memset(grids->J_21, 0., sizeof(float) * slab_n_real);
      memset(grids->Mvir_crit, 0, sizeof(float) * slab_n_real);
    }

    memset(grids->stars_filtered, 0, sizeof(fftwf_complex) * slab_n_complex);
    memset(grids->deltax_filtered, 0, sizeof(fftwf_complex) * slab_n_complex);
    memset(grids->sfr_filtered, 0, sizeof(fftwf_complex) * slab_n_complex);

    memset(grids->deltax, 0, sizeof(float) * slab_n_real * 2);
    memset(grids->stars, 0, sizeof(float) * slab_n_real * 2);
    memset(grids->sfr, 0, sizeof(float) * slab_n_real * 2);

    SID_log(" ...done", SID_LOG_CLOSE);
  }
}


void free_reionization_grids()
{
  SID_log("Freeing reionization grids...", SID_LOG_OPEN);

  tocf_grids_t *grids = &(run_globals.tocf_grids);

  SID_free(SID_FARG tocf_params.slab_n_complex);
  SID_free(SID_FARG tocf_params.slab_ix_start);
  SID_free(SID_FARG tocf_params.slab_nix);

  if (tocf_params.uvb_feedback)
  {
    fftwf_free(grids->J_21);
    fftwf_free(grids->J_21_at_ionization);
  }
  fftwf_free(grids->z_at_ionization);
  fftwf_free(grids->sfr_filtered);
  fftwf_free(grids->deltax_filtered);
  fftwf_free(grids->deltax);
  fftwf_free(grids->stars_filtered);
  fftwf_free(grids->xH);

  if (tocf_params.uvb_feedback)
    fftwf_free(grids->Mvir_crit);

  fftwf_free(grids->stars);
  fftwf_free(grids->sfr);
  fftwf_free(grids->buffer);

  SID_log(" ...done", SID_LOG_CLOSE);
}


int map_galaxies_to_slabs(int ngals)
{
  double box_size     = (double)(run_globals.params.BoxSize);
  int HII_dim         = tocf_params.HII_dim;

  // Loop through each valid galaxy and find what slab it sits in
  run_globals.tocf_grids.galaxy_to_slab_map = SID_malloc(sizeof(gal_to_slab_t) * ngals);
  gal_to_slab_t *galaxy_to_slab_map         = run_globals.tocf_grids.galaxy_to_slab_map;
  ptrdiff_t *slab_ix_start = tocf_params.slab_ix_start;

  galaxy_t *gal = run_globals.FirstGal;
  int gal_counter = 0;
  while (gal != NULL)
  {
    // TODO: Note that I am including ghosts here.  We will need to check the
    // validity of this.  By definition, if they are ghosts then their host
    // halo hasn't been identified at this time step and hence they haven't
    // been evolved.  Their properties (Sfr, StellarMass, etc.) will all have
    // been set when they were last identified.
    if (gal->Type < 3)
    {
      // TODO: for type 2 galaxies these positions will be set from the last
      // time they were identified.  If their host halo has moved significantly
      // since then, these positions won't reflect that and the satellites will
      // be spatially disconnected from their hosts.  We will need to fix this
      // at some point.

      // TODO: Get Greg to fix these positions to obey PBC!!
      if (gal->Pos[0] >= box_size)
        gal->Pos[0] -= box_size;
      else if (gal->Pos[0] < 0.0)
        gal->Pos[0] += box_size;
      ptrdiff_t ix = pos_to_cell(gal->Pos[0], box_size, HII_dim);

      assert((ix >= 0) && (ix < HII_dim));

      galaxy_to_slab_map[gal_counter].index = gal_counter;
      galaxy_to_slab_map[gal_counter].slab_ind = searchsorted(&ix, slab_ix_start, SID.n_proc, sizeof(ptrdiff_t), compare_ptrdiff, -1, -1);
      galaxy_to_slab_map[gal_counter++].galaxy = gal;
    }

    gal = gal->Next;
  }

  // sort the slab indices IN PLACE (n.b. compare_slab_assign is a stable comparison)
  qsort(galaxy_to_slab_map, gal_counter, sizeof(gal_to_slab_t), compare_slab_assign);

  return gal_counter;
}


void assign_Mvir_crit_to_galaxies(int ngals_in_slabs)
{

  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  gal_to_slab_t *galaxy_to_slab_map         = run_globals.tocf_grids.galaxy_to_slab_map;
  float *Mvir_crit = run_globals.tocf_grids.Mvir_crit;
  float *buffer = run_globals.tocf_grids.buffer;
  ptrdiff_t *slab_nix = tocf_params.slab_nix;
  ptrdiff_t *slab_ix_start = tocf_params.slab_ix_start;
  int HII_dim         = tocf_params.HII_dim;
  double box_size     = run_globals.params.BoxSize;

  int slab_map_offsets[SID.n_proc];
  for(int ii=0, i_gal=0; ii < SID.n_proc; ii++)
  {
    while((galaxy_to_slab_map[i_gal].slab_ind < ii) && (i_gal < ngals_in_slabs))
      i_gal++;

    if(galaxy_to_slab_map[i_gal].slab_ind == ii)
      slab_map_offsets[ii] = i_gal;
    else
      slab_map_offsets[ii] = -1;
  }

  for(int i_skip=0; i_skip < SID.n_proc; i_skip++)
  {
    int recv_from_rank = (SID.My_rank + i_skip) % SID.n_proc;
    int send_to_rank = (SID.My_rank - i_skip + SID.n_proc) % SID.n_proc;

    bool send_flag = false;
    bool recv_flag = (slab_map_offsets[recv_from_rank] > -1);
      
    if (i_skip > 0)
    {
      SID_Sendrecv(&recv_flag, sizeof(bool), SID_BYTE, recv_from_rank, 6393762,
          &send_flag, sizeof(bool), SID_BYTE, send_to_rank, 6393762, SID.COMM_WORLD);

      if(send_flag)
      {
        int n_cells = slab_nix[SID.My_rank]*HII_dim*HII_dim;
        SID_Send(Mvir_crit, n_cells, SID_FLOAT, send_to_rank, 793710, SID.COMM_WORLD);
      }

      if(recv_flag)
      {
        int n_cells = slab_nix[recv_from_rank]*HII_dim*HII_dim; 
        SID_Recv(buffer, n_cells, SID_FLOAT, recv_from_rank, 793710, SID.COMM_WORLD);
      }
    }
    else
    {
      int n_cells = slab_nix[recv_from_rank]*HII_dim*HII_dim; 
      memcpy(buffer, Mvir_crit, sizeof(float) * n_cells);
    }

    if(recv_flag)
    {
      int i_gal = slab_map_offsets[recv_from_rank];
      int ix_start = slab_ix_start[recv_from_rank];
      while((galaxy_to_slab_map[i_gal].slab_ind == recv_from_rank) && (i_gal < ngals_in_slabs))
      {
        // TODO: We should use the position of the FOF group here...
        galaxy_t *gal = galaxy_to_slab_map[i_gal].galaxy;
        int ix = pos_to_cell(gal->Pos[0], box_size, HII_dim) - ix_start;
        int iy = pos_to_cell(gal->Pos[1], box_size, HII_dim);
        int iz = pos_to_cell(gal->Pos[2], box_size, HII_dim);

        // Record the Mvir_crit (filtering mass) value
        gal->MvirCrit = (double)buffer[grid_index(ix, iy, iz, HII_dim, INDEX_REAL)];

        // increment counter
        i_gal++;
      }
    }

  }


}


void construct_stellar_grids(int snapshot, int ngals_in_slabs)
{
  double box_size     = (double)(run_globals.params.BoxSize);
  double Hubble_h     = run_globals.params.Hubble_h;
  float *stellar_grid = run_globals.tocf_grids.stars;
  float *sfr_grid     = run_globals.tocf_grids.sfr;
  int HII_dim         = tocf_params.HII_dim;
  run_units_t *units  = &(run_globals.units);
  double tHubble      = hubble_time(snapshot);

  gal_to_slab_t *galaxy_to_slab_map         = run_globals.tocf_grids.galaxy_to_slab_map;
  ptrdiff_t *slab_ix_start = tocf_params.slab_ix_start;
  int local_n_complex = (int)(tocf_params.slab_n_complex[SID.My_rank]);

  SID_log("Constructing stellar mass and sfr grids...", SID_LOG_OPEN | SID_LOG_TIMER);

  // init the grid
  for (int ii = 0; ii < local_n_complex; ii++)
  {
    *(stellar_grid + ii) = 0.0;
    *(sfr_grid + ii)     = 0.0;
  }

  // loop through each slab
  //
  // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
  // by slab index...
  int i_gal = 0;
  double cell_width = box_size / (double)HII_dim;
  ptrdiff_t *slab_nix = tocf_params.slab_nix;
  ptrdiff_t buffer_size = run_globals.tocf_grids.buffer_size;
  float *buffer = run_globals.tocf_grids.buffer;

  enum property { prop_stellar, prop_sfr };
  for(int prop = prop_stellar; prop <= prop_sfr; prop++)
  {
    for(int i_r=0; i_r < SID.n_proc; i_r++)
    {
      double min_xpos = (double)slab_ix_start[i_r] * cell_width;
      int nix = slab_nix[i_r];

      // init the buffer
      for(int ii=0; ii<buffer_size; ii++)
        buffer[ii] = 0.;

      // fill the local buffer for this slab
      while((i_gal < ngals_in_slabs) && (galaxy_to_slab_map[i_gal].slab_ind == i_r))
      {
        galaxy_t *gal = galaxy_to_slab_map[i_gal].galaxy;

        assert((galaxy_to_slab_map[i_gal].index >= 0) && (galaxy_to_slab_map[i_gal].index < ngals_in_slabs));
        assert((galaxy_to_slab_map[i_gal].slab_ind >= 0) && (galaxy_to_slab_map[i_gal].slab_ind < SID.n_proc));

        if (gal->Pos[0] >= box_size)
          gal->Pos[0] -= box_size;
        else if (gal->Pos[0] < 0.0)
          gal->Pos[0] += box_size;
        int ix = pos_to_cell(gal->Pos[0] - min_xpos, box_size, HII_dim);
        if (gal->Pos[1] >= box_size)
          gal->Pos[1] -= box_size;
        else if (gal->Pos[1] < 0.0)
          gal->Pos[1] += box_size;
        int iy = pos_to_cell(gal->Pos[1], box_size, HII_dim);
        if (gal->Pos[2] >= box_size)
          gal->Pos[2] -= box_size;
        else if (gal->Pos[2] < 0.0)
          gal->Pos[2] += box_size;
        int iz = pos_to_cell(gal->Pos[2], box_size, HII_dim);

        assert((ix < nix) && (ix >= 0));
        assert((iy < HII_dim) && (iy >= 0));
        assert((iz < HII_dim) && (iz >= 0));

        int ind = grid_index(ix, iy, iz, HII_dim, INDEX_REAL);

        assert((ind >=0) && (ind < nix*HII_dim*HII_dim));

        switch (prop) {
          case prop_stellar:
            buffer[ind] += gal->GrossStellarMass;
            break;

          case prop_sfr:
            buffer[ind] += gal->FescWeightedGSM;
            break;

          default:
            SID_log_error("Unrecognised property in slab creation.");
            ABORT(EXIT_FAILURE);
            break;
        }

        i_gal++;
      }

      // reduce on to the correct rank
      if(SID.My_rank == i_r)
        SID_Reduce(MPI_IN_PLACE, buffer, buffer_size, MPI_FLOAT, MPI_SUM, i_r, SID.COMM_WORLD);
      else
        SID_Reduce(buffer, buffer, buffer_size, MPI_FLOAT, MPI_SUM, i_r, SID.COMM_WORLD);

      if (SID.My_rank == i_r)
      {
        // copy the buffer into the real slab
        float *slab;
        switch (prop) {
          case prop_stellar:
            slab = stellar_grid;
            break;
          case prop_sfr:
            slab = sfr_grid;
            break;
        }

        int slab_size = slab_nix[i_r] * HII_dim * HII_dim;
        memcpy(slab, buffer, sizeof(float)*slab_size);

        // Do one final pass to put the grid in the correct (real) units (Msol or
        // Msol/s) and divide the sfr_grid by tHubble in order to convert the
        // stellar masses recorded into SFRs.
        for(int ii=0; ii < slab_size; ii++)
        {
          // put the values into the correct units
          switch (prop) {
            case prop_stellar:
              if (slab[ii] > 0)
                slab[ii] *= (1.e10 / Hubble_h);
              break;
            case prop_sfr:
              if (slab[ii] > 0)
                slab[ii] *= (units->UnitMass_in_g / units->UnitTime_in_s / SOLAR_MASS) / tHubble;
              break;
          }

          // check for underflow
          if (slab[ii] < 0)
            slab[ii] = 0;
        }
      }

    }
  }

  SID_log("done", SID_LOG_CLOSE);
}


void save_tocf_grids(hid_t parent_group_id, int snapshot)
{
  if (SID.My_rank == 0)
  {
    // Check if we even want to write anything...
    tocf_grids_t *grids = &(run_globals.tocf_grids);
    float epsilon = 0.0005;
    if ((grids->global_xH < epsilon) || (grids->global_xH > 1.0-epsilon))
      return;

    int   HII_dim       = tocf_params.HII_dim;
    hsize_t dims        = (hsize_t)pow(HII_dim, 3);
    float *grid;
    float *ps;
    int   ps_nbins;
    float average_deltaT;
    hid_t group_id;
    // double Hubble_h = run_globals.params.Hubble_h;
    
    // Save tocf grids
    // ----------------------------------------------------------------------------------------------------

    SID_log("Saving tocf grids...", SID_LOG_OPEN);

    group_id = H5Gcreate(parent_group_id, "Grids", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // float grids
    H5LTmake_dataset_float(group_id, "xH", 1, &dims, grids->xH);
    H5LTmake_dataset_float(group_id, "z_at_ionization", 1, &dims, grids->z_at_ionization);

    if (tocf_params.uvb_feedback)
    {
      H5LTmake_dataset_float(group_id, "J_21", 1, &dims, grids->J_21);
      H5LTmake_dataset_float(group_id, "J_21_at_ionization", 1, &dims, grids->J_21_at_ionization);
      H5LTmake_dataset_float(group_id, "Mvir_crit", 1, &dims, grids->Mvir_crit);
    }

    H5LTset_attribute_float(group_id, "xH", "global_xH", &(grids->global_xH), 1);

    // Save the escape fraction if we are using a redshift dependent escape fraction
    H5LTset_attribute_double(group_id, ".", "ReionEscapeFrac", &(run_globals.params.physics.ReionEscapeFrac), 1);

    // fftw padded grids
    grid = (float*)SID_calloc((int)dims * sizeof(float));
    for (int ii = 0; ii < HII_dim; ii++)
      for (int jj = 0; jj < HII_dim; jj++)
        for (int kk = 0; kk < HII_dim; kk++)
          grid[grid_index(ii, jj, kk, HII_dim, INDEX_REAL)] = ((float*)(grids->deltax))[grid_index(ii, jj, kk, HII_dim, INDEX_PADDED)];
    H5LTmake_dataset_float(group_id, "deltax", 1, &dims, grid);


    // Run delta_T_ps
    // ----------------------------------------------------------------------------------------------------

    SID_log("Calculating delta_T box and power spectrum...", SID_LOG_OPEN);

    memset((void*)grid, 0, sizeof(float) * (int)dims);

    delta_T_ps(
        run_globals.ZZ[snapshot],
        tocf_params.numcores,
        grids->xH,
        (float*)(grids->deltax),
        NULL,
        NULL,
        &average_deltaT,
        grid,
        &ps,
        &ps_nbins);

    H5LTmake_dataset_float(group_id, "delta_T", 1, &dims, grid);

    dims = ps_nbins * 3;
    H5LTmake_dataset_float(parent_group_id , "PowerSpectrum", 1               , &dims          , ps);
    H5LTset_attribute_int(parent_group_id  , "PowerSpectrum", "nbins"         , &ps_nbins      , 1);
    H5LTset_attribute_float(parent_group_id, "PowerSpectrum", "average_deltaT", &average_deltaT, 1);

    free(ps);

    SID_log("...done", SID_LOG_CLOSE);   // delta_T

    SID_free(SID_FARG grid);
    H5Gclose(group_id);

    SID_log("...done", SID_LOG_CLOSE);   // Saving tocf grids
  }
}


bool check_if_reionization_complete()
{
  bool complete = run_globals.tocf_grids.reion_complete;
  if (!complete)
  {
    if (SID.My_rank == 0)
    {
      complete = true;
      float *xH = run_globals.tocf_grids.xH;

      // If not all cells are ionised then reionization is still progressing...
      for (int ii=0; ii < HII_TOT_NUM_PIXELS; ii++)
      {
        if (xH[ii] != 0.0)
        {
          complete = false;
          break;
        }
      }
    }
    SID_Bcast(&complete, sizeof(bool), 0, SID.COMM_WORLD);
    run_globals.tocf_grids.reion_complete = complete;
  }
  return complete;
}
#endif
