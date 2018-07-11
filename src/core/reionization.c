#include "meraxes.h"
#include <assert.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <math.h>

void calculate_galaxy_fesc_vals(galaxy_t *gal, double new_stars, int snapshot)
{
    physics_params_t* params = &(run_globals.params.physics);

    float fesc_bh = params->EscapeFracBHNorm * (powf((1.0 + run_globals.ZZ[snapshot]) / 6.0, params->EscapeFracBHScaling));

    double fesc = params->EscapeFrac;

    if (params->EscapeFracDependency != 0) {
        double prop = 0.0;
        switch (params->EscapeFracDependency)
        {
            case 1:  // redshift
                prop = 1.0 + run_globals.ZZ[snapshot];
                break;
            case 2:  // stellar mass
                if (gal->StellarMass > 0.0)
                    prop = log10((gal->StellarMass + 1.0) / run_globals.params.Hubble_h) + 10.0;
                break;
            case 3: // star formation rate
                if (gal->Sfr > 0)
                    prop = log10((gal->Sfr + 1.0) * run_globals.units.UnitMass_in_g / run_globals.units.UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
            case 4: // cold gas density
                if (gal->ColdGas > 0.0)
                    prop = (log10((gal->ColdGas + 1.0)) + 10.) / (gal->DiskScaleLength * gal->DiskScaleLength * 100.) * run_globals.params.Hubble_h;
                break;
            case 5:  // halo mass
                if (gal->Mvir > 0.0)
                    prop = log10((gal->Mvir + 1.0) / run_globals.params.Hubble_h) + 10.0;
                break;
            default:
                mlog_error("Unrecognised EscapeFracDependency parameter value.");
        }

        double scaling = params->EscapeFracScaling / (params->EscapeFracP1*2) * 10.0;
        if (scaling > 10.)
            scaling = 10.;
        double x0 = pow(params->EscapeFracP0 - params->EscapeFracP1, scaling);
        double x1 = pow(params->EscapeFracP0 + params->EscapeFracP1, scaling);
        fesc = (pow(prop, scaling) - x0) / (x1 - x0);
    }
        
    if (fesc > 1.0)
        fesc = 1.0;
    else if (fesc < 0.0)
        fesc = 0.0;
    if (fesc_bh > 1.0)
        fesc_bh = 1.0;
    else if (fesc_bh < 0.0)
        fesc_bh = 0.0;

    gal->Fesc = fesc;
    gal->FescWeightedGSM += new_stars * fesc;

    // Here we just set the black hole escape fraction for use in the next
    // timestep when calculating previous_merger_driven_BH_growth().  It is in
    // this function that the EffectiveBHM will be updated.  For the galaxy
    // that we are working on now, the EffectiveBHM has already been calculated
    // using the FescBH value from the previous snapshot.
    // The upshot is that we don't need to do anything to the EffectiveBHM
    // here.  It's confusing I know.  I intend to re-write this to make things
    // more obvious at some point in the future.
    gal->FescBH = fesc_bh;
}

void set_quasar_fobs()
{
    physics_params_t* params = &(run_globals.params.physics);

    params->quasar_fobs = 1. - cos(params->quasar_open_angle / 180. * M_PI / 2.);
    mlog("Quasar radiation open angle is set to be %g, corresponding to an obscure fraction of %g",
        MLOG_MESG | MLOG_FLUSH, params->quasar_open_angle, params->quasar_fobs);
}

void set_ReionEfficiency()
{
    // Use the params passed to Meraxes via the input file to set the HII ionising efficiency factor
    physics_params_t* params = &(run_globals.params.physics);


    // The following is based on Sobacchi & Messinger (2013) eqn 7
    // with f_* removed and f_b added since we define f_coll as M_*/M_tot rather than M_vir/M_tot,
    // and also with the inclusion of the effects of the Helium fraction.
    params->ReionEfficiency = 1.0 / run_globals.params.BaryonFrac
        * params->ReionNionPhotPerBary / (1.0 - 0.75 * params->Y_He);

    // Account for instantaneous recycling factor so that stellar mass is cumulative
    if (params->Flag_IRA)
        params->ReionEfficiency /= params->SfRecycleFraction;

    mlog("Set value of run_globals.params.ReionEfficiency = %g", MLOG_MESG, params->ReionEfficiency);
}

void assign_slabs()
{
    mlog("Assigning slabs to MPI cores...", MLOG_OPEN);

    // Allocations made in this function are free'd in `free_reionization_grids`.
    fftwf_mpi_init();

    // Assign the slab size
    int n_rank = run_globals.mpi_size;
    int dim = run_globals.params.ReionGridDim;

    // Use fftw to find out what slab each rank should get
    ptrdiff_t local_nix, local_ix_start;
    ptrdiff_t local_n_complex = fftwf_mpi_local_size_3d(dim, dim, dim / 2 + 1, run_globals.mpi_comm, &local_nix, &local_ix_start);

    // let every core know...
    ptrdiff_t** slab_nix = &run_globals.reion_grids.slab_nix;
    *slab_nix = malloc(sizeof(ptrdiff_t) * n_rank); ///< array of number of x cells of every rank
    MPI_Allgather(&local_nix, sizeof(ptrdiff_t), MPI_BYTE, *slab_nix, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm);

    ptrdiff_t** slab_ix_start = &run_globals.reion_grids.slab_ix_start;
    *slab_ix_start = malloc(sizeof(ptrdiff_t) * n_rank); ///< array first x cell of every rank
    (*slab_ix_start)[0] = 0;
    for (int ii = 1; ii < n_rank; ii++)
        (*slab_ix_start)[ii] = (*slab_ix_start)[ii - 1] + (*slab_nix)[ii - 1];

    ptrdiff_t** slab_n_complex = &run_globals.reion_grids.slab_n_complex; ///< array of allocation counts for every rank
    *slab_n_complex = malloc(sizeof(ptrdiff_t) * n_rank); ///< array of allocation counts for every rank
    MPI_Allgather(&local_n_complex, sizeof(ptrdiff_t), MPI_BYTE, *slab_n_complex, sizeof(ptrdiff_t), MPI_BYTE, run_globals.mpi_comm);

    mlog("...done", MLOG_CLOSE);
}

void call_find_HII_bubbles(int snapshot, int nout_gals, timer_info* timer)
{
    // Thin wrapper round find_HII_bubbles

    int total_n_out_gals = 0;

    reion_grids_t* grids = &(run_globals.reion_grids);

    mlog("Getting ready to call find_HII_bubbles...", MLOG_OPEN);

    // Check to see if there are actually any galaxies at this snapshot
    MPI_Allreduce(&nout_gals, &total_n_out_gals, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    if (total_n_out_gals == 0) {
        mlog("No galaxies in the simulation - skipping...", MLOG_CLOSE);
        return;
    }

    // Construct the baryon grids
    construct_baryon_grids(snapshot, nout_gals);

    // Read in the dark matter density grid
    switch (run_globals.params.TreesID) {
        case VELOCIRAPTOR_TREES:
            read_dm_grid__velociraptor(snapshot, grids->deltax);
            break;
        case GBPTREES_TREES:
            read_dm_grid__gbptrees(snapshot, grids->deltax);
            break;
        default:
            mlog_error("Unrecognised input trees identifier (TreesID).");
            break;
    }

    // save the grids prior to doing FFTs to avoid precision loss and aliasing etc.
    for (int i_out = 0; i_out < run_globals.NOutputSnaps; i_out++)
        if (snapshot == run_globals.ListOutputSnaps[i_out] && run_globals.params.Flag_OutputGrids)
            save_reion_input_grids(snapshot);

    mlog("...done", MLOG_CLOSE);

    // Call find_HII_bubbles
    mlog("Calling find_HII_bubbles", MLOG_OPEN | MLOG_TIMERSTART);

    // Call find_HII_bubbles
    find_HII_bubbles(snapshot, timer);

    mlog("grids->volume_weighted_global_xH = %g", MLOG_MESG, grids->volume_weighted_global_xH);
    mlog("global mass weighted xHII = %g at z = %g", MLOG_MESG, 1.0 - grids->mass_weighted_global_xH, run_globals.ZZ[snapshot]);
    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

void init_reion_grids()
{
    reion_grids_t* grids = &(run_globals.reion_grids);
    int ReionGridDim = run_globals.params.ReionGridDim;
    ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
    ptrdiff_t slab_n_real = slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim; // TODO: NOT WORKING!!!
    ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

    mlog("Initialising grids...", MLOG_MESG);

    grids->volume_weighted_global_xH = 1.0;
    grids->mass_weighted_global_xH = 1.0;
    grids->started = 0;
    grids->finished = 0;

    for (int ii = 0; ii < slab_n_real; ii++) {
        grids->xH[ii] = 1.0;
        grids->z_at_ionization[ii] = -1;
        grids->r_bubble[ii] = 0.0;
    }

    for (int ii = 0; ii < slab_n_real; ii++)
        if (run_globals.params.ReionUVBFlag) {
            grids->J_21_at_ionization[ii] = 0.;
            grids->J_21[ii] = 0.;
            grids->Mvir_crit[ii] = 0;
        }

    for (int ii = 0; ii < slab_n_complex; ii++) {
        grids->stars_filtered[ii] = 0 + 0 * I;
        grids->deltax_filtered[ii] = 0 + 0 * I;
        grids->sfr_filtered[ii] = 0 + 0 * I;
    }

    for (int ii = 0; ii < slab_n_complex * 2; ii++) {
        grids->deltax[ii] = 0;
        grids->stars[ii] = 0;
        grids->sfr[ii] = 0;
    }

    mlog(" ...done", MLOG_CLOSE);
}

void malloc_reionization_grids()
{
    reion_grids_t* grids = &(run_globals.reion_grids);

    // run_globals.NStoreSnapshots is set in `initialize_halo_storage`
    run_globals.SnapshotDeltax = (float**)calloc(run_globals.NStoreSnapshots, sizeof(float*));

    grids->galaxy_to_slab_map = NULL;

    grids->xH = NULL;
    grids->stars = NULL;
    grids->stars_unfiltered = NULL;
    grids->stars_filtered = NULL;
    grids->deltax = NULL;
    grids->deltax_unfiltered = NULL;
    grids->deltax_filtered = NULL;
    grids->sfr = NULL;
    grids->sfr_unfiltered = NULL;
    grids->sfr_filtered = NULL;
    grids->z_at_ionization = NULL;
    grids->J_21_at_ionization = NULL;
    grids->J_21 = NULL;

    if (run_globals.params.Flag_PatchyReion) {
        assign_slabs();

        int ReionGridDim = run_globals.params.ReionGridDim;
        ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
        ptrdiff_t slab_n_real = slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim; // TODO: NOT WORKING!!!
        ptrdiff_t slab_n_complex = run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank];

        // create a buffer on each rank which is as large as the largest LOGICAL allocation on any single rank
        int max_cells = 0;

        for (int ii = 0; ii < run_globals.mpi_size; ii++)
            if (slab_nix[ii] > max_cells)
                max_cells = slab_nix[ii];

        max_cells *= ReionGridDim * ReionGridDim;
        grids->buffer_size = max_cells;

        grids->buffer = fftwf_alloc_real(max_cells);
        grids->stars = fftwf_alloc_real(slab_n_complex * 2); // padded for in-place FFT
        grids->stars_filtered = fftwf_alloc_complex(slab_n_complex);
        grids->deltax = fftwf_alloc_real(slab_n_complex * 2); // padded for in-place FFT
        grids->deltax_filtered = fftwf_alloc_complex(slab_n_complex);
        grids->sfr = fftwf_alloc_real(slab_n_complex * 2); // padded for in-place FFT
        grids->sfr_filtered = fftwf_alloc_complex(slab_n_complex);
        grids->xH = fftwf_alloc_real(slab_n_real);
        grids->z_at_ionization = fftwf_alloc_real(slab_n_real);
        grids->r_bubble = fftwf_alloc_real(slab_n_real);

        if (run_globals.params.ReionUVBFlag) {
            grids->J_21_at_ionization = fftwf_alloc_real(slab_n_real);
            grids->J_21 = fftwf_alloc_real(slab_n_real);
            grids->Mvir_crit = fftwf_alloc_real(slab_n_real);
        }

        init_reion_grids();
    }
}

void free_reionization_grids()
{
    mlog("Freeing reionization grids...", MLOG_OPEN);

    reion_grids_t* grids = &(run_globals.reion_grids);

    free(run_globals.reion_grids.slab_n_complex);
    free(run_globals.reion_grids.slab_ix_start);
    free(run_globals.reion_grids.slab_nix);

    if (run_globals.params.ReionUVBFlag) {
        fftwf_free(grids->J_21);
        fftwf_free(grids->J_21_at_ionization);
    }
    fftwf_free(grids->z_at_ionization);
    fftwf_free(grids->sfr_filtered);
    fftwf_free(grids->deltax_filtered);
    fftwf_free(grids->deltax);
    fftwf_free(grids->stars_filtered);
    fftwf_free(grids->xH);

    if (run_globals.params.ReionUVBFlag)
        fftwf_free(grids->Mvir_crit);

    fftwf_free(grids->stars);
    fftwf_free(grids->sfr);
    fftwf_free(grids->buffer);

    mlog(" ...done", MLOG_CLOSE);
}

int map_galaxies_to_slabs(int ngals)
{
    double box_size = run_globals.params.BoxSize;
    int ReionGridDim = run_globals.params.ReionGridDim;

    mlog("Mapping galaxies to slabs...", MLOG_OPEN);

    // Loop through each valid galaxy and find what slab it sits in
    if (ngals > 0)
        run_globals.reion_grids.galaxy_to_slab_map = malloc(sizeof(gal_to_slab_t) * ngals);
    else
        run_globals.reion_grids.galaxy_to_slab_map = NULL;

    gal_to_slab_t* galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
    ptrdiff_t* slab_ix_start = run_globals.reion_grids.slab_ix_start;

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

            ptrdiff_t ix = pos_to_ngp(gal->Pos[0], box_size, ReionGridDim);

            assert((ix >= 0) && (ix < ReionGridDim));

            galaxy_to_slab_map[gal_counter].index = gal_counter;
            galaxy_to_slab_map[gal_counter].slab_ind = searchsorted(&ix, slab_ix_start, run_globals.mpi_size, sizeof(ptrdiff_t), compare_ptrdiff, -1, -1);
            galaxy_to_slab_map[gal_counter++].galaxy = gal;
        }

        gal = gal->Next;
    }

    // sort the slab indices IN PLACE (n.b. compare_slab_assign is a stable comparison)
    if (galaxy_to_slab_map != NULL)
        qsort(galaxy_to_slab_map, (size_t)gal_counter, sizeof(gal_to_slab_t), compare_slab_assign);

    assert(gal_counter == ngals);

    mlog("...done.", MLOG_CLOSE);

    return gal_counter;
}

void assign_Mvir_crit_to_galaxies(int ngals_in_slabs)
{
    // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
    // by slab index...
    gal_to_slab_t* galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
    float* Mvir_crit = run_globals.reion_grids.Mvir_crit;
    float* buffer = run_globals.reion_grids.buffer;
    ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
    ptrdiff_t* slab_ix_start = run_globals.reion_grids.slab_ix_start;
    int ReionGridDim = run_globals.params.ReionGridDim;
    double box_size = run_globals.params.BoxSize;
    int total_assigned = 0;

    mlog("Assigning Mvir_crit to galaxies...", MLOG_OPEN);

    // Work out the index of the galaxy_to_slab_map where each slab begins.
    int slab_map_offsets[run_globals.mpi_size];
    for (int ii = 0, i_gal = 0; ii < run_globals.mpi_size; ii++) {
        if (galaxy_to_slab_map != NULL) {
            while ((galaxy_to_slab_map[i_gal].slab_ind < ii) && (i_gal < ngals_in_slabs))
                i_gal++;

            if (galaxy_to_slab_map[i_gal].slab_ind == ii)
                slab_map_offsets[ii] = i_gal;
            else
                slab_map_offsets[ii] = -1;
        }
        else
            // if this core has no galaxies then the offsets are -1 everywhere
            slab_map_offsets[ii] = -1;
    }

    // do a ring exchange of slabs between all cores
    for (int i_skip = 0; i_skip < run_globals.mpi_size; i_skip++) {
        int recv_from_rank = (run_globals.mpi_rank + i_skip) % run_globals.mpi_size;
        int send_to_rank = (run_globals.mpi_rank - i_skip + run_globals.mpi_size) % run_globals.mpi_size;

        bool send_flag = false;
        bool recv_flag = (slab_map_offsets[recv_from_rank] > -1);

        if (i_skip > 0) {
            MPI_Sendrecv(&recv_flag, sizeof(bool), MPI_BYTE, recv_from_rank, 6393762,
                &send_flag, sizeof(bool), MPI_BYTE, send_to_rank, 6393762, run_globals.mpi_comm, MPI_STATUS_IGNORE);

            // need to ensure sends and receives do not clash!
            if (send_to_rank > run_globals.mpi_rank) {
                if (send_flag) {
                    int n_cells = slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim;
                    MPI_Send(Mvir_crit, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
                if (recv_flag) {
                    int n_cells = slab_nix[recv_from_rank] * ReionGridDim * ReionGridDim;
                    MPI_Recv(buffer, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                }
            }
            else {
                if (recv_flag) {
                    int n_cells = slab_nix[recv_from_rank] * ReionGridDim * ReionGridDim;
                    MPI_Recv(buffer, n_cells, MPI_FLOAT, recv_from_rank, 793710, run_globals.mpi_comm, MPI_STATUS_IGNORE);
                }
                if (send_flag) {
                    int n_cells = slab_nix[run_globals.mpi_rank] * ReionGridDim * ReionGridDim;
                    MPI_Send(Mvir_crit, n_cells, MPI_FLOAT, send_to_rank, 793710, run_globals.mpi_comm);
                }
            }
        }
        else {
            int n_cells = slab_nix[recv_from_rank] * ReionGridDim * ReionGridDim;
            memcpy(buffer, Mvir_crit, sizeof(float) * n_cells);
        }

        // if this core has received a slab of Mvir_crit then assign values to the
        // galaxies which belong to this slab
        if (recv_flag) {
            int i_gal = slab_map_offsets[recv_from_rank];
            int ix_start = slab_ix_start[recv_from_rank];
            while ((galaxy_to_slab_map[i_gal].slab_ind == recv_from_rank) && (i_gal < ngals_in_slabs)) {
                // TODO: We should use the position of the FOF group here...
                galaxy_t* gal = galaxy_to_slab_map[i_gal].galaxy;
                int ix = pos_to_ngp(gal->Pos[0], box_size, ReionGridDim) - ix_start;
                int iy = pos_to_ngp(gal->Pos[1], box_size, ReionGridDim);
                int iz = pos_to_ngp(gal->Pos[2], box_size, ReionGridDim);

                assert(ix >= 0);
                assert(ix < slab_nix[recv_from_rank]);

                // Record the Mvir_crit (filtering mass) value
                gal->MvirCrit = (double)buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];

                // increment counters
                i_gal++;
                total_assigned++;
            }
        }
    }

    if (total_assigned != ngals_in_slabs)
        ABORT(EXIT_FAILURE);

    mlog("...done.", MLOG_CLOSE);
}

void construct_baryon_grids(int snapshot, int local_ngals)
{
    double box_size = run_globals.params.BoxSize;
    float* stellar_grid = run_globals.reion_grids.stars;
    float* sfr_grid = run_globals.reion_grids.sfr;
    int ReionGridDim = run_globals.params.ReionGridDim;
    double tHubble = hubble_time(snapshot);

    gal_to_slab_t* galaxy_to_slab_map = run_globals.reion_grids.galaxy_to_slab_map;
    ptrdiff_t* slab_ix_start = run_globals.reion_grids.slab_ix_start;
    int local_n_complex = (int)(run_globals.reion_grids.slab_n_complex[run_globals.mpi_rank]);

    mlog("Constructing stellar mass and sfr grids...", MLOG_OPEN | MLOG_TIMERSTART);

    // init the grid
    for (int ii = 0; ii < local_n_complex * 2; ii++) {
        stellar_grid[ii] = 0.0;
        sfr_grid[ii] = 0.0;
    }

    // loop through each slab
    //
    // N.B. We are assuming here that the galaxy_to_slab mapping has been sorted
    // by slab index...
    ptrdiff_t* slab_nix = run_globals.reion_grids.slab_nix;
    ptrdiff_t buffer_size = run_globals.reion_grids.buffer_size;
    float* buffer = run_globals.reion_grids.buffer;

    enum property { prop_stellar,
        prop_sfr };
    for (int prop = prop_stellar; prop <= prop_sfr; prop++) {
        int i_gal = 0;
        int skipped_gals = 0;
        long N_BlackHoleMassLimitReion = 0;

        for (int i_r = 0; i_r < run_globals.mpi_size; i_r++) {
            // init the buffer
            for (int ii = 0; ii < buffer_size; ii++)
                buffer[ii] = 0.;

            // if this core holds no galaxies then we don't need to fill the buffer
            if (local_ngals != 0)
                // fill the local buffer for this slab
                while (((i_gal - skipped_gals) < local_ngals) && (galaxy_to_slab_map[i_gal].slab_ind == i_r)) {
                    galaxy_t* gal = galaxy_to_slab_map[i_gal].galaxy;

                    // Dead galaxies should not be included here and are not in the
                    // local_ngals count.  They will, however, have been assigned to a
                    // slab so we will need to ignore them here...
                    if (gal->Type > 2) {
                        i_gal++;
                        skipped_gals++;
                        continue;
                    }

                    assert(galaxy_to_slab_map[i_gal].index >= 0);
                    assert((galaxy_to_slab_map[i_gal].slab_ind >= 0) && (galaxy_to_slab_map[i_gal].slab_ind < run_globals.mpi_size));

                    int ix = pos_to_ngp(gal->Pos[0], box_size, ReionGridDim) - slab_ix_start[i_r];
                    int iy = pos_to_ngp(gal->Pos[1], box_size, ReionGridDim);
                    int iz = pos_to_ngp(gal->Pos[2], box_size, ReionGridDim);

                    assert((ix < slab_nix[i_r]) && (ix >= 0));
                    assert((iy < ReionGridDim) && (iy >= 0));
                    assert((iz < ReionGridDim) && (iz >= 0));

                    int ind = grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL);

                    assert((ind >= 0) && (ind < slab_nix[i_r] * ReionGridDim * ReionGridDim));

                    // They are the same just now, but may be different in the future once the model is improved.
                    switch (prop) {
                    case prop_stellar:
                        buffer[ind] += gal->FescWeightedGSM;
                        // a trick to include quasar radiation using current 21cmFAST code
                        if (run_globals.params.physics.Flag_BHFeedback) {
                            if (gal->BlackHoleMass >= run_globals.params.physics.BlackHoleMassLimitReion)
                                buffer[ind] += gal->EffectiveBHM;
                            else
                                N_BlackHoleMassLimitReion += 1;
                        }
                        break;

                    case prop_sfr:
                        buffer[ind] += gal->FescWeightedGSM;
                        // for ionizing_source_formation_rate_grid, need further convertion due to different UV spectral index of quasar and stellar component
                        if (run_globals.params.physics.Flag_BHFeedback)
                            if (gal->BlackHoleMass >= run_globals.params.physics.BlackHoleMassLimitReion)
                                buffer[ind] += gal->EffectiveBHM * run_globals.params.physics.ReionAlphaUVBH / run_globals.params.physics.ReionAlphaUV;
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
                MPI_Reduce(MPI_IN_PLACE, buffer, buffer_size, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);
            else
                MPI_Reduce(buffer, buffer, buffer_size, MPI_FLOAT, MPI_SUM, i_r, run_globals.mpi_comm);

            if (run_globals.mpi_rank == i_r)

                // Do one final pass and divide the sfr_grid by tHubble
                // in order to convert the stellar masses recorded into SFRs before
                // finally copying the values into the appropriate slab.
                // TODO: Use a better timescale for SFR
                switch (prop) {
                case prop_sfr:
                    for (int ix = 0; ix < slab_nix[i_r]; ix++)
                        for (int iy = 0; iy < ReionGridDim; iy++)
                            for (int iz = 0; iz < ReionGridDim; iz++) {
                                double val = (double)buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];
                                val = (val > 0) ? val / tHubble : 0;
                                sfr_grid[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = (float)val;
                            }
                    break;

                case prop_stellar:
                    for (int ix = 0; ix < slab_nix[i_r]; ix++)
                        for (int iy = 0; iy < ReionGridDim; iy++)
                            for (int iz = 0; iz < ReionGridDim; iz++) {
                                float val = buffer[grid_index(ix, iy, iz, ReionGridDim, INDEX_REAL)];
                                if (val < 0)
                                    val = 0;
                                stellar_grid[grid_index(ix, iy, iz, ReionGridDim, INDEX_PADDED)] = val;
                            }
                    break;

                default:
                    mlog_error("Eh!?!");
                    ABORT(EXIT_FAILURE);
                }
        }
        MPI_Allreduce(MPI_IN_PLACE, &N_BlackHoleMassLimitReion, 1, MPI_DOUBLE, MPI_SUM, run_globals.mpi_comm);
        mlog("%d quasars are smaller than %g", MLOG_MESG, N_BlackHoleMassLimitReion, run_globals.params.physics.BlackHoleMassLimitReion);
    }

    mlog("done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

static void write_grid_float(const char* name, float* data, hid_t file_id, hid_t fspace_id, hid_t memspace_id, hid_t dcpl_id)
{
    // create the dataset
    hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, fspace_id,
        H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    // create the property list
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // write the dataset
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id, fspace_id, plist_id, data);

    // cleanup
    H5Pclose(plist_id);
    H5Dclose(dset_id);
}

void gen_grids_fname(int snapshot, char* name, bool relative)
{
    if (!relative)
        sprintf(name, "%s/%s_grids_%d.hdf5", run_globals.params.OutputDir, run_globals.params.FileNameGalaxies, snapshot);
    else
        sprintf(name, "%s_grids_%d.hdf5", run_globals.params.FileNameGalaxies, snapshot);
}

void save_reion_input_grids(int snapshot)
{
    reion_grids_t* grids = &(run_globals.reion_grids);
    int ReionGridDim = run_globals.params.ReionGridDim;
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);
    double UnitTime_in_s = run_globals.units.UnitTime_in_s;
    double UnitMass_in_g = run_globals.units.UnitMass_in_g;

    mlog("Saving tocf input grids...", MLOG_OPEN);

    char name[STRLEN];
    gen_grids_fname(snapshot, name, false);

    // create the file (in parallel)
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    // create the filespace
    hsize_t dims[3] = { ReionGridDim, ReionGridDim, ReionGridDim };
    hid_t fspace_id = H5Screate_simple(3, dims, NULL);

    // create the memspace
    hsize_t mem_dims[3] = { local_nix, ReionGridDim, ReionGridDim };
    hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

    // select a hyperslab in the filespace
    hsize_t start[3] = { run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank], 0, 0 };
    hsize_t count[3] = { local_nix, ReionGridDim, ReionGridDim };
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

    // set the dataset creation property list to use chunking along x-axis
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, (hsize_t[3]){ 1, ReionGridDim, ReionGridDim });
    
    // fftw padded grids
    float* grid = (float*)calloc(local_nix * ReionGridDim * ReionGridDim, sizeof(float));

    for (int ii = 0; ii < local_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++)
                grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] = (grids->deltax)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];
    write_grid_float("deltax", grid, file_id, fspace_id, memspace_id, dcpl_id);

    for (int ii = 0; ii < local_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++)
                grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] = (grids->stars)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)];
    write_grid_float("stars", grid, file_id, fspace_id, memspace_id, dcpl_id);

    for (int ii = 0; ii < local_nix; ii++)
        for (int jj = 0; jj < ReionGridDim; jj++)
            for (int kk = 0; kk < ReionGridDim; kk++)
                grid[grid_index(ii, jj, kk, ReionGridDim, INDEX_REAL)] = (grids->sfr)[grid_index(ii, jj, kk, ReionGridDim, INDEX_PADDED)]
                    * UnitMass_in_g / UnitTime_in_s
                    * SEC_PER_YEAR / SOLAR_MASS;
    write_grid_float("sfr", grid, file_id, fspace_id, memspace_id, dcpl_id);

    // tidy up
    free(grid);
    H5Pclose(dcpl_id);
    H5Sclose(memspace_id);
    H5Sclose(fspace_id);
    H5Fclose(file_id);

    mlog("...done", MLOG_CLOSE);
}

void save_reion_output_grids(int snapshot)
{
    reion_grids_t* grids = &(run_globals.reion_grids);
    int ReionGridDim = run_globals.params.ReionGridDim;
    int local_nix = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]);

    // float *ps;
    // int   ps_nbins;
    // float average_deltaT;
    // double Hubble_h = run_globals.params.Hubble_h;

    // Save tocf grids
    // ----------------------------------------------------------------------------------------------------

    mlog("Saving tocf output grids...", MLOG_OPEN);

    char name[STRLEN];
    gen_grids_fname(snapshot, name, false);

    // open the file (in parallel)
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);
    hid_t file_id = H5Fopen(name, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // create the filespace
    hsize_t dims[3] = { ReionGridDim, ReionGridDim, ReionGridDim };
    hid_t fspace_id = H5Screate_simple(3, dims, NULL);

    // create the memspace
    hsize_t mem_dims[3] = { local_nix, ReionGridDim, ReionGridDim };
    hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

    // select a hyperslab in the filespace
    hsize_t start[3] = { run_globals.reion_grids.slab_ix_start[run_globals.mpi_rank], 0, 0 };
    hsize_t count[3] = { local_nix, ReionGridDim, ReionGridDim };
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

    // set the dataset creation property list to use chunking along x-axis
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, (hsize_t[3]){ 1, ReionGridDim, ReionGridDim });

    // create and write the datasets
    write_grid_float("xH", grids->xH, file_id, fspace_id, memspace_id, dcpl_id);
    write_grid_float("z_at_ionization", grids->z_at_ionization, file_id, fspace_id, memspace_id, dcpl_id);
    write_grid_float("r_bubble", grids->r_bubble, file_id, fspace_id, memspace_id, dcpl_id);

    if (run_globals.params.ReionUVBFlag) {
        write_grid_float("J_21", grids->J_21, file_id, fspace_id, memspace_id, dcpl_id);
        write_grid_float("J_21_at_ionization", grids->J_21_at_ionization, file_id, fspace_id, memspace_id, dcpl_id);
        write_grid_float("Mvir_crit", grids->Mvir_crit, file_id, fspace_id, memspace_id, dcpl_id);
    }

    H5LTset_attribute_double(file_id, "xH", "volume_weighted_global_xH", &(grids->volume_weighted_global_xH), 1);
    H5LTset_attribute_double(file_id, "xH", "mass_weighted_global_xH", &(grids->mass_weighted_global_xH), 1);

    // // Run delta_T_ps
    // // ----------------------------------------------------------------------------------------------------

    // mlog("Calculating delta_T box and power spectrum...", MLOG_OPEN);

    // memset((void*)grid, 0, sizeof(float) * (int)dims);

    // delta_T_ps(
    //     run_globals.ZZ[snapshot],
    //     run_globals.params.numcores,
    //     grids->xH,
    //     (float*)(grids->deltax),
    //     &average_deltaT,
    //     grid,
    //     &ps,
    //     &ps_nbins);

    // H5LTmake_dataset_float(group_id, "delta_T", 1, &dims, grid);

    // dims = ps_nbins * 3;
    // H5LTmake_dataset_float(parent_group_id , "PowerSpectrum", 1               , &dims          , ps);
    // H5LTset_attribute_int(parent_group_id  , "PowerSpectrum", "nbins"         , &ps_nbins      , 1);
    // H5LTset_attribute_float(parent_group_id, "PowerSpectrum", "average_deltaT", &average_deltaT, 1);

    // free(ps);

    // mlog("...done", MLOG_CLOSE);   // delta_T

    // tidy up
    H5Pclose(dcpl_id);
    H5Sclose(memspace_id);
    H5Sclose(fspace_id);
    H5Fclose(file_id);

    mlog("...done", MLOG_CLOSE); // Saving tocf grids
}

bool check_if_reionization_ongoing()
{
    int started = run_globals.reion_grids.started;
    int finished = run_globals.reion_grids.finished;

    // First check if we've already finished on all cores.
    if (finished)
        return false;

    // Ok, so we haven't finished.  Have we started then?
    if (started) {
        // whether we want to continue even when reionization is finished
        // In order to keep outputting meraxes_grids_%d.hdf5
        if (run_globals.params.Flag_OutputGridsPostReion)
            return true;

        // So we have started, but have not previously found to be finished.  Have
        // we now finished though?
        float* xH = run_globals.reion_grids.xH;
        int ReionGridDim = run_globals.params.ReionGridDim;
        int slab_n_real = (int)(run_globals.reion_grids.slab_nix[run_globals.mpi_rank]) * ReionGridDim * ReionGridDim;

        // If not all cells are ionised then reionization is still progressing...
        finished = 1;
        for (int ii = 0; ii < slab_n_real; ii++)
            if (xH[ii] != 0.0) {
                finished = 0;
                break;
            }
    }
    else
        // Here we haven't finished or previously started.  Should we start then?
        if (run_globals.FirstGal != NULL)
            started = 1;

    // At this stage, `started` and `finished` should be set accordingly for each
    // individual core.  Now we need to combine them on all cores.
    MPI_Allreduce(MPI_IN_PLACE, &started, 1, MPI_INT, MPI_LOR, run_globals.mpi_comm);
    run_globals.reion_grids.started = started;
    MPI_Allreduce(MPI_IN_PLACE, &finished, 1, MPI_INT, MPI_LAND, run_globals.mpi_comm);
    run_globals.reion_grids.finished = finished;

    if (started && (!finished))
        return true;
    else
        return false;
}

void filter(fftwf_complex* box, int local_ix_start, int slab_nx, int grid_dim, float R)
{
    int filter_type = run_globals.params.ReionFilterType;
    int middle = grid_dim / 2;
    float box_size = run_globals.params.BoxSize;
    float delta_k = 2.0 * M_PI / box_size;

    // Loop through k-box
    for (int n_x = 0; n_x < slab_nx; n_x++) {
        float k_x;
        int n_x_global = n_x + local_ix_start;

        if (n_x_global > middle)
            k_x = (n_x_global - grid_dim) * delta_k;
        else
            k_x = n_x_global * delta_k;

        for (int n_y = 0; n_y < grid_dim; n_y++) {
            float k_y;

            if (n_y > middle)
                k_y = (n_y - grid_dim) * delta_k;
            else
                k_y = n_y * delta_k;

            for (int n_z = 0; n_z <= middle; n_z++) {
                float k_z = n_z * delta_k;

                float k_mag = sqrtf(k_x * k_x + k_y * k_y + k_z * k_z);

                float kR = k_mag * R; // Real space top-hat

                switch (filter_type) {
                case 0: // Real space top-hat
                    if (kR > 1e-4)
                        box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] *= (fftwf_complex)(3.0 * (sinf(kR) / powf(kR, 3) - cosf(kR) / powf(kR, 2)));
                    break;

                case 1: // k-space top hat
                    kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
                    if (kR > 1)
                        box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] = (fftwf_complex)0.0;
                    break;

                case 2: // Gaussian
                    kR *= 0.643; // Equates integrated volume to the real space top-hat
                    box[grid_index(n_x, n_y, n_z, grid_dim, INDEX_COMPLEX_HERM)] *= (fftwf_complex)(powf(M_E, -kR * kR / 2.0));
                    break;

                default:
                    if ((n_x == 0) && (n_y == 0) && (n_z == 0)) {
                        mlog_error("ReionFilterType.c: Warning, ReionFilterType type %d is undefined!", filter_type);
                        ABORT(EXIT_FAILURE);
                    }
                    break;
                }
            }
        }
    } // End looping through k box
}
