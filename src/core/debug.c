#include "meraxes.h"
#include <assert.h>
#include <unistd.h>

void check_mhysa_pointer()
{
    fprintf(stderr, "Addresses: run_globals.mhysa_self = %p\n", run_globals.mhysa_self);
    fprintf(stderr, "          &run_globals = %p\n", &run_globals);
    fprintf(stderr, "          &run_globals.FNameOut = %p\n", &run_globals.FNameOut);
    fprintf(stderr, "          &run_globals.AA = %p\n", &run_globals.AA);
    fprintf(stderr, "          &run_globals.ZZ = %p\n", &run_globals.ZZ);
    fprintf(stderr, "          &run_globals.LTTime = %p\n", &run_globals.LTTime);
    fprintf(stderr, "          &run_globals.RequestedForestId = %p\n", &run_globals.RequestedForestId);
    fprintf(stderr, "          &run_globals.ListOutputSnaps = %p\n", &run_globals.ListOutputSnaps);
    fprintf(stderr, "          &run_globals.SnapshotIndexLookup = %p\n", &run_globals.SnapshotIndexLookup);
    fprintf(stderr, "          &run_globals.SnapshotDeltax = %p\n", &run_globals.SnapshotDeltax);
    fprintf(stderr, "          &run_globals.FirstGal = %p\n", &run_globals.FirstGal);
    fprintf(stderr, "          &run_globals.LastGal = %p\n", &run_globals.LastGal);
    fprintf(stderr, "          &run_globals.mhysa_self = %p\n", &run_globals.mhysa_self);
}

void mpi_debug_here()
{
#ifdef DEBUG
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Task %d, PID %d on %s ready for attach\n", run_globals.mpi_rank, getpid(), hostname);
    printf("Once connected go up stack to 'sleep(5)' and 'set var i=7'\n");
    fflush(stdout);
    // cppcheck-suppress knownConditionTrueFalse
    while (0 == i)
        sleep(5);
#endif
}

static void find_missing_gals(fof_group_t* fof_group, int NFof, int flag)
{
    galaxy_t* gal = NULL;
    halo_t* halo = NULL;
    int counter = 0;
    int master_counter = 0;
    int missing_counter = 0;
    bool* gal_found = NULL;
    galaxy_t** missing_pointers;

    mlog("Running `find_missing_galaxies`...", MLOG_OPEN);

    // If flag =0 then we have more galaxies in the global linked list than when
    // we traverse the FOF groups.  If flag =1 then we have the opposite
    // situation...

    if (flag == 0) {
        // Count the number of galaxies
        gal = run_globals.FirstGal;
        while (gal != NULL) {
            gal->output_index = counter++;
            gal = gal->Next;
        }

        gal_found = calloc((size_t)counter, sizeof(bool));

        // Loop through each FOF halo and mark off each galaxy
        for (int i_fof = 0; i_fof < NFof; i_fof++) {
            halo = fof_group[i_fof].FirstHalo;
            while (halo != NULL) {
                gal = halo->Galaxy;
                while (gal != NULL) {
                    gal_found[gal->output_index] = true;
                    gal = gal->NextGalInHalo;
                }
                halo = halo->NextHaloInFOFGroup;
            }
        }
    } else if (flag == 1) {
        // Count the number of galaxies
        for (int i_fof = 0; i_fof < NFof; i_fof++) {
            halo = fof_group[i_fof].FirstHalo;
            while (halo != NULL) {
                gal = halo->Galaxy;
                while (gal != NULL) {
                    gal->output_index = counter++;
                    gal = gal->NextGalInHalo;
                }
                halo = halo->NextHaloInFOFGroup;
            }
        }

        mlog("I find counter=%d using FOF traversal...", MLOG_MESG, counter);

        gal_found = malloc(sizeof(bool) * counter);
        for (int ii = 0; ii < counter; ii++)
            gal_found[ii] = false;

        // Traverse the global linked list and mark off each galaxy
        gal = run_globals.FirstGal;
        master_counter = 0;
        while (gal != NULL) {
            if (!gal->ghost_flag)
                gal_found[gal->output_index] = true;
            master_counter++;
            gal = gal->Next;
        }
        mlog("I find %d gals traversing global list...", MLOG_MESG, master_counter);
    }

    // Now create an array which holds pointers to the missing galaxies
    for (int ii = 0; ii < counter; ii++)
        if (!gal_found[ii]) {
            mlog("ii = %d", MLOG_MESG, ii);
            missing_counter++;
        }

    master_counter = counter;

    // Check the number of gals with ghost_flag=true
    gal = run_globals.FirstGal;
    counter = 0;
    while (gal != NULL) {
        if (gal->ghost_flag)
            counter++;
        gal = gal->Next;
    }
    mlog("I find %d gals with ghost_flag=true", MLOG_MESG, counter);
    counter = 0;
    for (int i_fof = 0; i_fof < NFof; i_fof++) {
        halo = fof_group[i_fof].FirstHalo;
        while (halo != NULL) {
            gal = halo->Galaxy;
            while (gal != NULL) {
                if (gal->ghost_flag)
                    counter++;
                gal = gal->NextGalInHalo;
            }
            halo = halo->NextHaloInFOFGroup;
        }
    }
    mlog("I find %d gals with ghost_flag=true (FOF traversal)", MLOG_MESG, counter);

    missing_pointers = calloc((size_t)missing_counter, sizeof(galaxy_t*));

    // Loop through the galaxies and store the pointers of the missing ones
    counter = 0;
    if (flag == 0) {
        gal = run_globals.FirstGal;
        while (gal != NULL) {
            // Note that we only store non-ghost missing pointers here...
            if ((!gal_found[gal->output_index]) && (gal->SnapSkipCounter <= 0))
                missing_pointers[counter++] = gal;
            gal = gal->Next;
        }
    } else if (flag == 1) {
        for (int ii = 0; ii < master_counter; ii++)
            if (!gal_found[ii])
                for (int i_fof = 0; i_fof < NFof; i_fof++) {
                    halo = fof_group[i_fof].FirstHalo;
                    while (halo != NULL) {
                        gal = halo->Galaxy;
                        while (gal != NULL) {
                            if (gal->output_index == ii)
                                missing_pointers[counter++] = gal;
                            gal = gal->NextGalInHalo;
                        }
                        halo = halo->NextHaloInFOFGroup;
                    }
                }
    }

    mpi_debug_here();

    free(missing_pointers);
    free(gal_found);

    mlog("...done", MLOG_CLOSE);
}

void check_counts(fof_group_t* fof_group, int NGal, int NFof)
{
    int counter = 0;
    int gal_next_counter = 0;
    int halo_counter = 0;
    int halo_pop_count = 0;
    int total_NGal = 0;
    int total_NFof = 0;
    int total_NGhosts = 0;
    galaxy_t* gal = run_globals.FirstGal;
    halo_t* halo = NULL;

    mlog("Running counts check...", MLOG_OPEN | MLOG_TIMERSTART);

    MPI_Allreduce(&NFof, &total_NFof, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(&NGal, &total_NGal, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(&(run_globals.NGhosts), &total_NGhosts, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    mlog("NFof = %d", MLOG_MESG, total_NFof);
    mlog("NGal = %d", MLOG_MESG, total_NGal);
    mlog("NGhosts = %d", MLOG_MESG, total_NGhosts);

    counter = 0;
    while (gal != NULL) {
        counter++;
        gal = gal->Next;
    }
    MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    mlog("Counting using gal->Next gives %d gals (-%d ghosts = %d gals)",
        MLOG_MESG, counter, total_NGhosts,
        counter - total_NGhosts);
    gal_next_counter = counter;

    halo_pop_count = 0;
    counter = 0;
    halo_counter = 0;
    if (NGal > 0) {
        int ii;
        for (int i_fof = 0; i_fof < NFof; i_fof++) {
            int jj = 0;
            halo = fof_group[i_fof].FirstHalo;
            while (halo != NULL) {
                gal = halo->Galaxy;
                if (gal != NULL)
                    halo_pop_count++;
                ii = 0;
                while (gal != NULL) {
                    gal = gal->NextGalInHalo;
                    counter++;
                    ii++;
                    if (ii > 1e4)
                        ABORT(EXIT_FAILURE);
                }
                halo = halo->NextHaloInFOFGroup;
                halo_counter++;
                jj++;
                if (jj > 1e5)
                    ABORT(EXIT_FAILURE);
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &halo_counter, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE, &halo_pop_count, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    mlog("Counting using FOF groups gives %d gals in %d halos", MLOG_MESG, counter, halo_counter);
    mlog("%d halos are populated with at least one galaxy", MLOG_MESG, halo_pop_count);

    if ((gal_next_counter - total_NGhosts) != counter) {
        int flag;
        if ((gal_next_counter - total_NGhosts) > counter)
            flag = 0;
        else
            flag = 1;
        find_missing_gals(fof_group, NFof, flag);
        ABORT(EXIT_FAILURE);
    }

    mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);
}

void check_pointers(halo_t* halos, fof_group_t* fof_groups, trees_info_t* trees_info)
{
    galaxy_t *gal, gal_deref;
    halo_t* halo;
    int n_halos = trees_info->n_halos;
    int n_fof_groups = trees_info->n_fof_groups;

    mlog("Running pointers check.  Remember to run with Valgrind if you want to check the ->galaxy pointers.", MLOG_MESG);

    gal = run_globals.FirstGal;
    while (gal != NULL) {
        if (!gal->ghost_flag) {
            halo = gal->Halo;
            assert((halo - halos) < (long)n_halos);
        }
        galaxy_t* gal_pointer = gal->FirstGalInHalo;
        if (gal_pointer != NULL)
            gal_deref = *gal_pointer;
        gal_pointer = gal->NextGalInHalo;
        if (gal_pointer != NULL)
            gal_deref = *gal_pointer;
        gal_pointer = gal->Next;
        if (gal_pointer != NULL)
            gal_deref = *gal_pointer;
        if (gal->Type == 3) {
            gal_pointer = gal->MergerTarget;
            if (gal_pointer != NULL)
                gal_deref = *gal_pointer;
        }

        gal = gal->Next;
    }

    for (int ii = 0; ii < n_halos; ii++) {
        fof_group_t* fof_group = halos[ii].FOFGroup;
        // mlog("%ld < %ld\n", MLOG_MESG, fof_group-fof_groups, n_fof_groups);
        assert((fof_group - fof_groups) < (long)n_fof_groups);
        halo = halos[ii].NextHaloInFOFGroup;
        if (halo != NULL)
            // mlog("%llu < %llu", MLOG_MESG, halo-halos, (size_t)n_halos);
            assert((halo - halos) < (long)n_halos);
        gal = halos[ii].Galaxy;
        if (gal != NULL)
            gal_deref = *gal;
    }

    for (int ii = 0; ii < n_fof_groups; ii++) {
        halo = fof_groups[ii].FirstHalo;
        assert((halo - halos) < (long)n_halos);
        halo = fof_groups[ii].FirstOccupiedHalo;
        if (halo != NULL)
            assert((halo - halos) < (long)n_halos);
    }
}

void write_single_grid(const char* fname,
    float* grid,
    int local_ix_start,
    int local_nix,
    int dim,
    const char* grid_name,
    bool padded_flag,
    bool create_file_flag)
{
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

    H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);

    hid_t fd;
    if (create_file_flag)
        fd = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    else
        fd = H5Fopen(fname, H5F_ACC_RDWR, plist_id);

    H5Pclose(plist_id);

    float* grid_out;

    if (padded_flag) {
        grid_out = malloc(local_nix * dim * dim * sizeof(float));
        for (int ii = 0; ii < local_nix; ii++)
            for (int jj = 0; jj < dim; jj++)
                for (int kk = 0; kk < dim; kk++)
                    grid_out[grid_index(ii, jj, kk, dim, INDEX_REAL)] = grid[grid_index(ii, jj, kk, dim, INDEX_PADDED)];
    } else
        grid_out = grid;

    // create the filespace
    hsize_t dims[3] = { (hsize_t)dim, (hsize_t)dim, (hsize_t)dim };
    hid_t fspace_id = H5Screate_simple(3, dims, NULL);

    // create the memspace
    hsize_t mem_dims[3] = { (hsize_t)local_nix, (hsize_t)dim, (hsize_t)dim };
    hid_t memspace_id = H5Screate_simple(3, mem_dims, NULL);

    // select a hyperslab in the filespace
    hsize_t start[3] = { (hsize_t)local_ix_start, 0, 0 };
    hsize_t count[3] = { (hsize_t)local_nix, (hsize_t)dim, (hsize_t)dim };
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

    // crerate the dataset
    hid_t dset_id = H5Dcreate(fd, grid_name, H5T_NATIVE_FLOAT, fspace_id,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // write the dataset
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_id, fspace_id, plist_id, grid_out);

    // cleanup
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(memspace_id);
    H5Sclose(fspace_id);
    H5Fclose(fd);

    if (padded_flag)
        free(grid_out);
}
