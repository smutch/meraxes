#include "meraxes.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <limits.h>

typedef struct unique_t {
    int* count;
    int max_id;
    int func_switch;
} unique_t;

static herr_t iterate_func(hid_t grp, const char* name, const H5L_info_t* info,
    void* unique)
{
    if (strstr(name, "Snap") != NULL) {
        // this is a snapshot group
        char ds_name[100];
        sprintf(ds_name, "%s/ForestID", name);

        hsize_t n_forests;
        H5LTget_dataset_info(grp, ds_name, &n_forests, NULL, NULL);

        int* forest_ids = NULL;
        if (n_forests > 0) {
            forest_ids = malloc(sizeof(int) * n_forests);
            H5LTread_dataset_int(grp, ds_name, forest_ids);

            if (((unique_t*)unique)->func_switch == 1) {
                int* count = ((unique_t*)unique)->count;
                for (int ii = 0; ii < (int)n_forests; ii++)
                    count[forest_ids[ii] + 1]++;
            }

        } else if (((unique_t*)unique)->func_switch == 0) {
            n_forests = 2;
            forest_ids = malloc(sizeof(int) * n_forests);
            for (int ii = 0; ii < (int)n_forests; ii++)
                forest_ids[ii] = 0;
        }

        if (((unique_t*)unique)->func_switch == 0) {
            qsort(forest_ids, n_forests, sizeof(int), compare_ints);

            int* max_id = &(((unique_t*)unique)->max_id);
            if (forest_ids[n_forests - 1] > (int)(*max_id))
                *max_id = forest_ids[n_forests - 1];
        }

        if (forest_ids != NULL)
            free(forest_ids);
    }
    return 0;
}

int* read_forest_ids(hid_t fd)
{
    // open the file
    hid_t grp = H5Gopen(fd, "/", H5P_DEFAULT);

    // iterate through all present snapshots and find the maximum forest_id
    unique_t unique;
    unique.func_switch = 0;
    unique.max_id = INT_MIN;
    hsize_t idx = 0;
    H5Literate(grp, H5_INDEX_NAME, H5_ITER_INC, &idx, iterate_func, &unique);
    mlog("Maximum forest ID: %d", MLOG_MESG, unique.max_id);

    // now go back through and count the occurance of each ID
    int n_ids = unique.max_id + 2; // N.B. id=-1 values
    unique.count = malloc(sizeof(int) * n_ids);
    for (int ii = 0; ii < n_ids; ii++)
        unique.count[ii] = 0;
    idx = 0;
    unique.func_switch = 1;
    H5Literate(grp, H5_INDEX_NAME, H5_ITER_INC, &idx, iterate_func, &unique);

    int n_unique = 0;
    for (int ii = 0; ii < n_ids; ii++)
        if (unique.count[ii] > 0)
            n_unique++;

    mlog("N unique forest IDs: %d", MLOG_MESG, n_unique);

    int* forest_ids = malloc(sizeof(int) * n_unique);
    for (int ii = 0, jj = 0; ii < n_ids; ii++)
        if (unique.count[ii] > 0)
            forest_ids[jj++] = ii - 1;

    free(unique.count);

    // close the file
    H5Gclose(grp);

    return forest_ids;
}
