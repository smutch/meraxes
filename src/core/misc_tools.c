#include "meraxes.h"
#include <assert.h>
#include <math.h>

void myexit(int signum)
{
    fprintf(stderr, "Task: %d\tis exiting.\n\n\n", run_globals.mpi_rank);
    cn_quote();
    mpi_debug_here();
    cleanup();
    MPI_Finalize();
    exit(signum);
}

double calc_metallicity(double total_gas, double metals)
{
    double Z;

    if ((total_gas > 0) && (metals > 0))
        Z = metals / total_gas;
    else
        Z = 0.0;

    if (Z < 0)
        Z = 0.0;
    if (Z > 1)
        Z = 1.0;

    return Z;
}

int compare_ints(const void* a, const void* b)
{
    return *((int*)a) - *((int*)b);
}
int compare_floats(const void* a, const void* b) {
    float value = *(float*)a - *(float*)b;

    if (value > 0)
        return 1;
    else if (value < 0)
        return -1;
    else
        return 0;
}

int compare_ptrdiff(const void* a, const void* b)
{
    ptrdiff_t result = *(ptrdiff_t*)a - *(ptrdiff_t*)b;

    return (int)result;
}

int compare_slab_assign(const void* a, const void* b)
{
    int value = ((gal_to_slab_t*)a)->slab_ind - ((gal_to_slab_t*)b)->slab_ind;

    return value != 0 ? value : ((gal_to_slab_t*)a)->index - ((gal_to_slab_t*)b)->index;
}

static inline float apply_pbc_disp(float delta)
{
    float box_size = (float)(run_globals.params.BoxSize);

    if (fabs(delta - box_size) < fabs(delta))
        delta -= box_size;
    if (fabs(delta + box_size) < fabs(delta))
        delta += box_size;

    return delta;
}

float apply_pbc_pos(float x)
{
    float box_size = (float)(run_globals.params.BoxSize);

    if (x >= box_size)
        x -= box_size;
    else if (x < 0.0)
        x += box_size;

    return x;
}

int searchsorted(void* val,
    void* arr,
    int count,
    size_t size,
    int (*compare)(const void* a, const void* b),
    int imin,
    int imax)
{
    // check if we need to init imin and imax
    if ((imax < 0) && (imin < 0)) {
        imin = 0;
        imax = count - 1;
    }

    // test if we have found the result
    if ((imax - imin) < 0)
        return imax;
    else {
        // calculate midpoint to cut set in half
        int imid = imin + ((imax - imin) / 2);
        void* arr_val = (void*)(((char*)arr + imid * size));

        // three-way comparison
        if (compare(arr_val, val) > 0)
            // key is in lower subset
            return searchsorted(val, arr, count, size, compare, imin, imid - 1);
        else if (compare(arr_val, val) < 0)
            // key is in upper subset
            return searchsorted(val, arr, count, size, compare, imid + 1, imax);
        else
            // key has been found
            return imid;
    }
}

int pos_to_ngp(double x, double side, int nx)
{
    int ind = (int)nearbyint(x / side * (double)nx);

    if (ind > nx - 1)
        ind = 0;

    assert(ind > -1);

    return ind;
}

float comoving_distance(float a[3], float b[3])
{
    float dx = apply_pbc_disp(a[0] - b[0]);
    float dy = apply_pbc_disp(a[1] - b[1]);
    float dz = apply_pbc_disp(a[2] - b[2]);

    float dist = sqrtf(dx * dx + dy * dy + dz * dz);

    assert(dist <= (sqrtf(3.0) / 2.0 * run_globals.params.BoxSize));

    return dist;
}

double accurate_sumf(float* arr, int n)
{
    // inplace reorder and sum
    qsort(arr, n, sizeof(float), compare_floats);

    double total = 0;
    for (int ii = 0; ii < n; ii++)
        total += (double)(arr[ii]);

    return total;
}

int grid_index(int i, int j, int k, int dim, int type)
{
    int ind;

    switch (type) {
    case INDEX_PADDED:
        ind = k + (2 * (dim / 2 + 1)) * (j + dim * i);
        break;
    case INDEX_REAL:
        ind = k + dim * (j + dim * i);
        break;
    case INDEX_COMPLEX_HERM:
        ind = k + (dim / 2 + 1) * (j + dim * i);
        break;
    }

    return ind;
}

/// Numpy style isclose()
int isclosef(
    float a,
    float b,
    float rel_tol, ///< [in] = -1 for Numpy default
    float abs_tol) ///< [in] = -1 for Numpy default
{
    if (abs_tol < 0)
        abs_tol = 1e-8; ///< Numpy default
    if (rel_tol < 0)
        rel_tol = 1e-5; ///< Numpy default
    return fabs(a - b) <= (abs_tol + rel_tol * fabs(b));
}


int find_original_index(int index, int* lookup, int n_mappings)
{
    int new_index = -1;
    int* pointer = bsearch(&index, lookup, (size_t)n_mappings, sizeof(int), compare_ints);

    if (pointer)
        new_index = (int)(pointer - lookup);

    return new_index;
}

