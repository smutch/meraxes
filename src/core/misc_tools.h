#ifndef MISC_TOOLS_H
#define MISC_TOOLS_H

#include <stdbool.h>
#include <stdlib.h>

typedef enum index_type {
    INDEX_PADDED = 5674,
    INDEX_REAL,
    INDEX_COMPLEX_HERM,
} index_type;


#ifdef __cplusplus
extern "C" {
#endif

double calc_metallicity(double total_gas, double metals);
int compare_ints(const void* a, const void* b);
int compare_longs(const void* a, const void* b);
int compare_floats(const void* a, const void* b);
int compare_ptrdiff(const void* a, const void* b);
int compare_int_long(const void* a, const void* b);
int compare_slab_assign(const void* a, const void* b);
float apply_pbc_pos(float x);
int searchsorted(void* val,
    void* arr,
    int count,
    size_t size,
    int (*compare)(const void* a, const void* b),
    int imin,
    int imax);
int pos_to_ngp(double x, double side, int nx);
float comoving_distance(float a[3], float b[3]);
double accurate_sumf(float* arr, int n);
int grid_index(int i, int j, int k, int dim, index_type type);

int grid_index_smoothedSFR(int radii, int i, int j, int k, int filter_steps, int dim);
int grid_index_LC(int i, int j, int k, int dim, int dim_LC);

int isclosef(float a, float b, float rel_tol, float abs_tol);
int find_original_index(int index, int* lookup, int n_mappings);
double interp(double xp, double *x, double *y, int nPts);
double trapz_table(double *y, double *x, int nPts, double a, double b);
bool check_for_flag(int flag, int tree_flags);

#ifdef __cplusplus
}
#endif

#endif
