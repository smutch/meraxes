#ifndef MERGERS_H
#define MERGERS_H

#ifdef __cplusplus
extern "C" {
#endif

double calculate_merging_time(struct galaxy_t* gal, int snapshot);
void merge_with_target(struct galaxy_t* gal, int* dead_gals, int snapshot);

#ifdef __cplusplus
}
#endif

#endif
