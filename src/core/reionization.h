#ifndef REIONIZATION_H
#define REIONIZATION_H

#include <fftw3-mpi.h>
#include <stdbool.h>

#include "meraxes.h"
#include "utils.h"

// Constants relevant for spin temperature (taken from 21cmFAST)
#define SIGMA_HI (double)(6.3e-18) /* HI ionization  cross section at 13.6 eV in cm^-2 */
#define TINY (double)(1e-30)

// Note these have the hubble (little h) factor included (taken from 21cmFAST)
#define RHOcrit                                                                                                        \
  (double)((3.0 * HUBBLE * HUBBLE * run_globals.params.Hubble_h * run_globals.params.Hubble_h /                        \
            (8.0 * M_PI * GRAVITY)) *                                                                                  \
           (MPC * MPC * MPC) / SOLAR_MASS) /* Msun Mpc^-3 */ /* at z=0 */
#define RHOcrit_cgs                                                                                                    \
  (double)(3.0 * HUBBLE * HUBBLE * run_globals.params.Hubble_h * run_globals.params.Hubble_h /                         \
           (8.0 * M_PI * GRAVITY)) /* g pcm^-3 */ /* at z=0 */
#define OMb (run_globals.params.BaryonFrac * run_globals.params.OmegaM)
#define No                                                                                                             \
  (double)(RHOcrit_cgs * OMb * (1 - run_globals.params.physics.Y_He) /                                                 \
           PROTONMASS) /*  current hydrogen number density estimate  (#/cm^3)  ~1.92e-7*/
#define He_No                                                                                                          \
  (double)(RHOcrit_cgs * OMb * run_globals.params.physics.Y_He /                                                       \
           (4.0 * PROTONMASS))              /*  current helium number density estimate */
#define f_H (double)(No / (No + He_No))     /* hydrogen number fraction */
#define f_He (double)(He_No / (No + He_No)) /* helium number fraction */
#define FRACT_FLOAT_ERR (double)(1e-7)      /* fractional floating point error */
#define N_b0 (double)(No + He_No)           /* present-day baryon num density, H + He */

#define N_RSD_STEPS (int)(50)

// Parameters taken from 21cmFAST
#define MAX_TK (float)5e4
#define L_FACTOR 0.620350491 // Factor relating cube length to filter radius = (4PI/3)^(-1/3)
#define MAX_DVDR (float)(0.2)

#define alphaB_10k (double)(2.59e-13) /* taken from osterbrock for T=10000 */

typedef struct gal_to_slab_t
{
  int index;
  struct galaxy_t* galaxy;
  int slab_ind;
} gal_to_slab_t;

#ifdef __cplusplus
extern "C"
{
#endif

  void update_galaxy_fesc_vals(struct galaxy_t* gal, double new_stars, int snapshot);
  void set_quasar_fobs(void);
  void set_ReionEfficiency(void);
  void assign_slabs(void);
  void call_find_HII_bubbles(int snapshot, int nout_gals, timer_info* timer);
  void call_ComputeTs(int snapshot, int nout_gals, timer_info* timer);
  void init_reion_grids(void);
  void malloc_reionization_grids(void);
  void free_reionization_grids(void);
  int map_galaxies_to_slabs(int ngals);
  void assign_Mvir_crit_to_galaxies(int ngals_in_slabs, int flag_feed);
  void construct_baryon_grids(int snapshot, int ngals);
  void gen_grids_fname(const int snapshot, char* name, const bool relative);
  void save_reion_input_grids(int snapshot);
  void save_reion_output_grids(int snapshot);
  bool check_if_reionization_ongoing(int snapshot);
  void filter(fftwf_complex* box, int local_ix_start, int slab_nx, int grid_dim, float R, int filter_type);
  void velocity_gradient(fftwf_complex* box, int local_ix_start, int slab_nx, int grid_dim);

#ifdef __cplusplus
}
#endif

#endif
