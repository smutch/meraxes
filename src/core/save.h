#ifndef SAVE_H
#define SAVE_H

#include "meraxes.h"

typedef struct galaxy_output_t
{
  long long HaloID;

  // Unique ID for the galaxy
  unsigned long ID;

#ifdef CALC_MAGS
  float Mags[MAGS_N_BANDS];
  float DustyMags[MAGS_N_BANDS];
#if USE_MINI_HALOS
  float MagsIII[MAGS_N_BANDS];
#endif
#endif

  int Type;
  int CentralGal;
  int GhostFlag;
  int Len;
  int MaxLen;

  float Pos[3];
  float Vel[3];
  float Spin;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float FOFMvir;

  // baryonic reservoirs
  float HotGas;
  float MetalsHotGas;
  float ColdGas;
  float MetalsColdGas;
  float H2Frac;
  float H2Mass;
  float HIMass;
  float Mcool;
  float DiskScaleLength;
  float StellarMass;
  float GrossStellarMass;
  float Fesc;
  float FescWeightedGSM;
  float MetalsStellarMass;
  float Sfr;
  float EjectedGas;
  float MetalsEjectedGas;
  float BlackHoleMass;
  float FescBH;
  float BHemissivity;
  float EffectiveBHM;
  float BlackHoleAccretedHotMass;
  float BlackHoleAccretedColdMass;

  int Galaxy_Population; // You need it also if you are not disentangling PopIII/PopII (when Mini_halos is off, this is
                         // = 2)
#if USE_MINI_HALOS
  float GrossStellarMassIII;
  float FescIII;
  float FescIIIWeightedGSM;

  float RmetalBubble;
  int Flag_ExtMetEnr;
  float Metal_Probability;
  float GalMetal_Probability;
  float StellarMass_II;
  float StellarMass_III;
  float Remnant_Mass;
  float MvirCrit_MC;
#endif

  // misc
  float Rcool;
  float Cos_Inc;
  float MergTime;
  float MergerStartRadius;
  float BaryonFracModifier;
  float FOFMvirModifier;
  float MvirCrit;
  float dt;
  float MergerBurstMass;

  // baryonic histories
  float MWMSA; // Mass weighted mean stellar age
  float NewStars[N_HISTORY_SNAPS];
#if USE_MINI_HALOS
  float NewStars_II[N_HISTORY_SNAPS];
  float NewStars_III[N_HISTORY_SNAPS];
#endif
} galaxy_output_t;

#ifdef __cplusplus
extern "C"
{
#endif

  void prepare_galaxy_for_output(struct galaxy_t gal, galaxy_output_t* galout, int i_snap);
  void calc_hdf5_props(void);
  void prep_hdf5_file(void);
  void create_master_file(void);
  void write_snapshot(int n_write, int i_out, int* last_n_write);

#ifdef __cplusplus
}
#endif

#endif
