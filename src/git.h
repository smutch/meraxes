/*
            * This file has been automatically generated.  Any modifications will be
            * overwritten during compilation.
            */

#define GITREF_STR "2cdf12c83c10b08edcd9710f7146e3d0bba72170"

#define GITDIFF_STR "diff --git a/src/core/save.c b/src/core/save.c\n"\
"index 9d0fd55..ae6a901 100644\n"\
"--- a/src/core/save.c\n"\
"+++ b/src/core/save.c\n"\
"@@ -419,6 +419,7 @@ void calc_hdf5_props()\n"\
"     h5props->field_types[i++]   = h5props->array_nhist_f_tid;\n"\
" \n"\
"     // Blackhole or Emissivity related\n"\
"+    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Stellaremissivity);\n"\
"     h5props->dst_field_sizes[i] = sizeof(galout.Stellaremissivity);\n"\
"     h5props->field_names[i]     = \"Stellaremissivity\";\n"\
"     h5props->field_units[i]     = \"1e60 photons\";\n"\
"diff --git a/src/git.h b/src/git.h\n"\
"index a70501e..6445b3d 100644\n"\
"--- a/src/git.h\n"\
"+++ b/src/git.h\n"\
"@@ -3,691 +3,717 @@\n"\
"             * overwritten during compilation.\n"\
"             */\n"\
" \n"\
"-#define GITREF_STR \"536fe2ee35656425b9183cf74ae4614c0bfdd60e\"\n"\
"+#define GITREF_STR \"2cdf12c83c10b08edcd9710f7146e3d0bba72170\"\n"\
" \n"\
"-#define GITDIFF_STR \"diff --git a/include/meraxes.h b/include/meraxes.h\\n\"\\\n"\
"-\"index b700754..e3749ba 100644\\n\"\\\n"\
"-\"--- a/include/meraxes.h\\n\"\\\n"\
"-\"+++ b/include/meraxes.h\\n\"\\\n"\
"-\"@@ -1,7 +1,6 @@\\n\"\\\n"\
"-\" // These defines were added at compilation time\\n\"\\\n"\
"-\" // (see SConstruct prepend_user_defines())\\n\"\\\n"\
"-\"-#define NDEBUG\\n\"\\\n"\
"-\"-#define USE_MPI 1\\n\"\\\n"\
"-\"+#define U\\n\"\\\n"\
"-\" // -------------------------------------------\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"@@ -108,7 +107,12 @@ typedef struct physics_params_t {\\n\"\\\n"\
"-\"   double EnergyPerSN;\\n\"\\\n"\
"-\"   double IMFNormConst;\\n\"\\\n"\
"-\"   double RadioModeEff;\\n\"\\\n"\
"-\"+  double QuasarModeEff;\\n\"\\\n"\
"-\"   double BlackHoleGrowthRate;\\n\"\\\n"\
"-\"+  double EddingtonRatio;\\n\"\\\n"\
"-\"+  double quasar_mode_scaling;\\n\"\\\n"\
"-\"+  double quasar_open_angel;\\n\"\\\n"\
"-\"+  double quasar_fobs;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   double ThreshMajorMerger;\\n\"\\\n"\
"-\"   double MinMergerStellarMass;\\n\"\\\n"\
"-\"@@ -120,12 +124,16 @@ typedef struct physics_params_t {\\n\"\\\n"\
"-\"   // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value\\n\"\\\n"\
"-\"   double ReionEfficiency;\\n\"\\\n"\
"-\"   double ReionNionPhotPerBary;\\n\"\\\n"\
"-\"-  double ReionEscapeFrac;\\n\"\\\n"\
"-\"+  double ReionEscapeFrac;  \\n\"\\\n"\
"-\"+  double ReionEscapeFracBH;\\n\"\\\n"\
"-\"+  double BlackHoleSeed;\\n\"\\\n"\
"-\"+  double BlackHoleMassLimitReion;\\n\"\\\n"\
"-\"   double ReionTcool;\\n\"\\\n"\
"-\"   double Y_He;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   double  ReionGammaHaloBias;\\n\"\\\n"\
"-\"   double  ReionAlphaUV;\\n\"\\\n"\
"-\"+  double  ReionAlphaUVBH;\\n\"\\\n"\
"-\"   double  ReionRBubbleMin;\\n\"\\\n"\
"-\"   double  ReionRBubbleMax;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"@@ -147,9 +155,13 @@ typedef struct physics_params_t {\\n\"\\\n"\
"-\"   double  ReionSMParam_d;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   // Flags\\n\"\\\n"\
"-\"-  int Flag_RedshiftDepEscFrac;\\n\"\\\n"\
"-\"+  double RedshiftDepEscFracNorm;\\n\"\\\n"\
"-\"+  double RedshiftDepEscFracScaling;\\n\"\\\n"\
"-\"+  double RedshiftDepEscFracBHNorm;\\n\"\\\n"\
"-\"+  double RedshiftDepEscFracBHScaling;\\n\"\\\n"\
"-\"   int Flag_ReionizationModifier;\\n\"\\\n"\
"-\"   int Flag_BHFeedback;\\n\"\\\n"\
"-\"+  int Flag_BHReion;\\n\"\\\n"\
"-\"   int Flag_IRA;\\n\"\\\n"\
"-\"   int Flag_FixDiskRadiusOnInfall;\\n\"\\\n"\
"-\"   int Flag_FixVmaxOnInfall;\\n\"\\\n"\
"-\"@@ -176,6 +188,8 @@ typedef struct run_params_t {\\n\"\\\n"\
"-\"   char             MagBands[STRLEN];\\n\"\\\n"\
"-\"   char             ForestIDFile[STRLEN];\\n\"\\\n"\
"-\"   char             MvirCritFile[STRLEN];\\n\"\\\n"\
"-\"+  char             MassRatioModifier[STRLEN];\\n\"\\\n"\
"-\"+  char             BaryonFracModifier[STRLEN];\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   physics_params_t physics;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"@@ -214,6 +228,7 @@ typedef struct run_params_t {\\n\"\\\n"\
"-\"   int              FlagReadDumpFile;\\n\"\\\n"\
"-\"   int              FlagMCMC;\\n\"\\\n"\
"-\"   int              Flag_PatchyReion;\\n\"\\\n"\
"-\"+  int              Flag_output_grids;\\n\"\\\n"\
"-\" } run_params_t;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"@@ -373,12 +388,20 @@ typedef struct galaxy_t {\\n\"\\\n"\
"-\"   double Mcool;\\n\"\\\n"\
"-\"   double StellarMass;\\n\"\\\n"\
"-\"   double GrossStellarMass;\\n\"\\\n"\
"-\"+  double FescWeightedGSM;\\n\"\\\n"\
"-\"+  double Stellaremissivity;\\n\"\\\n"\
"-\"+  double MergerSemissivity;\\n\"\\\n"\
"-\"   double MetalsStellarMass;\\n\"\\\n"\
"-\"   double DiskScaleLength;\\n\"\\\n"\
"-\"   double Sfr;\\n\"\\\n"\
"-\"   double EjectedGas;\\n\"\\\n"\
"-\"   double MetalsEjectedGas;\\n\"\\\n"\
"-\"   double BlackHoleMass;\\n\"\\\n"\
"-\"+  double BHemissivity;\\n\"\\\n"\
"-\"+  double EffectiveBHM;\\n\"\\\n"\
"-\"+  double BlackHoleAccretedHotMass;\\n\"\\\n"\
"-\"+  double BlackHoleAccretedColdMass;\\n\"\\\n"\
"-\"+  double BlackHoleAccretingColdMass;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   // baryonic hostories\\n\"\\\n"\
"-\"   double mwmsa_num;\\n\"\\\n"\
"-\"@@ -444,11 +467,18 @@ typedef struct galaxy_output_t {\\n\"\\\n"\
"-\"   float DiskScaleLength;\\n\"\\\n"\
"-\"   float StellarMass;\\n\"\\\n"\
"-\"   float GrossStellarMass;\\n\"\\\n"\
"-\"+  float Stellaremissivity;\\n\"\\\n"\
"-\"+  float MergerSemissivity;\\n\"\\\n"\
"-\"+  float FescWeightedGSM;\\n\"\\\n"\
"-\"   float MetalsStellarMass;\\n\"\\\n"\
"-\"   float Sfr;\\n\"\\\n"\
"-\"   float EjectedGas;\\n\"\\\n"\
"-\"   float MetalsEjectedGas;\\n\"\\\n"\
"-\"   float BlackHoleMass;\\n\"\\\n"\
"-\"+  float BHemissivity;\\n\"\\\n"\
"-\"+  float EffectiveBHM;\\n\"\\\n"\
"-\"+  float BlackHoleAccretedHotMass;\\n\"\\\n"\
"-\"+  float BlackHoleAccretedColdMass;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   // misc\\n\"\\\n"\
"-\"   float Rcool;\\n\"\\\n"\
"-\"@@ -457,6 +487,7 @@ typedef struct galaxy_output_t {\\n\"\\\n"\
"-\"   float MergerStartRadius;\\n\"\\\n"\
"-\"   float BaryonFracModifier;\\n\"\\\n"\
"-\"   float MvirCrit;\\n\"\\\n"\
"-\"+  float dt;\\n\"\\\n"\
"-\"   float MergerBurstMass;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   // baryonic histories\\n\"\\\n"\
"-\"@@ -511,6 +542,17 @@ typedef struct catalog_halo_t {\\n\"\\\n"\
"-\"   char      padding[8];                //!< Alignment padding\\n\"\\\n"\
"-\" } catalog_halo_t;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"+typedef struct Modifier\\n\"\\\n"\
"-\"+{\\n\"\\\n"\
"-\"+    float logMmin;\\n\"\\\n"\
"-\"+    float logMmax;\\n\"\\\n"\
"-\"+    float mass_mean;\\n\"\\\n"\
"-\"+    float mass_errl;\\n\"\\\n"\
"-\"+    float mass_erru;\\n\"\\\n"\
"-\"+    float ratio;\\n\"\\\n"\
"-\"+    float ratio_errl;\\n\"\\\n"\
"-\"+    float ratio_erru;\\n\"\\\n"\
"-\"+} Modifier;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" //! Global variables which will will be passed around\\n\"\\\n"\
"-\" typedef struct run_globals_t {\\n\"\\\n"\
"-\"@@ -523,7 +565,10 @@ typedef struct run_globals_t {\\n\"\\\n"\
"-\"   double             *AA;\\n\"\\\n"\
"-\"   double             *ZZ;\\n\"\\\n"\
"-\"   double             *LTTime;\\n\"\\\n"\
"-\"+  double             *mass_weighted_xHII;\\n\"\\\n"\
"-\"   int                *RequestedForestId;\\n\"\\\n"\
"-\"+  int                RequestedMassRatioModifier;\\n\"\\\n"\
"-\"+  int                RequestedBaryonFracModifier;\\n\"\\\n"\
"-\"   int                *ListOutputSnaps;\\n\"\\\n"\
"-\"   halo_t            **SnapshotHalo;\\n\"\\\n"\
"-\"   fof_group_t       **SnapshotFOFGroup;\\n\"\\\n"\
"-\"@@ -550,6 +595,8 @@ typedef struct run_globals_t {\\n\"\\\n"\
"-\"   int                 NStoreSnapshots;\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"   bool                SelectForestsSwitch;\\n\"\\\n"\
"-\"+  Modifier            *mass_ratio_modifier;\\n\"\\\n"\
"-\"+  Modifier            *baryon_frac_modifier;\\n\"\\\n"\
"-\" } run_globals_t;\\n\"\\\n"\
"-\" #ifdef _MAIN\\n\"\\\n"\
"-\" run_globals_t run_globals;\\n\"\\\n"\
"-\"@@ -615,7 +662,7 @@ float        apply_pbc_pos(float x);\\n\"\\\n"\
"-\" double       accurate_sumf(float *arr, int n);\\n\"\\\n"\
"-\" int          grid_index(int i, int j, int k, int dim, int type);\\n\"\\\n"\
"-\" void         mpi_debug_here(void);\\n\"\\\n"\
"-\"-int 				 isclosef(float a, float b, float rel_tol, float abs_tol);\\n\"\\\n"\
"-\"+int          isclosef(float a, float b, float rel_tol, float abs_tol);\\n\"\\\n"\
"-\" void         printProgress (double percentage);\\n\"\\\n"\
"-\" void         check_counts(fof_group_t *fof_group, int NGal, int NFof);\\n\"\\\n"\
"-\" void         cn_quote(void);\\n\"\\\n"\
"-\"@@ -626,14 +673,19 @@ double       calculate_Mvir(double Mvir, int len);\\n\"\\\n"\
"-\" double       calculate_Rvir(double Mvir, int snapshot);\\n\"\\\n"\
"-\" double       calculate_Vvir(double Mvir, double Rvir);\\n\"\\\n"\
"-\" double       calculate_spin_param(halo_t *halo);\\n\"\\\n"\
"-\"+void         read_mass_ratio_modifiers(int snapshot);\\n\"\\\n"\
"-\"+void         read_baryon_frac_modifiers(int snapshot);\\n\"\\\n"\
"-\"+double       interpolate_modifier(Modifier *modifier_data, double logM);\\n\"\\\n"\
"-\" void         read_cooling_functions(void);\\n\"\\\n"\
"-\" double       interpolate_cooling_rate(double logTemp, double logZ);\\n\"\\\n"\
"-\" double       gas_cooling(galaxy_t *gal);\\n\"\\\n"\
"-\" void         cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass);\\n\"\\\n"\
"-\" double       calc_metallicity(double total_gas, double metals);\\n\"\\\n"\
"-\"+double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass, double x);\\n\"\\\n"\
"-\"+void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio, int snapshot);\\n\"\\\n"\
"-\"+void         previous_merger_driven_BH_growth(galaxy_t *gal);\\n\"\\\n"\
"-\"+double       calculate_BHemissivity(double BlackHoleMass, double accreted_mass);\\n\"\\\n"\
"-\" void         reincorporate_ejected_gas(galaxy_t *gal);\\n\"\\\n"\
"-\"-double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass);\\n\"\\\n"\
"-\"-void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio);\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" // Magnitude related\\n\"\\\n"\
"-\" void   init_luminosities(galaxy_t *gal);\\n\"\\\n"\
"-\"@@ -653,6 +705,7 @@ void   assign_slabs(void);\\n\"\\\n"\
"-\" void   init_reion_grids(void);\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" void   filter(fftwf_complex *box, int local_ix_start, int slab_nx, int grid_dim, float R);\\n\"\\\n"\
"-\"+void   set_quasar_fobs(void);\\n\"\\\n"\
"-\" void   find_HII_bubbles(float redshift);\\n\"\\\n"\
"-\" double tocf_modifier(galaxy_t *gal, double Mvir);\\n\"\\\n"\
"-\" void   set_ReionEfficiency(void);\\n\"\\\n"\
"-\"diff --git a/src/SConstruct b/src/SConstruct\\n\"\\\n"\
"-\"index 26fc42d..3ff3daa 100644\\n\"\\\n"\
"-\"--- a/src/SConstruct\\n\"\\\n"\
"-\"+++ b/src/SConstruct\\n\"\\\n"\
"-\"@@ -236,6 +236,7 @@ if not GetOption('help'):\\n\"\\\n"\
"-\"                 env.AppendUnique(CCFLAGS = ['-I'+dep['inclp']])\\n\"\\\n"\
"-\"             if 'libp' in dep and dep['libp'] is not None:\\n\"\\\n"\
"-\"                 env.AppendUnique(LIBPATH = [dep['libp']])\\n\"\\\n"\
"-\"+                env.AppendUnique(RPATH = [dep['libp']])\\n\"\\\n"\
"-\"             if 'lib' in dep and dep['lib'] is not None:\\n\"\\\n"\
"-\"                 env.AppendUnique(LIBS = [dep['lib']])\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"diff --git a/src/core/modifiers.c b/src/core/modifiers.c\\n\"\\\n"\
"-\"index 9cfc074..71d5025 100644\\n\"\\\n"\
"-\"--- a/src/core/modifiers.c\\n\"\\\n"\
"-\"+++ b/src/core/modifiers.c\\n\"\\\n"\
"-\"@@ -4,11 +4,12 @@\\n\"\\\n"\
"-\" #include <hdf5_hl.h>\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" #define NFIELDS 8\\n\"\\\n"\
"-\"-#define N_START 9\\n\"\\\n"\
"-\"+#define N_START 0\\n\"\\\n"\
"-\" #define N_LOGMS 31\\n\"\\\n"\
"-\" #define DELTA_M 1.0\\n\"\\\n"\
"-\" #define MIN_LOGM 7.5\\n\"\\\n"\
"-\" #define MAX_LOGM 11.5\\n\"\\\n"\
"-\"+#define M_OFFSET 0.5\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" void read_mass_ratio_modifiers(int snapshot){\\n\"\\\n"\
"-\"   if (strlen(run_globals.params.MassRatioModifier) == 0){\\n\"\\\n"\
"-\"diff --git a/src/core/save.c b/src/core/save.c\\n\"\\\n"\
"-\"index 54e57bb..9d0fd55 100644\\n\"\\\n"\
"+#define GITDIFF_STR \"diff --git a/src/core/save.c b/src/core/save.c\\n\"\\\n"\
"+\"index 9d0fd55..ae6a901 100644\\n\"\\\n"\
" \"--- a/src/core/save.c\\n\"\\\n"\
" \"+++ b/src/core/save.c\\n\"\\\n"\
"-\"@@ -418,7 +418,7 @@ void calc_hdf5_props()\\n\"\\\n"\
"-\"     h5props->field_h_conv[i]    = \\\"v/h\\\";\\n\"\\\n"\
"+\"@@ -419,6 +419,7 @@ void calc_hdf5_props()\\n\"\\\n"\
" \"     h5props->field_types[i++]   = h5props->array_nhist_f_tid;\\n\"\\\n"\
" \" \\n\"\\\n"\
"-\"-	// Blackhole or Emissivity related\\n\"\\\n"\
"-\"+    // Blackhole or Emissivity related\\n\"\\\n"\
"+\"     // Blackhole or Emissivity related\\n\"\\\n"\
"+\"+    h5props->dst_offsets[i]     = HOFFSET(galaxy_output_t, Stellaremissivity);\\n\"\\\n"\
" \"     h5props->dst_field_sizes[i] = sizeof(galout.Stellaremissivity);\\n\"\\\n"\
" \"     h5props->field_names[i]     = \\\"Stellaremissivity\\\";\\n\"\\\n"\
" \"     h5props->field_units[i]     = \\\"1e60 photons\\\";\\n\"\\\n"\
"-\"diff --git a/src/deps/gstar.py b/src/deps/gstar.py\\n\"\\\n"\
"-\"index 669a023..3361b75 100644\\n\"\\\n"\
"-\"--- a/src/deps/gstar.py\\n\"\\\n"\
"-\"+++ b/src/deps/gstar.py\\n\"\\\n"\
"-\"@@ -4,7 +4,7 @@ deps = {}\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" deps['exec'] = {\\n\"\\\n"\
"-\"     'mpicc' : '/usr/local/x86_64/gnu/openmpi-1.10.2-psm/bin/mpicc',\\n\"\\\n"\
"-\"-    'cc' : '/usr/local/gcc-5.1.0/bin/gcc',\\n\"\\\n"\
"-\"+	'cc' : '/usr/local/intel-15.3.0/composer_xe_2015.3.187/bin/intel64/icc',\\n\"\\\n"\
"-\"     'git' : '/usr/bin/git',\\n\"\\\n"\
"-\" }\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"@@ -15,8 +15,8 @@ deps['gsl'] = {\\n\"\\\n"\
"-\" }\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\" deps['hdf5'] = {\\n\"\\\n"\
"-\"-    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/include',\\n\"\\\n"\
"-\"-    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/lib',\\n\"\\\n"\
"-\"+    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/include',\\n\"\\\n"\
"-\"+    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/lib',\\n\"\\\n"\
"-\"     'lib'   : ['hdf5', 'hdf5_hl', 'z'],\\n\"\\\n"\
"-\" }\\n\"\\\n"\
"-\" \\n\"\\\n"\
" \"diff --git a/src/git.h b/src/git.h\\n\"\\\n"\
"-\"index 3d521ae..0b21f2a 100644\\n\"\\\n"\
"-\"--- a/src/git.h\\n\"\\\n"\
"-\"+++ b/src/git.h\\n\"\\\n"\
"-\"@@ -3,6 +3,349 @@\\n\"\\\n"\
"-\"             * overwritten during compilation.\\n\"\\\n"\
"-\"             */\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"-#define GITREF_STR \\\"7cd7348aab75d26b87669795b4080bf8b0f6cf0b\\\"\\n\"\\\n"\
"-\"+#define GITREF_STR \\\"536fe2ee35656425b9183cf74ae4614c0bfdd60e\\\"\\n\"\\\n"\
"-\"+\\n\"\\\n"\
"-\"+#define GITDIFF_STR \\\"diff --git a/include/meraxes.h b/include/meraxes.h\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index b700754..e3749ba 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/include/meraxes.h\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ b/include/meraxes.h\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -1,7 +1,6 @@\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" // These defines were added at compilation time\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" // (see SConstruct prepend_user_defines())\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-#define NDEBUG\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-#define USE_MPI 1\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+#define U\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" // -------------------------------------------\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -108,7 +107,12 @@ typedef struct physics_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double EnergyPerSN;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double IMFNormConst;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double RadioModeEff;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double QuasarModeEff;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double BlackHoleGrowthRate;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double EddingtonRatio;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double quasar_mode_scaling;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double quasar_open_angel;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double quasar_fobs;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double ThreshMajorMerger;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double MinMergerStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -120,12 +124,16 @@ typedef struct physics_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double ReionEfficiency;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double ReionNionPhotPerBary;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-  double ReionEscapeFrac;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double ReionEscapeFrac;  \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double ReionEscapeFracBH;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double BlackHoleSeed;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double BlackHoleMassLimitReion;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double ReionTcool;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double Y_He;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double  ReionGammaHaloBias;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double  ReionAlphaUV;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double  ReionAlphaUVBH;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double  ReionRBubbleMin;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double  ReionRBubbleMax;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -147,9 +155,13 @@ typedef struct physics_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double  ReionSMParam_d;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   // Flags\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-  int Flag_RedshiftDepEscFrac;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double RedshiftDepEscFracNorm;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double RedshiftDepEscFracScaling;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double RedshiftDepEscFracBHNorm;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double RedshiftDepEscFracBHScaling;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int Flag_ReionizationModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int Flag_BHFeedback;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  int Flag_BHReion;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int Flag_IRA;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int Flag_FixDiskRadiusOnInfall;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int Flag_FixVmaxOnInfall;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -176,6 +188,8 @@ typedef struct run_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   char             MagBands[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   char             ForestIDFile[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   char             MvirCritFile[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  char             MassRatioModifier[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  char             BaryonFracModifier[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   physics_params_t physics;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -214,6 +228,7 @@ typedef struct run_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int              FlagReadDumpFile;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int              FlagMCMC;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int              Flag_PatchyReion;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  int              Flag_output_grids;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" } run_params_t;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -373,12 +388,20 @@ typedef struct galaxy_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double Mcool;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double StellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double GrossStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double FescWeightedGSM;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double Stellaremissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double MergerSemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double MetalsStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double DiskScaleLength;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double Sfr;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double EjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double MetalsEjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double BlackHoleMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double BHemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double EffectiveBHM;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double BlackHoleAccretedHotMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double BlackHoleAccretedColdMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double BlackHoleAccretingColdMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   // baryonic hostories\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double mwmsa_num;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -444,11 +467,18 @@ typedef struct galaxy_output_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float DiskScaleLength;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float StellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float GrossStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float Stellaremissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float MergerSemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float FescWeightedGSM;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float MetalsStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float Sfr;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float EjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float MetalsEjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float BlackHoleMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float BHemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float EffectiveBHM;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float BlackHoleAccretedHotMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float BlackHoleAccretedColdMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   // misc\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float Rcool;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -457,6 +487,7 @@ typedef struct galaxy_output_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float MergerStartRadius;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float BaryonFracModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float MvirCrit;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  float dt;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   float MergerBurstMass;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   // baryonic histories\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -511,6 +542,17 @@ typedef struct catalog_halo_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   char      padding[8];                //!< Alignment padding\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" } catalog_halo_t;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+typedef struct Modifier\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+{\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float logMmin;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float logMmax;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float mass_mean;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float mass_errl;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float mass_erru;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float ratio;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float ratio_errl;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    float ratio_erru;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+} Modifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" //! Global variables which will will be passed around\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" typedef struct run_globals_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -523,7 +565,10 @@ typedef struct run_globals_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double             *AA;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double             *ZZ;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   double             *LTTime;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  double             *mass_weighted_xHII;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int                *RequestedForestId;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  int                RequestedMassRatioModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  int                RequestedBaryonFracModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int                *ListOutputSnaps;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   halo_t            **SnapshotHalo;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   fof_group_t       **SnapshotFOFGroup;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -550,6 +595,8 @@ typedef struct run_globals_t {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   int                 NStoreSnapshots;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   bool                SelectForestsSwitch;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  Modifier            *mass_ratio_modifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+  Modifier            *baryon_frac_modifier;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" } run_globals_t;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #ifdef _MAIN\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" run_globals_t run_globals;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -615,7 +662,7 @@ float        apply_pbc_pos(float x);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       accurate_sumf(float *arr, int n);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" int          grid_index(int i, int j, int k, int dim, int type);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         mpi_debug_here(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-int 				 isclosef(float a, float b, float rel_tol, float abs_tol);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+int          isclosef(float a, float b, float rel_tol, float abs_tol);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         printProgress (double percentage);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         check_counts(fof_group_t *fof_group, int NGal, int NFof);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         cn_quote(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -626,14 +673,19 @@ double       calculate_Mvir(double Mvir, int len);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       calculate_Rvir(double Mvir, int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       calculate_Vvir(double Mvir, double Rvir);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       calculate_spin_param(halo_t *halo);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+void         read_mass_ratio_modifiers(int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+void         read_baryon_frac_modifiers(int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+double       interpolate_modifier(Modifier *modifier_data, double logM);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         read_cooling_functions(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       interpolate_cooling_rate(double logTemp, double logZ);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       gas_cooling(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double       calc_metallicity(double total_gas, double metals);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass, double x);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio, int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+void         previous_merger_driven_BH_growth(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+double       calculate_BHemissivity(double BlackHoleMass, double accreted_mass);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void         reincorporate_ejected_gas(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" // Magnitude related\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void   init_luminosities(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -653,6 +705,7 @@ void   assign_slabs(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void   init_reion_grids(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void   filter(fftwf_complex *box, int local_ix_start, int slab_nx, int grid_dim, float R);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+void   set_quasar_fobs(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void   find_HII_bubbles(float redshift);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" double tocf_modifier(galaxy_t *gal, double Mvir);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void   set_ReionEfficiency(void);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/src/SConstruct b/src/SConstruct\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 26fc42d..3ff3daa 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/src/SConstruct\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ b/src/SConstruct\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -236,6 +236,7 @@ if not GetOption('help'):\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"                 env.AppendUnique(CCFLAGS = ['-I'+dep['inclp']])\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"             if 'libp' in dep and dep['libp'] is not None:\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"                 env.AppendUnique(LIBPATH = [dep['libp']])\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+                env.AppendUnique(RPATH = [dep['libp']])\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"             if 'lib' in dep and dep['lib'] is not None:\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"                 env.AppendUnique(LIBS = [dep['lib']])\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/src/core/modifiers.c b/src/core/modifiers.c\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 9cfc074..71d5025 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/src/core/modifiers.c\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ b/src/core/modifiers.c\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -4,11 +4,12 @@\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #include <hdf5_hl.h>\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #define NFIELDS 8\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-#define N_START 9\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+#define N_START 0\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #define N_LOGMS 31\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #define DELTA_M 1.0\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #define MIN_LOGM 7.5\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #define MAX_LOGM 11.5\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+#define M_OFFSET 0.5\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" void read_mass_ratio_modifiers(int snapshot){\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"   if (strlen(run_globals.params.MassRatioModifier) == 0){\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/src/core/save.c b/src/core/save.c\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 54e57bb..9d0fd55 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/src/core/save.c\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ b/src/core/save.c\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -418,7 +418,7 @@ void calc_hdf5_props()\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     h5props->field_h_conv[i]    = \\\\\\\"v/h\\\\\\\";\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     h5props->field_types[i++]   = h5props->array_nhist_f_tid;\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-	// Blackhole or Emissivity related\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    // Blackhole or Emissivity related\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     h5props->dst_field_sizes[i] = sizeof(galout.Stellaremissivity);\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     h5props->field_names[i]     = \\\\\\\"Stellaremissivity\\\\\\\";\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     h5props->field_units[i]     = \\\\\\\"1e60 photons\\\\\\\";\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/src/deps/gstar.py b/src/deps/gstar.py\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 669a023..3361b75 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/src/deps/gstar.py\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ b/src/deps/gstar.py\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -4,7 +4,7 @@ deps = {}\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" deps['exec'] = {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     'mpicc' : '/usr/local/x86_64/gnu/openmpi-1.10.2-psm/bin/mpicc',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-    'cc' : '/usr/local/gcc-5.1.0/bin/gcc',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+	'cc' : '/usr/local/intel-15.3.0/composer_xe_2015.3.187/bin/intel64/icc',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     'git' : '/usr/bin/git',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" }\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -15,8 +15,8 @@ deps['gsl'] = {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" }\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" deps['hdf5'] = {\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/include',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/lib',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/include',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/lib',\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"     'lib'   : ['hdf5', 'hdf5_hl', 'z'],\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" }\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/src/git.h b/src/git.h\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"deleted file mode 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 3d521ae..0000000\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/src/git.h\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ /dev/null\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -1,8 +0,0 @@\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-/*\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-            * This file has been automatically generated.  Any modifications will be\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-            * overwritten during compilation.\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-            */\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-#define GITREF_STR \\\\\\\"7cd7348aab75d26b87669795b4080bf8b0f6cf0b\\\\\\\"\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-#define GITDIFF_STR \\\\\\\"\\\\\\\"\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/src/setup_module.sh b/src/setup_module.sh\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 8d60e62..53d1d5a 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/src/setup_module.sh\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ b/src/setup_module.sh\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -1,4 +1,4 @@\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" #!/bin/bash\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\" module purge\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.10.0-openmpi-1.10.2-psm \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\\\\\":\\\\\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\\\\\"'\\\\\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.8.17-openmpi-1.10.2-psm \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+#for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\\\\\":\\\\\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\\\\\"'\\\\\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"diff --git a/tmp b/tmp\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"deleted file mode 100644\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"index 38c69fd..0000000\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"--- a/tmp\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"+++ /dev/null\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"@@ -1,44 +0,0 @@\\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/SConstruct => src/SConstruct                                                                   \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/SConstruct                                                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/SConstruct                                                          \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/cleanup.c => src/core/cleanup.c                                                           \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/cleanup.c                                                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/dracarys.c => src/core/dracarys.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/dracarys.c                                                                                    \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/core/dracarys.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/find_HII_bubbles.c => src/core/find_HII_bubbles.c                                         \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/find_HII_bubbles.c                                                                            \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/core/find_HII_bubbles.c                                             \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/galaxies.c => src/core/galaxies.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/galaxies.c                                                                                    \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/core/galaxies.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/init.c => src/core/init.c                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/init.c                                                                                        \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/core/init.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/read_halos.c => src/core/read_halos.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/read_halos.c                                                                                  \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/read_params.c => src/core/read_params.c                                                   \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/read_params.c                                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/reionization.c => src/core/reionization.c                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/reionization.c                                                                                \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/core/reionization.c                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/core/save.c => src/core/save.c                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/core/save.c                                                                                        \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/core/save.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/deps/gstar.py => src/deps/gstar.py                                                             \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/deps/gstar.py                                                                                      \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/meraxes.h => src/meraxes.h                                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/meraxes.h                                                                                          \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/meraxes.h                                                           \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/physics/cooling.c => src/physics/cooling.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/physics/cooling.c                                                                                  \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/physics/mergers.c => src/physics/mergers.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/physics/mergers.c                                                                                  \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-CONFLICT (rename/modify): Merge conflict in src/physics/mergers.c                                                   \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Renaming meraxes/src/physics/star_formation.c => src/physics/star_formation.c                                       \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Auto-merging src/physics/star_formation.c                                                                           \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Removing meraxes/README.md                                                                                          \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Removing meraxes/include/.gitignore                                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Removing meraxes/input/Mini_Tiamat.par                                                                              \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Removing meraxes/input/Tiamat_MR.par                                                                                \\\\n\\\"\\\\\\n\"\\\n"\
"-\"+\\\"-Removing meraxes/input/Tiny_Tiamat.par                                                                              \\\\n\\\"\\\\\\n\"\\\n"\
"-\" \\n\"\\\n"\
"-\"-#define GITDIFF_STR \\\"\\\"\\n\"\\\n"\
"-\"diff --git a/src/setup_module.sh b/src/setup_module.sh\\n\"\\\n"\
"-\"index 8d60e62..53d1d5a 100644\\n\"\\\n"\
"-\"--- a/src/setup_module.sh\\n\"\\\n"\
"-\"+++ b/src/setup_module.sh\\n\"\\\n"\
"-\"@@ -1,4 +1,4 @@\\n\"\\\n"\
"-\" #!/bin/bash\\n\"\\\n"\
"-\" module purge\\n\"\\\n"\
"-\"-module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.10.0-openmpi-1.10.2-psm \\n\"\\\n"\
"-\"-for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\":\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\"'\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\n\"\\\n"\
"-\"+module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.8.17-openmpi-1.10.2-psm \\n\"\\\n"\
"-\"+#for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\":\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\"'\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\n\"\\\n"\
"-\"diff --git a/tmp b/tmp\\n\"\\\n"\
" \"deleted file mode 100644\\n\"\\\n"\
"-\"index 38c69fd..0000000\\n\"\\\n"\
"-\"--- a/tmp\\n\"\\\n"\
"+\"index a70501e..0000000\\n\"\\\n"\
"+\"--- a/src/git.h\\n\"\\\n"\
" \"+++ /dev/null\\n\"\\\n"\
"-\"@@ -1,44 +0,0 @@\\n\"\\\n"\
"-\"-Renaming meraxes/src/SConstruct => src/SConstruct                                                                   \\n\"\\\n"\
"-\"-Auto-merging src/SConstruct                                                                                         \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/SConstruct                                                          \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/cleanup.c => src/core/cleanup.c                                                           \\n\"\\\n"\
"-\"-Auto-merging src/core/cleanup.c                                                                                     \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/dracarys.c => src/core/dracarys.c                                                         \\n\"\\\n"\
"-\"-Auto-merging src/core/dracarys.c                                                                                    \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/core/dracarys.c                                                     \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/find_HII_bubbles.c => src/core/find_HII_bubbles.c                                         \\n\"\\\n"\
"-\"-Auto-merging src/core/find_HII_bubbles.c                                                                            \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/core/find_HII_bubbles.c                                             \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/galaxies.c => src/core/galaxies.c                                                         \\n\"\\\n"\
"-\"-Auto-merging src/core/galaxies.c                                                                                    \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/core/galaxies.c                                                     \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/init.c => src/core/init.c                                                                 \\n\"\\\n"\
"-\"-Auto-merging src/core/init.c                                                                                        \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/core/init.c                                                         \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/read_halos.c => src/core/read_halos.c                                                     \\n\"\\\n"\
"-\"-Auto-merging src/core/read_halos.c                                                                                  \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/read_params.c => src/core/read_params.c                                                   \\n\"\\\n"\
"-\"-Auto-merging src/core/read_params.c                                                                                 \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/reionization.c => src/core/reionization.c                                                 \\n\"\\\n"\
"-\"-Auto-merging src/core/reionization.c                                                                                \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/core/reionization.c                                                 \\n\"\\\n"\
"-\"-Renaming meraxes/src/core/save.c => src/core/save.c                                                                 \\n\"\\\n"\
"-\"-Auto-merging src/core/save.c                                                                                        \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/core/save.c                                                         \\n\"\\\n"\
"-\"-Renaming meraxes/src/deps/gstar.py => src/deps/gstar.py                                                             \\n\"\\\n"\
"-\"-Auto-merging src/deps/gstar.py                                                                                      \\n\"\\\n"\
"-\"-Renaming meraxes/src/meraxes.h => src/meraxes.h                                                                     \\n\"\\\n"\
"-\"-Auto-merging src/meraxes.h                                                                                          \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/meraxes.h                                                           \\n\"\\\n"\
"-\"-Renaming meraxes/src/physics/cooling.c => src/physics/cooling.c                                                     \\n\"\\\n"\
"-\"-Auto-merging src/physics/cooling.c                                                                                  \\n\"\\\n"\
"-\"-Renaming meraxes/src/physics/mergers.c => src/physics/mergers.c                                                     \\n\"\\\n"\
"-\"-Auto-merging src/physics/mergers.c                                                                                  \\n\"\\\n"\
"-\"-CONFLICT (rename/modify): Merge conflict in src/physics/mergers.c                                                   \\n\"\\\n"\
"-\"-Renaming meraxes/src/physics/star_formation.c => src/physics/star_formation.c                                       \\n\"\\\n"\
"-\"-Auto-merging src/physics/star_formation.c                                                                           \\n\"\\\n"\
"-\"-Removing meraxes/README.md                                                                                          \\n\"\\\n"\
"-\"-Removing meraxes/include/.gitignore                                                                                 \\n\"\\\n"\
"-\"-Removing meraxes/input/Mini_Tiamat.par                                                                              \\n\"\\\n"\
"-\"-Removing meraxes/input/Tiamat_MR.par                                                                                \\n\"\\\n"\
"-\"-Removing meraxes/input/Tiny_Tiamat.par                                                                              \\n\"\\\n"\
"+\"@@ -1,693 +0,0 @@\\n\"\\\n"\
"+\"-/*\\n\"\\\n"\
"+\"-            * This file has been automatically generated.  Any modifications will be\\n\"\\\n"\
"+\"-            * overwritten during compilation.\\n\"\\\n"\
"+\"-            */\\n\"\\\n"\
"+\"-\\n\"\\\n"\
"+\"-#define GITREF_STR \\\"536fe2ee35656425b9183cf74ae4614c0bfdd60e\\\"\\n\"\\\n"\
"+\"-\\n\"\\\n"\
"+\"-#define GITDIFF_STR \\\"diff --git a/include/meraxes.h b/include/meraxes.h\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index b700754..e3749ba 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/include/meraxes.h\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/include/meraxes.h\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -1,7 +1,6 @@\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" // These defines were added at compilation time\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" // (see SConstruct prepend_user_defines())\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-#define NDEBUG\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-#define USE_MPI 1\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+#define U\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" // -------------------------------------------\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -108,7 +107,12 @@ typedef struct physics_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double EnergyPerSN;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double IMFNormConst;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double RadioModeEff;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double QuasarModeEff;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double BlackHoleGrowthRate;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double EddingtonRatio;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double quasar_mode_scaling;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double quasar_open_angel;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double quasar_fobs;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double ThreshMajorMerger;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double MinMergerStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -120,12 +124,16 @@ typedef struct physics_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double ReionEfficiency;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double ReionNionPhotPerBary;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-  double ReionEscapeFrac;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double ReionEscapeFrac;  \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double ReionEscapeFracBH;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double BlackHoleSeed;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double BlackHoleMassLimitReion;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double ReionTcool;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double Y_He;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double  ReionGammaHaloBias;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double  ReionAlphaUV;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double  ReionAlphaUVBH;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double  ReionRBubbleMin;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double  ReionRBubbleMax;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -147,9 +155,13 @@ typedef struct physics_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double  ReionSMParam_d;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   // Flags\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-  int Flag_RedshiftDepEscFrac;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double RedshiftDepEscFracNorm;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double RedshiftDepEscFracScaling;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double RedshiftDepEscFracBHNorm;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double RedshiftDepEscFracBHScaling;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int Flag_ReionizationModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int Flag_BHFeedback;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  int Flag_BHReion;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int Flag_IRA;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int Flag_FixDiskRadiusOnInfall;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int Flag_FixVmaxOnInfall;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -176,6 +188,8 @@ typedef struct run_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   char             MagBands[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   char             ForestIDFile[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   char             MvirCritFile[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  char             MassRatioModifier[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  char             BaryonFracModifier[STRLEN];\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   physics_params_t physics;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -214,6 +228,7 @@ typedef struct run_params_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int              FlagReadDumpFile;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int              FlagMCMC;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int              Flag_PatchyReion;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  int              Flag_output_grids;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" } run_params_t;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -373,12 +388,20 @@ typedef struct galaxy_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double Mcool;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double StellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double GrossStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double FescWeightedGSM;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double Stellaremissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double MergerSemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double MetalsStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double DiskScaleLength;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double Sfr;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double EjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double MetalsEjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double BlackHoleMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double BHemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double EffectiveBHM;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double BlackHoleAccretedHotMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double BlackHoleAccretedColdMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double BlackHoleAccretingColdMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   // baryonic hostories\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double mwmsa_num;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -444,11 +467,18 @@ typedef struct galaxy_output_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float DiskScaleLength;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float StellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float GrossStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float Stellaremissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float MergerSemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float FescWeightedGSM;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float MetalsStellarMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float Sfr;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float EjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float MetalsEjectedGas;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float BlackHoleMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float BHemissivity;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float EffectiveBHM;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float BlackHoleAccretedHotMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float BlackHoleAccretedColdMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   // misc\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float Rcool;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -457,6 +487,7 @@ typedef struct galaxy_output_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float MergerStartRadius;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float BaryonFracModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float MvirCrit;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  float dt;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   float MergerBurstMass;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   // baryonic histories\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -511,6 +542,17 @@ typedef struct catalog_halo_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   char      padding[8];                //!< Alignment padding\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" } catalog_halo_t;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+typedef struct Modifier\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+{\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float logMmin;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float logMmax;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float mass_mean;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float mass_errl;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float mass_erru;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float ratio;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float ratio_errl;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    float ratio_erru;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+} Modifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" //! Global variables which will will be passed around\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" typedef struct run_globals_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -523,7 +565,10 @@ typedef struct run_globals_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double             *AA;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double             *ZZ;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   double             *LTTime;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  double             *mass_weighted_xHII;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int                *RequestedForestId;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  int                RequestedMassRatioModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  int                RequestedBaryonFracModifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int                *ListOutputSnaps;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   halo_t            **SnapshotHalo;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   fof_group_t       **SnapshotFOFGroup;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -550,6 +595,8 @@ typedef struct run_globals_t {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   int                 NStoreSnapshots;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   bool                SelectForestsSwitch;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  Modifier            *mass_ratio_modifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+  Modifier            *baryon_frac_modifier;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" } run_globals_t;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #ifdef _MAIN\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" run_globals_t run_globals;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -615,7 +662,7 @@ float        apply_pbc_pos(float x);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       accurate_sumf(float *arr, int n);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" int          grid_index(int i, int j, int k, int dim, int type);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         mpi_debug_here(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-int 				 isclosef(float a, float b, float rel_tol, float abs_tol);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+int          isclosef(float a, float b, float rel_tol, float abs_tol);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         printProgress (double percentage);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         check_counts(fof_group_t *fof_group, int NGal, int NFof);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         cn_quote(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -626,14 +673,19 @@ double       calculate_Mvir(double Mvir, int len);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       calculate_Rvir(double Mvir, int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       calculate_Vvir(double Mvir, double Rvir);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       calculate_spin_param(halo_t *halo);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+void         read_mass_ratio_modifiers(int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+void         read_baryon_frac_modifiers(int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+double       interpolate_modifier(Modifier *modifier_data, double logM);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         read_cooling_functions(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       interpolate_cooling_rate(double logTemp, double logZ);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       gas_cooling(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double       calc_metallicity(double total_gas, double metals);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass, double x);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio, int snapshot);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+void         previous_merger_driven_BH_growth(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+double       calculate_BHemissivity(double BlackHoleMass, double accreted_mass);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void         reincorporate_ejected_gas(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" // Magnitude related\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void   init_luminosities(galaxy_t *gal);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -653,6 +705,7 @@ void   assign_slabs(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void   init_reion_grids(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void   filter(fftwf_complex *box, int local_ix_start, int slab_nx, int grid_dim, float R);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+void   set_quasar_fobs(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void   find_HII_bubbles(float redshift);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" double tocf_modifier(galaxy_t *gal, double Mvir);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void   set_ReionEfficiency(void);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/src/SConstruct b/src/SConstruct\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 26fc42d..3ff3daa 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/src/SConstruct\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/src/SConstruct\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -236,6 +236,7 @@ if not GetOption('help'):\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"                 env.AppendUnique(CCFLAGS = ['-I'+dep['inclp']])\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"             if 'libp' in dep and dep['libp'] is not None:\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"                 env.AppendUnique(LIBPATH = [dep['libp']])\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+                env.AppendUnique(RPATH = [dep['libp']])\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"             if 'lib' in dep and dep['lib'] is not None:\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"                 env.AppendUnique(LIBS = [dep['lib']])\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/src/core/modifiers.c b/src/core/modifiers.c\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 9cfc074..71d5025 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/src/core/modifiers.c\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/src/core/modifiers.c\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -4,11 +4,12 @@\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #include <hdf5_hl.h>\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #define NFIELDS 8\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-#define N_START 9\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+#define N_START 0\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #define N_LOGMS 31\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #define DELTA_M 1.0\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #define MIN_LOGM 7.5\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #define MAX_LOGM 11.5\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+#define M_OFFSET 0.5\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" void read_mass_ratio_modifiers(int snapshot){\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"   if (strlen(run_globals.params.MassRatioModifier) == 0){\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/src/core/save.c b/src/core/save.c\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 54e57bb..9d0fd55 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/src/core/save.c\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/src/core/save.c\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -418,7 +418,7 @@ void calc_hdf5_props()\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     h5props->field_h_conv[i]    = \\\\\\\"v/h\\\\\\\";\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     h5props->field_types[i++]   = h5props->array_nhist_f_tid;\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-	// Blackhole or Emissivity related\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    // Blackhole or Emissivity related\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     h5props->dst_field_sizes[i] = sizeof(galout.Stellaremissivity);\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     h5props->field_names[i]     = \\\\\\\"Stellaremissivity\\\\\\\";\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     h5props->field_units[i]     = \\\\\\\"1e60 photons\\\\\\\";\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/src/deps/gstar.py b/src/deps/gstar.py\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 669a023..3361b75 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/src/deps/gstar.py\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/src/deps/gstar.py\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -4,7 +4,7 @@ deps = {}\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" deps['exec'] = {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     'mpicc' : '/usr/local/x86_64/gnu/openmpi-1.10.2-psm/bin/mpicc',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-    'cc' : '/usr/local/gcc-5.1.0/bin/gcc',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+	'cc' : '/usr/local/intel-15.3.0/composer_xe_2015.3.187/bin/intel64/icc',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     'git' : '/usr/bin/git',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" }\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -15,8 +15,8 @@ deps['gsl'] = {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" }\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" deps['hdf5'] = {\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/include',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/lib',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/include',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/lib',\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"     'lib'   : ['hdf5', 'hdf5_hl', 'z'],\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" }\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/src/git.h b/src/git.h\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 3d521ae..0b21f2a 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/src/git.h\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/src/git.h\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -3,6 +3,349 @@\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"             * overwritten during compilation.\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"             */\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-#define GITREF_STR \\\\\\\"7cd7348aab75d26b87669795b4080bf8b0f6cf0b\\\\\\\"\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+#define GITREF_STR \\\\\\\"536fe2ee35656425b9183cf74ae4614c0bfdd60e\\\\\\\"\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+#define GITDIFF_STR \\\\\\\"diff --git a/include/meraxes.h b/include/meraxes.h\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index b700754..e3749ba 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/include/meraxes.h\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ b/include/meraxes.h\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -1,7 +1,6 @@\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" // These defines were added at compilation time\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" // (see SConstruct prepend_user_defines())\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-#define NDEBUG\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-#define USE_MPI 1\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+#define U\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" // -------------------------------------------\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -108,7 +107,12 @@ typedef struct physics_params_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double EnergyPerSN;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double IMFNormConst;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double RadioModeEff;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double QuasarModeEff;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double BlackHoleGrowthRate;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double EddingtonRatio;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double quasar_mode_scaling;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double quasar_open_angel;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double quasar_fobs;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double ThreshMajorMerger;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double MinMergerStellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -120,12 +124,16 @@ typedef struct physics_params_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double ReionEfficiency;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double ReionNionPhotPerBary;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-  double ReionEscapeFrac;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double ReionEscapeFrac;  \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double ReionEscapeFracBH;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double BlackHoleSeed;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double BlackHoleMassLimitReion;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double ReionTcool;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double Y_He;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double  ReionGammaHaloBias;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double  ReionAlphaUV;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double  ReionAlphaUVBH;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double  ReionRBubbleMin;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double  ReionRBubbleMax;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -147,9 +155,13 @@ typedef struct physics_params_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double  ReionSMParam_d;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   // Flags\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-  int Flag_RedshiftDepEscFrac;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double RedshiftDepEscFracNorm;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double RedshiftDepEscFracScaling;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double RedshiftDepEscFracBHNorm;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double RedshiftDepEscFracBHScaling;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int Flag_ReionizationModifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int Flag_BHFeedback;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  int Flag_BHReion;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int Flag_IRA;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int Flag_FixDiskRadiusOnInfall;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int Flag_FixVmaxOnInfall;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -176,6 +188,8 @@ typedef struct run_params_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   char             MagBands[STRLEN];\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   char             ForestIDFile[STRLEN];\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   char             MvirCritFile[STRLEN];\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  char             MassRatioModifier[STRLEN];\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  char             BaryonFracModifier[STRLEN];\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   physics_params_t physics;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -214,6 +228,7 @@ typedef struct run_params_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int              FlagReadDumpFile;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int              FlagMCMC;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int              Flag_PatchyReion;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  int              Flag_output_grids;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" } run_params_t;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -373,12 +388,20 @@ typedef struct galaxy_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double Mcool;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double StellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double GrossStellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double FescWeightedGSM;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double Stellaremissivity;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double MergerSemissivity;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double MetalsStellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double DiskScaleLength;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double Sfr;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double EjectedGas;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double MetalsEjectedGas;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double BlackHoleMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double BHemissivity;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double EffectiveBHM;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double BlackHoleAccretedHotMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double BlackHoleAccretedColdMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double BlackHoleAccretingColdMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   // baryonic hostories\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double mwmsa_num;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -444,11 +467,18 @@ typedef struct galaxy_output_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float DiskScaleLength;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float StellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float GrossStellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float Stellaremissivity;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float MergerSemissivity;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float FescWeightedGSM;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float MetalsStellarMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float Sfr;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float EjectedGas;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float MetalsEjectedGas;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float BlackHoleMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float BHemissivity;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float EffectiveBHM;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float BlackHoleAccretedHotMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float BlackHoleAccretedColdMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   // misc\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float Rcool;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -457,6 +487,7 @@ typedef struct galaxy_output_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float MergerStartRadius;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float BaryonFracModifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float MvirCrit;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  float dt;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   float MergerBurstMass;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   // baryonic histories\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -511,6 +542,17 @@ typedef struct catalog_halo_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   char      padding[8];                //!< Alignment padding\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" } catalog_halo_t;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+typedef struct Modifier\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+{\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float logMmin;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float logMmax;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float mass_mean;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float mass_errl;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float mass_erru;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float ratio;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float ratio_errl;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    float ratio_erru;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+} Modifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" //! Global variables which will will be passed around\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" typedef struct run_globals_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -523,7 +565,10 @@ typedef struct run_globals_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double             *AA;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double             *ZZ;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   double             *LTTime;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  double             *mass_weighted_xHII;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int                *RequestedForestId;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  int                RequestedMassRatioModifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  int                RequestedBaryonFracModifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int                *ListOutputSnaps;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   halo_t            **SnapshotHalo;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   fof_group_t       **SnapshotFOFGroup;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -550,6 +595,8 @@ typedef struct run_globals_t {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   int                 NStoreSnapshots;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   bool                SelectForestsSwitch;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  Modifier            *mass_ratio_modifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+  Modifier            *baryon_frac_modifier;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" } run_globals_t;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #ifdef _MAIN\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" run_globals_t run_globals;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -615,7 +662,7 @@ float        apply_pbc_pos(float x);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       accurate_sumf(float *arr, int n);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" int          grid_index(int i, int j, int k, int dim, int type);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         mpi_debug_here(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-int 				 isclosef(float a, float b, float rel_tol, float abs_tol);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+int          isclosef(float a, float b, float rel_tol, float abs_tol);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         printProgress (double percentage);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         check_counts(fof_group_t *fof_group, int NGal, int NFof);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         cn_quote(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -626,14 +673,19 @@ double       calculate_Mvir(double Mvir, int len);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       calculate_Rvir(double Mvir, int snapshot);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       calculate_Vvir(double Mvir, double Rvir);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       calculate_spin_param(halo_t *halo);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+void         read_mass_ratio_modifiers(int snapshot);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+void         read_baryon_frac_modifiers(int snapshot);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+double       interpolate_modifier(Modifier *modifier_data, double logM);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         read_cooling_functions(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       interpolate_cooling_rate(double logTemp, double logZ);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       gas_cooling(galaxy_t *gal);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         cool_gas_onto_galaxy(galaxy_t *gal, double cooling_mass);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double       calc_metallicity(double total_gas, double metals);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass, double x);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio, int snapshot);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+void         previous_merger_driven_BH_growth(galaxy_t *gal);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+double       calculate_BHemissivity(double BlackHoleMass, double accreted_mass);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void         reincorporate_ejected_gas(galaxy_t *gal);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-double       radio_mode_BH_heating(galaxy_t *gal, double cooling_mass);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-void         merger_driven_BH_growth(galaxy_t *gal, double merger_ratio);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" // Magnitude related\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void   init_luminosities(galaxy_t *gal);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -653,6 +705,7 @@ void   assign_slabs(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void   init_reion_grids(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void   filter(fftwf_complex *box, int local_ix_start, int slab_nx, int grid_dim, float R);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+void   set_quasar_fobs(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void   find_HII_bubbles(float redshift);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" double tocf_modifier(galaxy_t *gal, double Mvir);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void   set_ReionEfficiency(void);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/src/SConstruct b/src/SConstruct\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 26fc42d..3ff3daa 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/src/SConstruct\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ b/src/SConstruct\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -236,6 +236,7 @@ if not GetOption('help'):\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"                 env.AppendUnique(CCFLAGS = ['-I'+dep['inclp']])\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"             if 'libp' in dep and dep['libp'] is not None:\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"                 env.AppendUnique(LIBPATH = [dep['libp']])\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+                env.AppendUnique(RPATH = [dep['libp']])\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"             if 'lib' in dep and dep['lib'] is not None:\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"                 env.AppendUnique(LIBS = [dep['lib']])\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/src/core/modifiers.c b/src/core/modifiers.c\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 9cfc074..71d5025 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/src/core/modifiers.c\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ b/src/core/modifiers.c\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -4,11 +4,12 @@\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #include <hdf5_hl.h>\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #define NFIELDS 8\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-#define N_START 9\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+#define N_START 0\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #define N_LOGMS 31\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #define DELTA_M 1.0\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #define MIN_LOGM 7.5\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #define MAX_LOGM 11.5\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+#define M_OFFSET 0.5\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" void read_mass_ratio_modifiers(int snapshot){\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"   if (strlen(run_globals.params.MassRatioModifier) == 0){\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/src/core/save.c b/src/core/save.c\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 54e57bb..9d0fd55 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/src/core/save.c\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ b/src/core/save.c\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -418,7 +418,7 @@ void calc_hdf5_props()\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     h5props->field_h_conv[i]    = \\\\\\\\\\\\\\\"v/h\\\\\\\\\\\\\\\";\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     h5props->field_types[i++]   = h5props->array_nhist_f_tid;\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-	// Blackhole or Emissivity related\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    // Blackhole or Emissivity related\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     h5props->dst_field_sizes[i] = sizeof(galout.Stellaremissivity);\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     h5props->field_names[i]     = \\\\\\\\\\\\\\\"Stellaremissivity\\\\\\\\\\\\\\\";\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     h5props->field_units[i]     = \\\\\\\\\\\\\\\"1e60 photons\\\\\\\\\\\\\\\";\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/src/deps/gstar.py b/src/deps/gstar.py\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 669a023..3361b75 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/src/deps/gstar.py\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ b/src/deps/gstar.py\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -4,7 +4,7 @@ deps = {}\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" deps['exec'] = {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     'mpicc' : '/usr/local/x86_64/gnu/openmpi-1.10.2-psm/bin/mpicc',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-    'cc' : '/usr/local/gcc-5.1.0/bin/gcc',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+	'cc' : '/usr/local/intel-15.3.0/composer_xe_2015.3.187/bin/intel64/icc',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     'git' : '/usr/bin/git',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" }\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -15,8 +15,8 @@ deps['gsl'] = {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" }\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" deps['hdf5'] = {\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/include',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/lib',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/include',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.8.17-openmpi-1.10.2-psm/lib',\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"     'lib'   : ['hdf5', 'hdf5_hl', 'z'],\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" }\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/src/git.h b/src/git.h\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"deleted file mode 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 3d521ae..0000000\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/src/git.h\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ /dev/null\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -1,8 +0,0 @@\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-/*\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-            * This file has been automatically generated.  Any modifications will be\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-            * overwritten during compilation.\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-            */\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-#define GITREF_STR \\\\\\\\\\\\\\\"7cd7348aab75d26b87669795b4080bf8b0f6cf0b\\\\\\\\\\\\\\\"\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-#define GITDIFF_STR \\\\\\\\\\\\\\\"\\\\\\\\\\\\\\\"\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/src/setup_module.sh b/src/setup_module.sh\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 8d60e62..53d1d5a 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/src/setup_module.sh\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ b/src/setup_module.sh\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -1,4 +1,4 @@\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" #!/bin/bash\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\" module purge\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.10.0-openmpi-1.10.2-psm \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\\\\\\\\\\\\\":\\\\\\\\\\\\\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\\\\\\\\\\\\\"'\\\\\\\\\\\\\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.8.17-openmpi-1.10.2-psm \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+#for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\\\\\\\\\\\\\":\\\\\\\\\\\\\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\\\\\\\\\\\\\"'\\\\\\\\\\\\\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"diff --git a/tmp b/tmp\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"deleted file mode 100644\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"index 38c69fd..0000000\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"--- a/tmp\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"+++ /dev/null\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"@@ -1,44 +0,0 @@\\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/SConstruct => src/SConstruct                                                                   \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/SConstruct                                                                                         \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/SConstruct                                                          \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/cleanup.c => src/core/cleanup.c                                                           \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/cleanup.c                                                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/dracarys.c => src/core/dracarys.c                                                         \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/dracarys.c                                                                                    \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/core/dracarys.c                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/find_HII_bubbles.c => src/core/find_HII_bubbles.c                                         \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/find_HII_bubbles.c                                                                            \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/core/find_HII_bubbles.c                                             \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/galaxies.c => src/core/galaxies.c                                                         \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/galaxies.c                                                                                    \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/core/galaxies.c                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/init.c => src/core/init.c                                                                 \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/init.c                                                                                        \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/core/init.c                                                         \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/read_halos.c => src/core/read_halos.c                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/read_halos.c                                                                                  \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/read_params.c => src/core/read_params.c                                                   \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/read_params.c                                                                                 \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/reionization.c => src/core/reionization.c                                                 \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/reionization.c                                                                                \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/core/reionization.c                                                 \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/core/save.c => src/core/save.c                                                                 \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/core/save.c                                                                                        \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/core/save.c                                                         \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/deps/gstar.py => src/deps/gstar.py                                                             \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/deps/gstar.py                                                                                      \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/meraxes.h => src/meraxes.h                                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/meraxes.h                                                                                          \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/meraxes.h                                                           \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/physics/cooling.c => src/physics/cooling.c                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/physics/cooling.c                                                                                  \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/physics/mergers.c => src/physics/mergers.c                                                     \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/physics/mergers.c                                                                                  \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-CONFLICT (rename/modify): Merge conflict in src/physics/mergers.c                                                   \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Renaming meraxes/src/physics/star_formation.c => src/physics/star_formation.c                                       \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Auto-merging src/physics/star_formation.c                                                                           \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Removing meraxes/README.md                                                                                          \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Removing meraxes/include/.gitignore                                                                                 \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Removing meraxes/input/Mini_Tiamat.par                                                                              \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Removing meraxes/input/Tiamat_MR.par                                                                                \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+\\\\\\\"-Removing meraxes/input/Tiny_Tiamat.par                                                                              \\\\\\\\n\\\\\\\"\\\\\\\\\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-#define GITDIFF_STR \\\\\\\"\\\\\\\"\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/src/setup_module.sh b/src/setup_module.sh\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 8d60e62..53d1d5a 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/src/setup_module.sh\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ b/src/setup_module.sh\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -1,4 +1,4 @@\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" #!/bin/bash\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\" module purge\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.10.0-openmpi-1.10.2-psm \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\\\\\":\\\\\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\\\\\"'\\\\\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.8.17-openmpi-1.10.2-psm \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+#for libp in `grep 'libp' ~/bitbucket/new_meraxes/src/deps/gstar.py| awk 'BEGIN{FS=\\\\\\\":\\\\\\\"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d \\\\\\\"'\\\\\\\"|tr -d ','`:$LD_LIBRARY_PATH; done\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"diff --git a/tmp b/tmp\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"deleted file mode 100644\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"index 38c69fd..0000000\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"--- a/tmp\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"+++ /dev/null\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"@@ -1,44 +0,0 @@\\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/SConstruct => src/SConstruct                                                                   \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/SConstruct                                                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/SConstruct                                                          \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/cleanup.c => src/core/cleanup.c                                                           \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/cleanup.c                                                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/dracarys.c => src/core/dracarys.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/dracarys.c                                                                                    \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/core/dracarys.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/find_HII_bubbles.c => src/core/find_HII_bubbles.c                                         \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/find_HII_bubbles.c                                                                            \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/core/find_HII_bubbles.c                                             \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/galaxies.c => src/core/galaxies.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/galaxies.c                                                                                    \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/core/galaxies.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/init.c => src/core/init.c                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/init.c                                                                                        \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/core/init.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/read_halos.c => src/core/read_halos.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/read_halos.c                                                                                  \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/read_params.c => src/core/read_params.c                                                   \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/read_params.c                                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/reionization.c => src/core/reionization.c                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/reionization.c                                                                                \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/core/reionization.c                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/core/save.c => src/core/save.c                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/core/save.c                                                                                        \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/core/save.c                                                         \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/deps/gstar.py => src/deps/gstar.py                                                             \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/deps/gstar.py                                                                                      \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/meraxes.h => src/meraxes.h                                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/meraxes.h                                                                                          \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/meraxes.h                                                           \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/physics/cooling.c => src/physics/cooling.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/physics/cooling.c                                                                                  \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/physics/mergers.c => src/physics/mergers.c                                                     \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/physics/mergers.c                                                                                  \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-CONFLICT (rename/modify): Merge conflict in src/physics/mergers.c                                                   \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Renaming meraxes/src/physics/star_formation.c => src/physics/star_formation.c                                       \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Auto-merging src/physics/star_formation.c                                                                           \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Removing meraxes/README.md                                                                                          \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Removing meraxes/include/.gitignore                                                                                 \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Removing meraxes/input/Mini_Tiamat.par                                                                              \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Removing meraxes/input/Tiamat_MR.par                                                                                \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\\"-Removing meraxes/input/Tiny_Tiamat.par                                                                              \\\\n\\\"\\\\\\n\"\\\n"\
"+\"-\\n\"\\\n"\
" \n"\

