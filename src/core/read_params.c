#include "meraxes.h"
#include "parse_paramfile.h"
#include <assert.h>

static void check_problem_params(run_params_t* run_params)
{
  if (run_params->NSteps != 1) {
    mlog_error("The current version of the code only works if NSteps = 1. Sorry! Exiting...");
    ABORT(EXIT_FAILURE);
  }

  if (run_params->Flag_SeparateQSOXrays != 0) {
    mlog_error("The current version of the code only works if Flag_SeparateQSOXrays = 0. Sorry! Exiting...");
    ABORT(EXIT_FAILURE);
  }

  if (strlen(run_globals.params.ForestIDFile) != 0) {
    mlog("*** YOU HAVE PROVIDED A REQUESTED FORESTID FILE. THIS FEATURE HAS NOT BE WELL TESTED. YMMV! ***", MLOG_MESG);
  }

#ifdef USE_CUDA
  if (run_params->Flag_IncludeSpinTemp != 0) {
    mlog_error(
      "Spin temperature features are not currently available in the GPU version of find_HII_bubbles!  Exiting...");
    ABORT(EXIT_FAILURE);
  }
#endif
}

static void store_params(entry_t entry[123],
                         int n_entries,
                         char** params_tag,
                         int n_param,
                         int used_tag[123],
                         const int* params_type,
                         void** params_addr)
{
  int level = 0;
  char prefix[16] = "\0";
  char key[STRLEN + 64];

  for (int i_entry = 0; i_entry < n_entries; i_entry++) {
    // DEBUG
    // mlog("Checking %s", MLOG_MESG, entry[i_entry].key);

    // reset prefix if we have descended an indentation level
    if (entry[i_entry].level < level)
      *prefix = '\0';

    strncpy(key, prefix, STRLEN);
    strncat(key, entry[i_entry].key, STRLEN);
    level = entry[i_entry].level;

    // DEBUG
    // mlog("level = %d :: prefix = %s", MLOG_MESG, level, prefix);

    int tag_index = -1;
    for (int ii = 0; ii < n_param; ii++)
      if (strcmp(key, params_tag[ii]) == 0) {
        tag_index = ii;
        break;
      }

    if (tag_index < 0) {
      mlog("<WARNING> %s is an unrecognised parameter (prefix='%s').", MLOG_MESG, entry[i_entry].key, prefix);
      ABORT(EXIT_FAILURE);
    }

    if (used_tag[tag_index] == 1)
      continue;

    switch (params_type[tag_index]) {
      case PARAM_TYPE_DOUBLE:
        *((double*)params_addr[tag_index]) = atof(entry[i_entry].value);
        break;

      case PARAM_TYPE_FLOAT:
        *((float*)params_addr[tag_index]) = (float)atof(entry[i_entry].value);
        break;

      case PARAM_TYPE_STRING:
        strncpy(params_addr[tag_index], entry[i_entry].value, STRLEN);
        break;

      case PARAM_TYPE_INT:
        *((int*)params_addr[tag_index]) = atoi(entry[i_entry].value);
        break;

      case PARAM_TYPE_LONGLONG:
        *((long long*)params_addr[tag_index]) = atoll(entry[i_entry].value);
        break;

      default:
        mlog_error("Unknown param type.");
        break;
    }
    used_tag[tag_index] = 1;
  }
}

void read_parameter_file(char* fname, int mode)
{
  // mode = 0 : for single runs
  // mode = 1 : for interactive runs (do not remalloc arrays)

  run_params_t* run_params = &(run_globals.params);

  if (run_globals.mpi_rank == 0) {
    int ii;
    int used_tag[PARAM_MAX_ENTRIES], required_tag[PARAM_MAX_ENTRIES];
    entry_t entry[PARAM_MAX_ENTRIES];

    int n_entries;
    hdf5_output_t* hdf5props = &(run_globals.hdf5props);
    int* params_type = hdf5props->params_type;
    void** params_addr = hdf5props->params_addr;
    char** params_tag = hdf5props->params_tag;
    int n_param = hdf5props->params_count;

    mlog("\nreading parameter file:\n", MLOG_MESG);

    // even if this is an interactive run, we have to init required_tag
    // (used_tag is initialised later)
    for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
      required_tag[ii] = 0;

    // malloc global arrays and init param properties
    if (mode == 0) {
      const int tag_length = 128;
      hdf5props->params_tag = malloc(sizeof(char*) * PARAM_MAX_ENTRIES);
      for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
        hdf5props->params_tag[ii] = malloc(sizeof(char) * tag_length);
      params_tag = hdf5props->params_tag;
      hdf5props->params_type = malloc(sizeof(int) * PARAM_MAX_ENTRIES);
      params_type = hdf5props->params_type;
      hdf5props->params_addr = malloc(sizeof(void*) * PARAM_MAX_ENTRIES);
      params_addr = hdf5props->params_addr;

      // Initialise values and arrays
      n_param = 0;
      for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++) {
        required_tag[ii] = 0;
        hdf5props->params_type[ii] = PARAM_TYPE_UNUSED;
      }

      strncpy(params_tag[n_param], "DefaultsFile", tag_length);
      params_addr[n_param] = run_params->DefaultsFile;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "OutputDir", tag_length);
      params_addr[n_param] = run_params->OutputDir;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "OutputSnapshots", tag_length);
      params_addr[n_param] = run_params->OutputSnapsString;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "PhotometricTablesDir", tag_length);
      params_addr[n_param] = run_params->PhotometricTablesDir;
#ifndef CALC_MAGS
      required_tag[n_param] = 0;
#else
      required_tag[n_param] = 1;
#endif
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "TargetSnaps");
      params_addr[n_param] = run_params->TargetSnaps;
#ifndef CALC_MAGS
      required_tag[n_param] = 0;
#else
      required_tag[n_param] = 1;
#endif
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "BetaBands", tag_length);
      params_addr[n_param] = run_params->BetaBands;
#ifndef CALC_MAGS
      required_tag[n_param] = 0;
#else
      required_tag[n_param] = 1;
#endif
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "RestBands");
      params_addr[n_param] = run_params->RestBands;
#ifndef CALC_MAGS
      required_tag[n_param] = 0;
#else
      required_tag[n_param] = 1;
#endif
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "BirthCloudLifetime");
      params_addr[n_param] = &(run_params->BirthCloudLifetime);
#ifndef CALC_MAGS
      required_tag[n_param] = 0;
#else
      required_tag[n_param] = 1;
#endif
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "TablesForXHeatingDir", tag_length);
      params_addr[n_param] = run_params->TablesForXHeatingDir;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "CoolingFuncsDir");
      params_addr[n_param] = run_params->CoolingFuncsDir;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strcpy(params_tag[n_param], "StellarFeedbackDir");
      params_addr[n_param] = run_params->StellarFeedbackDir;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "FileNameGalaxies", tag_length);
      params_addr[n_param] = run_params->FileNameGalaxies;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "SimName", tag_length);
      params_addr[n_param] = run_params->SimName;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "TreesID", tag_length);
      params_addr[n_param] = &(run_params->TreesID);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "SimulationDir", tag_length);
      params_addr[n_param] = run_params->SimulationDir;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "CatalogFilePrefix", tag_length);
      params_addr[n_param] = run_params->CatalogFilePrefix;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;

      strncpy(params_tag[n_param], "NSteps", tag_length);
      params_addr[n_param] = &(run_params->NSteps);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "BoxSize", tag_length);
      params_addr[n_param] = &(run_params->BoxSize);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "VolumeFactor", tag_length);
      params_addr[n_param] = &(run_params->VolumeFactor);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "RandomSeed", tag_length);
      params_addr[n_param] = &(run_params->RandomSeed);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ForestIDFile", tag_length);
      params_addr[n_param] = &(run_params->ForestIDFile);
      required_tag[n_param] = 0;
      params_type[n_param++] = PARAM_TYPE_STRING;
      *(run_params->ForestIDFile) = '\0';

      strncpy(params_tag[n_param], "MvirCritFile", tag_length);
      params_addr[n_param] = &(run_params->MvirCritFile);
      required_tag[n_param] = 0;
      params_type[n_param++] = PARAM_TYPE_STRING;
      *(run_params->MvirCritFile) = '\0';

      strncpy(params_tag[n_param], "MvirCritMCFile", tag_length);
      params_addr[n_param] = &(run_params->MvirCritMCFile);
      required_tag[n_param] = 0;
      params_type[n_param++] = PARAM_TYPE_STRING;
      *(run_params->MvirCritMCFile) = '\0';

      strncpy(params_tag[n_param], "MassRatioModifier", tag_length);
      params_addr[n_param] = &(run_params->MassRatioModifier);
      required_tag[n_param] = 0;
      params_type[n_param++] = PARAM_TYPE_STRING;
      *(run_params->MassRatioModifier) = '\0';

      strncpy(params_tag[n_param], "BaryonFracModifier", tag_length);
      params_addr[n_param] = &(run_params->BaryonFracModifier);
      required_tag[n_param] = 0;
      params_type[n_param++] = PARAM_TYPE_STRING;
      *(run_params->BaryonFracModifier) = '\0';

      strncpy(params_tag[n_param], "FFTW3WisdomDir", tag_length);
      params_addr[n_param] = &(run_params->FFTW3WisdomDir);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_STRING;
      *(run_params->FFTW3WisdomDir) = '\0';

      strncpy(params_tag[n_param], "UnitVelocity_in_cm_per_s", tag_length);
      params_addr[n_param] = &(run_globals.units.UnitVelocity_in_cm_per_s);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "UnitLength_in_cm", tag_length);
      params_addr[n_param] = &(run_globals.units.UnitLength_in_cm);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "UnitMass_in_g", tag_length);
      params_addr[n_param] = &(run_globals.units.UnitMass_in_g);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "Hubble_h", tag_length);
      params_addr[n_param] = &(run_params->Hubble_h);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "BaryonFrac", tag_length);
      params_addr[n_param] = &(run_params->BaryonFrac);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "OmegaM", tag_length);
      params_addr[n_param] = &(run_params->OmegaM);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "OmegaK", tag_length);
      params_addr[n_param] = &(run_params->OmegaK);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "OmegaR", tag_length);
      params_addr[n_param] = &(run_params->OmegaR);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "OmegaLambda", tag_length);
      params_addr[n_param] = &(run_params->OmegaLambda);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "Sigma8", tag_length);
      params_addr[n_param] = &(run_params->Sigma8);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "wLambda", tag_length);
      params_addr[n_param] = &(run_params->wLambda);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SpectralIndex", tag_length);
      params_addr[n_param] = &(run_params->SpectralIndex);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "PartMass", tag_length);
      params_addr[n_param] = &(run_params->PartMass);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "NPart", tag_length);
      params_addr[n_param] = &(run_params->NPart);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_LONGLONG;

      strncpy(params_tag[n_param], "MergerTimeFactor", tag_length);
      params_addr[n_param] = &(run_params->physics.MergerTimeFactor);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "FlagSubhaloVirialProps", tag_length);
      params_addr[n_param] = &(run_params->FlagSubhaloVirialProps);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "FlagInteractive", tag_length);
      params_addr[n_param] = &(run_params->FlagInteractive);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "FlagMCMC", tag_length);
      params_addr[n_param] = &(run_params->FlagMCMC);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "FlagIgnoreProgIndex", tag_length);
      params_addr[n_param] = &(run_params->FlagIgnoreProgIndex);
      required_tag[n_param] = 0;
      params_type[n_param++] = PARAM_TYPE_INT;
      run_params->FlagIgnoreProgIndex = 0;

      // Physics params

      strncpy(params_tag[n_param], "EscapeFracDependency", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracDependency;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "SfDiskVelOpt", tag_length);
      params_addr[n_param] = &(run_params->physics).SfDiskVelOpt;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "SfPrescription", tag_length);
      params_addr[n_param] = &(run_params->physics).SfPrescription;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_ReionizationModifier", tag_length);
      params_addr[n_param] = &(run_params->physics).Flag_ReionizationModifier;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_BHFeedback", tag_length);
      params_addr[n_param] = &(run_params->physics).Flag_BHFeedback;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_IRA", tag_length);
      params_addr[n_param] = &(run_params->physics).Flag_IRA;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_FixDiskRadiusOnInfall", tag_length);
      params_addr[n_param] = &(run_params->physics).Flag_FixDiskRadiusOnInfall;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_FixVmaxOnInfall", tag_length);
      params_addr[n_param] = &(run_params->physics).Flag_FixVmaxOnInfall;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_ReheatToFOFGroupTemp", tag_length);
      params_addr[n_param] = &(run_params->physics).Flag_ReheatToFOFGroupTemp;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "SfEfficiency", tag_length);
      params_addr[n_param] = &(run_params->physics).SfEfficiency;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SfEfficiencyScaling", tag_length);
      params_addr[n_param] = &(run_params->physics).SfEfficiencyScaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SfCriticalSDNorm", tag_length);
      params_addr[n_param] = &(run_params->physics).SfCriticalSDNorm;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SfRecycleFraction", tag_length);
      params_addr[n_param] = &(run_params->physics).SfRecycleFraction;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnModel", tag_length);
      params_addr[n_param] = &(run_params->physics).SnModel;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "SnEjectionRedshiftDep", tag_length);
      params_addr[n_param] = &(run_params->physics).SnEjectionRedshiftDep;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnEjectionEff", tag_length);
      params_addr[n_param] = &(run_params->physics).SnEjectionEff;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnEjectionScaling", tag_length);
      params_addr[n_param] = &(run_params->physics).SnEjectionScaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnEjectionScaling2", tag_length);
      params_addr[n_param] = &(run_params->physics).SnEjectionScaling2;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnEjectionNorm", tag_length);
      params_addr[n_param] = &(run_params->physics).SnEjectionNorm;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnReheatRedshiftDep", tag_length);
      params_addr[n_param] = &(run_params->physics).SnReheatRedshiftDep;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnReheatEff", tag_length);
      params_addr[n_param] = &(run_params->physics).SnReheatEff;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnReheatLimit", tag_length);
      params_addr[n_param] = &(run_params->physics).SnReheatLimit;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnReheatScaling", tag_length);
      params_addr[n_param] = &(run_params->physics).SnReheatScaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnReheatScaling2", tag_length);
      params_addr[n_param] = &(run_params->physics).SnReheatScaling2;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "SnReheatNorm", tag_length);
      params_addr[n_param] = &(run_params->physics).SnReheatNorm;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReincorporationModel", tag_length);
      params_addr[n_param] = &(run_params->physics).ReincorporationModel;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ReincorporationEff", tag_length);
      params_addr[n_param] = &(run_params->physics).ReincorporationEff;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "MaxCoolingMassFactor", tag_length);
      params_addr[n_param] = &(run_params->physics).MaxCoolingMassFactor;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "Yield", tag_length);
      params_addr[n_param] = &(run_params->physics).Yield;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ThreshMajorMerger", tag_length);
      params_addr[n_param] = &((run_params->physics).ThreshMajorMerger);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "MinMergerStellarMass", tag_length);
      params_addr[n_param] = &((run_params->physics).MinMergerStellarMass);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "MinMergerRatioForBurst", tag_length);
      params_addr[n_param] = &((run_params->physics).MinMergerRatioForBurst);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "MergerBurstScaling", tag_length);
      params_addr[n_param] = &((run_params->physics).MergerBurstScaling);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "MergerBurstFactor", tag_length);
      params_addr[n_param] = &((run_params->physics).MergerBurstFactor);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "RadioModeEff", tag_length);
      params_addr[n_param] = &(run_params->physics).RadioModeEff;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "QuasarModeEff", tag_length);
      params_addr[n_param] = &(run_params->physics).QuasarModeEff;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "BlackHoleSeed", tag_length);
      params_addr[n_param] = &(run_params->physics).BlackHoleSeed;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "BlackHoleGrowthRate", tag_length);
      params_addr[n_param] = &(run_params->physics).BlackHoleGrowthRate;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EddingtonRatio", tag_length);
      params_addr[n_param] = &(run_params->physics).EddingtonRatio;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "quasar_mode_scaling", tag_length);
      params_addr[n_param] = &(run_params->physics).quasar_mode_scaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "quasar_open_angle", tag_length);
      params_addr[n_param] = &(run_params->physics).quasar_open_angle;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSobacchi_Zre", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSobacchi_Zre;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSobacchi_DeltaZre", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSobacchi_DeltaZre;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSobacchi_DeltaZsc", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSobacchi_DeltaZsc;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSobacchi_T0", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSobacchi_T0;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionTcool", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionTcool;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionNionPhotPerBary", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionNionPhotPerBary;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "BlackHoleMassLimitReion", tag_length);
      params_addr[n_param] = &(run_params->physics).BlackHoleMassLimitReion;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionGnedin_z0", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionGnedin_z0;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionGnedin_zr", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionGnedin_zr;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "Flag_PatchyReion", tag_length);
      params_addr[n_param] = &(run_params->Flag_PatchyReion);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_IncludeLymanWerner", tag_length);
      params_addr[n_param] = &(run_params->Flag_IncludeLymanWerner);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_IncludeSpinTemp", tag_length);
      params_addr[n_param] = &(run_params->Flag_IncludeSpinTemp);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_IncludeRecombinations", tag_length);
      params_addr[n_param] = &(run_params->Flag_IncludeRecombinations);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_Compute21cmBrightTemp", tag_length);
      params_addr[n_param] = &(run_params->Flag_Compute21cmBrightTemp);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "TsNumFilterSteps", tag_length);
      params_addr[n_param] = &(run_params->TsNumFilterSteps);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_ComputePS", tag_length);
      params_addr[n_param] = &(run_params->Flag_ComputePS);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_IncludePecVelsFor21cm", tag_length);
      params_addr[n_param] = &(run_params->Flag_IncludePecVelsFor21cm);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "TsVelocityComponent", tag_length);
      params_addr[n_param] = &(run_params->TsVelocityComponent);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_ConstructLightcone", tag_length);
      params_addr[n_param] = &(run_params->Flag_ConstructLightcone);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ReionSfrTimescale", tag_length);
      params_addr[n_param] = &(run_params->ReionSfrTimescale);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EndRedshiftLightcone", tag_length);
      params_addr[n_param] = &(run_params->EndRedshiftLightcone);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      //            strncpy(params_tag[n_param], "CurrentLCPos", tag_length);
      //            params_addr[n_param] = &(run_params->CurrentLCPos);
      //            required_tag[n_param] = 1;
      //            params_type[n_param++] = PARAM_TYPE_LONGLONG;

      strncpy(params_tag[n_param], "Flag_SeparateQSOXrays", tag_length);
      params_addr[n_param] = &(run_params->Flag_SeparateQSOXrays);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_OutputGrids", tag_length);
      params_addr[n_param] = &(run_params->Flag_OutputGrids);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Flag_OutputGridsPostReion", tag_length);
      params_addr[n_param] = &(run_params->Flag_OutputGridsPostReion);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ReionGridDim", tag_length);
      params_addr[n_param] = &(run_params->ReionGridDim);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ReionRBubbleMin", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionRBubbleMin;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionRBubbleMax", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionRBubbleMax;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionRBubbleMaxRecomb", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionRBubbleMaxRecomb;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionDeltaRFactor", tag_length);
      params_addr[n_param] = &(run_params->ReionDeltaRFactor);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionGammaHaloBias", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionGammaHaloBias;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EscapeFracNorm", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracNorm;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EscapeFracRedshiftOffset", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracRedshiftOffset;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EscapeFracRedshiftScaling", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracRedshiftScaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EscapeFracPropScaling", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracPropScaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EscapeFracBHNorm", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracBHNorm;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "EscapeFracBHScaling", tag_length);
      params_addr[n_param] = &(run_params->physics).EscapeFracBHScaling;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSMParam_m0", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSMParam_m0;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSMParam_a", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSMParam_a;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSMParam_b", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSMParam_b;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSMParam_c", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSMParam_c;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionSMParam_d", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionSMParam_d;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionUVBFlag", tag_length);
      params_addr[n_param] = &(run_params->ReionUVBFlag);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ReionFilterType", tag_length);
      params_addr[n_param] = &(run_params->ReionFilterType);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "TsHeatingFilterType", tag_length);
      params_addr[n_param] = &(run_params->TsHeatingFilterType);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "ReionPowerSpecDeltaK", tag_length);
      params_addr[n_param] = &(run_params->ReionPowerSpecDeltaK);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionAlphaUV", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionAlphaUV;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionAlphaUVBH", tag_length);
      params_addr[n_param] = &(run_params->physics).ReionAlphaUVBH;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strncpy(params_tag[n_param], "ReionRtoMFilterType", tag_length);
      params_addr[n_param] = &(run_params->ReionRtoMFilterType);
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_INT;

      strncpy(params_tag[n_param], "Y_He", tag_length);
      params_addr[n_param] = &(run_params->physics).Y_He;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "ReionMaxHeatingRedshift");
      params_addr[n_param] = &(run_params->physics).ReionMaxHeatingRedshift;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "LXrayGal");
      params_addr[n_param] = &(run_params->physics).LXrayGal;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "NuXrayGalThreshold");
      params_addr[n_param] = &(run_params->physics).NuXrayGalThreshold;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SpecIndexXrayGal");
      params_addr[n_param] = &(run_params->physics).SpecIndexXrayGal;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "LXrayQSO");
      params_addr[n_param] = &(run_params->physics).LXrayQSO;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "NuXrayQSOThreshold");
      params_addr[n_param] = &(run_params->physics).NuXrayQSOThreshold;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "SpecIndexXrayQSO");
      params_addr[n_param] = &(run_params->physics).SpecIndexXrayQSO;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "NuXraySoftCut");
      params_addr[n_param] = &(run_params->physics).NuXraySoftCut;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      strcpy(params_tag[n_param], "NuXrayMax");
      params_addr[n_param] = &(run_params->physics).NuXrayMax;
      required_tag[n_param] = 1;
      params_type[n_param++] = PARAM_TYPE_DOUBLE;

      hdf5props->params_count = n_param;
    }

    // Initialise used_tag
    for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
      used_tag[ii] = 0;

    // Parse the user parameter file first
    n_entries = parse_paramfile(fname, entry);
    store_params(entry, n_entries, params_tag, n_param, used_tag, params_type, params_addr);

    // Now parse the default parameter file
    n_entries = parse_paramfile(run_params->DefaultsFile, entry);
    store_params(entry, n_entries, params_tag, n_param, used_tag, params_type, params_addr);

    // Check to see if we are missing any required parameters
    for (ii = 0; ii < n_param; ii++)
      if ((used_tag[ii] == 0) && (required_tag[ii] == 1)) {
        mlog_error("I miss a value for tag '%s' in parameter file '%s'.", params_tag[ii], fname);
        ABORT(EXIT_FAILURE);
      }

    for (ii = 0; ii < n_param; ii++) {
      if (used_tag[ii] == 1) {
        mlog("%35s\t", MLOG_MESG, params_tag[ii]);

        switch (params_type[ii]) {
          case PARAM_TYPE_DOUBLE:
            mlog("%g", MLOG_MESG | MLOG_CONT, *((double*)(params_addr[ii])));
            break;

          case PARAM_TYPE_FLOAT:
            mlog("%g", MLOG_MESG | MLOG_CONT, *((float*)(params_addr[ii])));
            break;

          case PARAM_TYPE_STRING:
            mlog("%s", MLOG_MESG | MLOG_CONT, (char*)params_addr[ii]);
            break;

          case PARAM_TYPE_INT:
            mlog("%d", MLOG_MESG | MLOG_CONT, *((int*)(params_addr[ii])));
            break;

          case PARAM_TYPE_LONGLONG:
            mlog("%lld", MLOG_MESG | MLOG_CONT, *((long long*)(params_addr[ii])));
            break;

          default:
            mlog_error("Unknown param type.");
            break;
        }
      }
    }

    ii = (int)strlen(run_params->OutputDir);
    if (ii > 0)
      if (run_params->OutputDir[ii - 1] != '/')
        strcat(run_params->OutputDir, "/");

    check_problem_params(run_params);
  } // END if(run_globals.mpi_rank==0)

  // If running mpi then broadcast the run parameters to all cores
  MPI_Bcast(run_params, sizeof(run_params_t), MPI_BYTE, 0, run_globals.mpi_comm);
  MPI_Bcast(&(run_globals.units), sizeof(run_units_t), MPI_BYTE, 0, run_globals.mpi_comm);
}
