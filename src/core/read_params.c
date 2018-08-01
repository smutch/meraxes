#include "meraxes.h"
#include "parse_paramfile.h"

static void check_problem_params(run_params_t* run_params)
{
    if (run_params->NSteps != 1) {
        mlog_error("The current version of the code only works if NSteps = 1. Sorry! Exiting...");
        ABORT(EXIT_FAILURE);
    }
}

static void store_params(entry_t entry[123],
    int n_entries,
    char** params_tag,
    int n_param,
    int used_tag[123],
    int* params_type,
    void** params_addr)
{
    int level = 0;
    char prefix[16] = "\0";
    char key[STRLEN];

    for (int i_entry = 0; i_entry < n_entries; i_entry++) {
        // DEBUG
        // mlog("Checking %s", MLOG_MESG, entry[i_entry].key);

        // reset prefix if we have descended an indentation level
        if (entry[i_entry].level < level)
            *prefix = '\0';

        strlcpy(key, prefix, STRLEN);
        strlcat(key, entry[i_entry].key, STRLEN);
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
            strlcpy(params_addr[tag_index], entry[i_entry].value, STRLEN);
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
            hdf5props->params_tag = malloc(sizeof(char*) * PARAM_MAX_ENTRIES);
            for (ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
                hdf5props->params_tag[ii] = malloc(sizeof(char) * 128);
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

            strlcpy(params_tag[n_param], "DefaultsFile", STRLEN);
            params_addr[n_param] = run_params->DefaultsFile;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "OutputDir", STRLEN);
            params_addr[n_param] = run_params->OutputDir;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "PhotometricTablesDir", STRLEN);
            params_addr[n_param] = run_params->PhotometricTablesDir;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "CoolingFuncsDir", STRLEN);
            params_addr[n_param] = run_params->CoolingFuncsDir;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "SSPModel", STRLEN);
            params_addr[n_param] = run_params->SSPModel;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "IMF", STRLEN);
            params_addr[n_param] = run_params->IMF;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "MagSystem", STRLEN);
            params_addr[n_param] = run_params->MagSystem;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "MagBands", STRLEN);
            params_addr[n_param] = run_params->MagBands;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "FileNameGalaxies", STRLEN);
            params_addr[n_param] = run_params->FileNameGalaxies;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "SimName", STRLEN);
            params_addr[n_param] = run_params->SimName;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "TreesID", STRLEN);
            params_addr[n_param] = &(run_params->TreesID);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "SimulationDir", STRLEN);
            params_addr[n_param] = run_params->SimulationDir;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "CatalogFilePrefix", STRLEN);
            params_addr[n_param] = run_params->CatalogFilePrefix;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "FileWithOutputSnaps", STRLEN);
            params_addr[n_param] = run_params->FileWithOutputSnaps;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_STRING;

            strlcpy(params_tag[n_param], "NSteps", STRLEN);
            params_addr[n_param] = &(run_params->NSteps);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "BoxSize", STRLEN);
            params_addr[n_param] = &(run_params->BoxSize);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "VolumeFactor", STRLEN);
            params_addr[n_param] = &(run_params->VolumeFactor);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "RandomSeed", STRLEN);
            params_addr[n_param] = &(run_params->RandomSeed);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "ForestIDFile", STRLEN);
            params_addr[n_param] = &(run_params->ForestIDFile);
            required_tag[n_param] = 0;
            params_type[n_param++] = PARAM_TYPE_STRING;
            *(run_params->ForestIDFile) = '\0';

            strlcpy(params_tag[n_param], "MvirCritFile", STRLEN);
            params_addr[n_param] = &(run_params->MvirCritFile);
            required_tag[n_param] = 0;
            params_type[n_param++] = PARAM_TYPE_STRING;
            *(run_params->MvirCritFile) = '\0';

            strlcpy(params_tag[n_param], "MassRatioModifier", STRLEN);
            params_addr[n_param] = &(run_params->MassRatioModifier);
            required_tag[n_param] = 0;
            params_type[n_param++] = PARAM_TYPE_STRING;
            *(run_params->MassRatioModifier) = '\0';

            strlcpy(params_tag[n_param], "BaryonFracModifier", STRLEN);
            params_addr[n_param] = &(run_params->BaryonFracModifier);
            required_tag[n_param] = 0;
            params_type[n_param++] = PARAM_TYPE_STRING;
            *(run_params->BaryonFracModifier) = '\0';

            strlcpy(params_tag[n_param], "UnitVelocity_in_cm_per_s", STRLEN);
            params_addr[n_param] = &(run_globals.units.UnitVelocity_in_cm_per_s);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "UnitLength_in_cm", STRLEN);
            params_addr[n_param] = &(run_globals.units.UnitLength_in_cm);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "UnitMass_in_g", STRLEN);
            params_addr[n_param] = &(run_globals.units.UnitMass_in_g);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "Hubble_h", STRLEN);
            params_addr[n_param] = &(run_params->Hubble_h);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "BaryonFrac", STRLEN);
            params_addr[n_param] = &(run_params->BaryonFrac);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "OmegaM", STRLEN);
            params_addr[n_param] = &(run_params->OmegaM);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "OmegaK", STRLEN);
            params_addr[n_param] = &(run_params->OmegaK);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "OmegaR", STRLEN);
            params_addr[n_param] = &(run_params->OmegaR);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "OmegaLambda", STRLEN);
            params_addr[n_param] = &(run_params->OmegaLambda);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "Sigma8", STRLEN);
            params_addr[n_param] = &(run_params->Sigma8);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "wLambda", STRLEN);
            params_addr[n_param] = &(run_params->wLambda);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SpectralIndex", STRLEN);
            params_addr[n_param] = &(run_params->SpectralIndex);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "PartMass", STRLEN);
            params_addr[n_param] = &(run_params->PartMass);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "NPart", STRLEN);
            params_addr[n_param] = &(run_params->NPart);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_LONGLONG;

            strlcpy(params_tag[n_param], "MergerTimeFactor", STRLEN);
            params_addr[n_param] = &(run_params->physics.MergerTimeFactor);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "FlagSubhaloVirialProps", STRLEN);
            params_addr[n_param] = &(run_params->FlagSubhaloVirialProps);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "FlagInteractive", STRLEN);
            params_addr[n_param] = &(run_params->FlagInteractive);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "FlagMCMC", STRLEN);
            params_addr[n_param] = &(run_params->FlagMCMC);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            // Physics params

            strlcpy(params_tag[n_param], "EscapeFracDependency", STRLEN);
            params_addr[n_param] = &(run_params->physics).EscapeFracDependency;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "SfDiskVelOpt", STRLEN);
            params_addr[n_param] = &(run_params->physics).SfDiskVelOpt;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "SfPrescription", STRLEN);
            params_addr[n_param] = &(run_params->physics).SfPrescription;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_ReionizationModifier", STRLEN);
            params_addr[n_param] = &(run_params->physics).Flag_ReionizationModifier;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_BHFeedback", STRLEN);
            params_addr[n_param] = &(run_params->physics).Flag_BHFeedback;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_IRA", STRLEN);
            params_addr[n_param] = &(run_params->physics).Flag_IRA;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_FixDiskRadiusOnInfall", STRLEN);
            params_addr[n_param] = &(run_params->physics).Flag_FixDiskRadiusOnInfall;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_FixVmaxOnInfall", STRLEN);
            params_addr[n_param] = &(run_params->physics).Flag_FixVmaxOnInfall;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_ReheatToFOFGroupTemp", STRLEN);
            params_addr[n_param] = &(run_params->physics).Flag_ReheatToFOFGroupTemp;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "SfEfficiency", STRLEN);
            params_addr[n_param] = &(run_params->physics).SfEfficiency;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SfEfficiencyScaling", STRLEN);
            params_addr[n_param] = &(run_params->physics).SfEfficiencyScaling;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SfCriticalSDNorm", STRLEN);
            params_addr[n_param] = &(run_params->physics).SfCriticalSDNorm;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SfRecycleFraction", STRLEN);
            params_addr[n_param] = &(run_params->physics).SfRecycleFraction;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnEjectionEff", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnEjectionEff;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnEjectionScaling", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnEjectionScaling;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnEjectionNorm", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnEjectionNorm;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnReheatEff", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnReheatEff;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnReheatLimit", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnReheatLimit;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnReheatScaling", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnReheatScaling;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "SnReheatNorm", STRLEN);
            params_addr[n_param] = &(run_params->physics).SnReheatNorm;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReincorporationEff", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReincorporationEff;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "MaxCoolingMassFactor", STRLEN);
            params_addr[n_param] = &(run_params->physics).MaxCoolingMassFactor;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "Yield", STRLEN);
            params_addr[n_param] = &(run_params->physics).Yield;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "IMFSlope", STRLEN);
            params_addr[n_param] = &(run_params->physics).IMFSlope;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EnergyPerSN", STRLEN);
            params_addr[n_param] = &(run_params->physics).EnergyPerSN;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "IMFNormConst", STRLEN);
            params_addr[n_param] = &(run_params->physics).IMFNormConst;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "eta_SNII", STRLEN);
            params_addr[n_param] = &(run_params->physics).eta_SNII;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "frac_mass_SSP_above_SNII", STRLEN);
            params_addr[n_param] = &(run_params->physics).frac_mass_SSP_above_SNII;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ThreshMajorMerger", STRLEN);
            params_addr[n_param] = &((run_params->physics).ThreshMajorMerger);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "MinMergerStellarMass", STRLEN);
            params_addr[n_param] = &((run_params->physics).MinMergerStellarMass);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "MinMergerRatioForBurst", STRLEN);
            params_addr[n_param] = &((run_params->physics).MinMergerRatioForBurst);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "MergerBurstScaling", STRLEN);
            params_addr[n_param] = &((run_params->physics).MergerBurstScaling);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "MergerBurstFactor", STRLEN);
            params_addr[n_param] = &((run_params->physics).MergerBurstFactor);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "RadioModeEff", STRLEN);
            params_addr[n_param] = &(run_params->physics).RadioModeEff;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "QuasarModeEff", STRLEN);
            params_addr[n_param] = &(run_params->physics).QuasarModeEff;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "BlackHoleSeed", STRLEN);
            params_addr[n_param] = &(run_params->physics).BlackHoleSeed;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "BlackHoleGrowthRate", STRLEN);
            params_addr[n_param] = &(run_params->physics).BlackHoleGrowthRate;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EddingtonRatio", STRLEN);
            params_addr[n_param] = &(run_params->physics).EddingtonRatio;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "quasar_mode_scaling", STRLEN);
            params_addr[n_param] = &(run_params->physics).quasar_mode_scaling;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "quasar_open_angle", STRLEN);
            params_addr[n_param] = &(run_params->physics).quasar_open_angle;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSobacchi_Zre", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSobacchi_Zre;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSobacchi_DeltaZre", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSobacchi_DeltaZre;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSobacchi_DeltaZsc", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSobacchi_DeltaZsc;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSobacchi_T0", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSobacchi_T0;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionTcool", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionTcool;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionNionPhotPerBary", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionNionPhotPerBary;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "BlackHoleMassLimitReion", STRLEN);
            params_addr[n_param] = &(run_params->physics).BlackHoleMassLimitReion;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionGnedin_z0", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionGnedin_z0;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionGnedin_zr", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionGnedin_zr;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "Flag_PatchyReion", STRLEN);
            params_addr[n_param] = &(run_params->Flag_PatchyReion);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_OutputGrids", STRLEN);
            params_addr[n_param] = &(run_params->Flag_OutputGrids);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Flag_OutputGridsPostReion", STRLEN);
            params_addr[n_param] = &(run_params->Flag_OutputGridsPostReion);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "ReionGridDim", STRLEN);
            params_addr[n_param] = &(run_params->ReionGridDim);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "ReionRBubbleMin", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionRBubbleMin;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionRBubbleMax", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionRBubbleMax;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionDeltaRFactor", STRLEN);
            params_addr[n_param] = &(run_params->ReionDeltaRFactor);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionGammaHaloBias", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionGammaHaloBias;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EscapeFracNorm", STRLEN);
            params_addr[n_param] = &(run_params->physics).EscapeFracNorm;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EscapeFracOffset", STRLEN);
            params_addr[n_param] = &(run_params->physics).EscapeFracOffset;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EscapeFracScaling", STRLEN);
            params_addr[n_param] = &(run_params->physics).EscapeFracScaling;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EscapeFracBHNorm", STRLEN);
            params_addr[n_param] = &(run_params->physics).EscapeFracBHNorm;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "EscapeFracBHScaling", STRLEN);
            params_addr[n_param] = &(run_params->physics).EscapeFracBHScaling;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSMParam_m0", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSMParam_m0;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSMParam_a", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSMParam_a;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSMParam_b", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSMParam_b;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSMParam_c", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSMParam_c;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionSMParam_d", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionSMParam_d;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionUVBFlag", STRLEN);
            params_addr[n_param] = &(run_params->ReionUVBFlag);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "ReionFilterType", STRLEN);
            params_addr[n_param] = &(run_params->ReionFilterType);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "ReionPowerSpecDeltaK", STRLEN);
            params_addr[n_param] = &(run_params->ReionPowerSpecDeltaK);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionAlphaUV", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionAlphaUV;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionAlphaUVBH", STRLEN);
            params_addr[n_param] = &(run_params->physics).ReionAlphaUVBH;
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_DOUBLE;

            strlcpy(params_tag[n_param], "ReionRtoMFilterType", STRLEN);
            params_addr[n_param] = &(run_params->ReionRtoMFilterType);
            required_tag[n_param] = 1;
            params_type[n_param++] = PARAM_TYPE_INT;

            strlcpy(params_tag[n_param], "Y_He", STRLEN);
            params_addr[n_param] = &(run_params->physics).Y_He;
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

        if (run_globals.mpi_rank == 0) {
            for (ii = 0; ii < n_param; ii++)
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
