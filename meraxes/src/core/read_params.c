#include "meraxes.h"
#include "parse_paramfile.h"

static void inline store_params(
  entry_t entry[PARAM_MAX_ENTRIES],
  int   n_entries,
  char  params_tag[PARAM_MAX_ENTRIES][128],
  int   n_param,
  int   used_tag[PARAM_MAX_ENTRIES],
  int   params_type[PARAM_MAX_ENTRIES],
  void *params_addr[PARAM_MAX_ENTRIES],
  run_params_t *run_params)
{
  int level;
  int tag_index;
  int temp;
  char prefix[16] = "\0";
  char key[STRLEN];

  for(int i_entry=0; i_entry<n_entries; i_entry++)
  {
    // DEBUG
    // SID_log("Checking %s", SID_LOG_COMMENT, entry[i_entry].key);

    // reset prefix if we have descended an indentation level
    if (entry[i_entry].level < level)
      *prefix = '\0';

    strcpy(key, prefix);
    strcat(key, entry[i_entry].key);
    level = entry[i_entry].level;

    // DEBUG
    // SID_log("level = %d :: prefix = %s", SID_LOG_COMMENT, level, prefix);

    tag_index=-1;
    for(int ii=0; ii<n_param; ii++)
    {
      if(strcmp(key, params_tag[ii])==0)
      {
        tag_index = ii;
        break;
      }
    }

    if(strcmp(key, "TOCF_Flag")==0)
    {
      temp = atoi(entry[i_entry].value);
      if(used_tag[tag_index]==0)
      {
        *((int *) params_addr[tag_index]) = atoi(entry[i_entry].value);
        used_tag[tag_index] = 1;
      }
      if(temp!=1)
      {
        SID_log("Skipping TOCF params block...", SID_LOG_COMMENT);
        temp = entry[i_entry++].level;
        while(entry[i_entry++].level>temp);
        i_entry-=2;
      } else
        sprintf(prefix, "TOCF_");
    }

    if (tag_index<0)
    {
      SID_log_warning("%s is an unrecognised parameter. Ignoring...", SID_LOG_COMMENT, entry[i_entry].key);
      continue;
    }

    if(used_tag[tag_index]==1)
      continue;

    switch(params_type[tag_index])
    {
      case PARAM_TYPE_DOUBLE:
        *((double *) params_addr[tag_index]) = atof(entry[i_entry].value);
        break;
      case PARAM_TYPE_FLOAT:
        *((float *) params_addr[tag_index]) = atof(entry[i_entry].value);
        break;
      case PARAM_TYPE_STRING:
        strcpy(params_addr[tag_index], entry[i_entry].value);
        break;
      case PARAM_TYPE_INT:
        *((int *) params_addr[tag_index]) = atoi(entry[i_entry].value);
        break;
    }
    used_tag[tag_index] = 1;
  }
}


void read_parameter_file(run_globals_t *run_globals, char *fname)
{
  run_params_t *run_params = &(run_globals->params);

  if(SID.My_rank == 0)
  {
    int i, n_param;
    int used_tag[PARAM_MAX_ENTRIES], required_tag[PARAM_MAX_ENTRIES];
    char defaults_file[STRLEN];
    int    params_type[PARAM_MAX_ENTRIES];
    void   *params_addr[PARAM_MAX_ENTRIES];
    char   params_tag[PARAM_MAX_ENTRIES][128];
    entry_t entry[PARAM_MAX_ENTRIES];
    int    n_entries;


    // Initialise values and arrays
    n_param = 0;
    for (i = 0; i < PARAM_MAX_ENTRIES; i++) {
      used_tag[i] = 0;
      required_tag[i] = 0;
    }

    printf("\nreading parameter file:\n\n");

    strcpy(params_tag[n_param], "DefaultsFile");
    params_addr[n_param] = defaults_file;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "OutputDir");
    params_addr[n_param] = run_params->OutputDir;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "PhotometricTablesDir");
    params_addr[n_param] = run_params->PhotometricTablesDir;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "CoolingFuncsDir");
    params_addr[n_param] = run_params->CoolingFuncsDir;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "SSPModel");
    params_addr[n_param] = run_params->SSPModel;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "IMF");
    params_addr[n_param] = run_params->IMF;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "MagSystem");
    params_addr[n_param] = run_params->MagSystem;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "MagBands");
    params_addr[n_param] = run_params->MagBands;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "FileNameGalaxies");
    params_addr[n_param] = run_params->FileNameGalaxies;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "SimName");
    params_addr[n_param] = run_params->SimName;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "SimulationDir");
    params_addr[n_param] = run_params->SimulationDir;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "CatalogFilePrefix");
    params_addr[n_param] = run_params->CatalogFilePrefix;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "NEverySnap");
    params_addr[n_param] = &(run_params->NEverySnap);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "NScanSnap");
    params_addr[n_param] = &(run_params->NScanSnap);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "FileWithOutputSnaps");
    params_addr[n_param] = run_params->FileWithOutputSnaps;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "TotalSimSnaps");
    params_addr[n_param] = &(run_params->TotalSimSnaps);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "FirstFile");
    params_addr[n_param] = &(run_params->FirstFile);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "LastFile");
    params_addr[n_param] = &(run_params->LastFile);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "BoxSize");
    params_addr[n_param] = &(run_params->BoxSize);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "VolumeFactor");
    params_addr[n_param] = &(run_params->VolumeFactor);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "RandomSeed");
    params_addr[n_param] = &(run_params->RandomSeed);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "ForestIDFile");
    params_addr[n_param] = &(run_params->ForestIDFile);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_STRING;
    *(run_params->ForestIDFile) = '\0';

    strcpy(params_tag[n_param], "MultipleRunsParamFile");
    params_addr[n_param] = &(run_params->MultipleRunsParamFile);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_STRING;
    *(run_params->MultipleRunsParamFile) = '\0';

    strcpy(params_tag[n_param], "RandomSeed");
    params_addr[n_param] = &(run_params->RandomSeed);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "ThreshMajorMerger");
    params_addr[n_param] = &(run_params->ThreshMajorMerger);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "RecycleFraction");
    params_addr[n_param] = &(run_params->RecycleFraction);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "UnitVelocity_in_cm_per_s");
    params_addr[n_param] = &(run_globals->units.UnitVelocity_in_cm_per_s);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "UnitLength_in_cm");
    params_addr[n_param] = &(run_globals->units.UnitLength_in_cm);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "UnitMass_in_g");
    params_addr[n_param] = &(run_globals->units.UnitMass_in_g);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "Hubble_h");
    params_addr[n_param] = &(run_params->Hubble_h);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "BaryonFrac");
    params_addr[n_param] = &(run_params->BaryonFrac);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "OmegaM");
    params_addr[n_param] = &(run_params->OmegaM);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "OmegaK");
    params_addr[n_param] = &(run_params->OmegaK);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "OmegaR");
    params_addr[n_param] = &(run_params->OmegaR);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "OmegaLambda");
    params_addr[n_param] = &(run_params->OmegaLambda);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "Sigma8");
    params_addr[n_param] = &(run_params->Sigma8);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "wLambda");
    params_addr[n_param] = &(run_params->wLambda);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "SpectralIndex");
    params_addr[n_param] = &(run_params->SpectralIndex);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "PartMass");
    params_addr[n_param] = &(run_params->PartMass);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "MergerTimeFactor");
    params_addr[n_param] = &(run_params->MergerTimeFactor);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    // Physics params

    strcpy(params_tag[n_param], "SfEfficiency");
    params_addr[n_param] = &(run_params->physics).SfEfficiency;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "SfRecycleFraction");
    params_addr[n_param] = &(run_params->physics).SfRecycleFraction;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_z_re");
    params_addr[n_param] = &(run_params->physics).reion_z_re;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_delta_z_re");
    params_addr[n_param] = &(run_params->physics).reion_delta_z_re;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_delta_z_sc");
    params_addr[n_param] = &(run_params->physics).reion_delta_z_sc;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_T0");
    params_addr[n_param] = &(run_params->physics).reion_T0;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_Tcool");
    params_addr[n_param] = &(run_params->physics).reion_Tcool;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_Nion_phot_per_bary");
    params_addr[n_param] = &(run_params->physics).reion_Nion_phot_per_bary;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_escape_frac");
    params_addr[n_param] = &(run_params->physics).reion_escape_frac;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "reion_mean_n_rec");
    params_addr[n_param] = &(run_params->physics).reion_mean_n_rec;
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "TOCF_Flag");
    params_addr[n_param] = &(run_params->TOCF_Flag);
    required_tag[n_param] = 1;
    params_type[n_param++] = PARAM_TYPE_INT;

#ifdef USE_TOCF
    strcpy(params_tag[n_param], "TOCF_LogFileDir");
    params_addr[n_param] = &(tocf_params.logfile_dir);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_STRING;

    strcpy(params_tag[n_param], "TOCF_dim");
    params_addr[n_param] = &(tocf_params.dim);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "TOCF_HII_dim");
    params_addr[n_param] = &(tocf_params.HII_dim);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "TOCF_numcores");
    params_addr[n_param] = &(tocf_params.numcores);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "TOCF_ram");
    params_addr[n_param] = &(tocf_params.ram);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_ion_tvir_min");
    params_addr[n_param] = &(tocf_params.ion_tvir_min);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_DOUBLE;

    strcpy(params_tag[n_param], "TOCF_HII_eff_factor");
    params_addr[n_param] = &(tocf_params.HII_eff_factor);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_r_bubble_min");
    params_addr[n_param] = &(tocf_params.r_bubble_min);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_r_bubble_max");
    params_addr[n_param] = &(tocf_params.r_bubble_max);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_delta_r_HII_factor");
    params_addr[n_param] = &(tocf_params.delta_r_HII_factor);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_gamma_halo_bias");
    params_addr[n_param] = &(tocf_params.gamma_halo_bias);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_m_0_sm");
    params_addr[n_param] = &(tocf_params.m_0_sm);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_a_sm");
    params_addr[n_param] = &(tocf_params.a_sm);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_b_sm");
    params_addr[n_param] = &(tocf_params.b_sm);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_c_sm");
    params_addr[n_param] = &(tocf_params.c_sm);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_d_sm");
    params_addr[n_param] = &(tocf_params.d_sm);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_FLOAT;

    strcpy(params_tag[n_param], "TOCF_uvb_feedback");
    params_addr[n_param] = &(tocf_params.uvb_feedback);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "TOCF_compute_mfp");
    params_addr[n_param] = &(tocf_params.compute_mfp);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_INT;

    strcpy(params_tag[n_param], "TOCF_HII_filter");
    params_addr[n_param] = &(tocf_params.HII_filter);
    required_tag[n_param] = 0;
    params_type[n_param++] = PARAM_TYPE_INT;
#endif


    // Parse the user parameter file first
    n_entries = parse_paramfile(fname, entry);
    store_params(entry, n_entries, params_tag, n_param, used_tag, params_type, params_addr, run_params);

    // Now parse the default parameter file
    n_entries = parse_paramfile(defaults_file, entry);
    store_params(entry, n_entries, params_tag, n_param, used_tag, params_type, params_addr, run_params);

    // Check to see if we are missing any required parameters
    for(i = 0; i < n_param; i++)
    {
      if((used_tag[i]==0) && (required_tag[i]==1))
      {
        SID_log_error("I miss a value for tag '%s' in parameter file '%s'.", params_tag[i], fname);
      }
    }

    if(SID.My_rank == 0)
    {
      for (i = 0; i < n_param; i++) {

        if(used_tag[i]==1)
        {
          printf("%35s\t", params_tag[i]);

          switch (params_type[i])
          {
            case PARAM_TYPE_DOUBLE:
              printf("%g\n", *((double *) (params_addr[i])));
              break;
            case PARAM_TYPE_FLOAT:
              printf("%g\n", *((float *) (params_addr[i])));
              break;
            case PARAM_TYPE_STRING:
              printf("%s\n", (char *) params_addr[i]);
              break;
            case PARAM_TYPE_INT:
              printf("%d\n", *((int *) (params_addr[i])));
              break;
          }

          // if (i==(n_param-12))
          //   printf("\t\t%35s\n", "--- physics parameters ---");

          // if (i==(n_param-4))
          //   printf("\t\t%35s\n", "--- TOCF parameters ---");
        }
      }
    }

    i = strlen(run_params->OutputDir);
    if(i > 0)
      if(run_params->OutputDir[i-1] != '/')
        strcat(run_params->OutputDir, "/");

  }  // END if(SID.My_rank==0)

  // If running mpi then broadcast the run parameters to all cores
  SID_Bcast(run_params, sizeof(run_params_t), 0, SID.COMM_WORLD);
  SID_Bcast(&(run_globals->units), sizeof(run_units_t), 0, SID.COMM_WORLD);
#ifdef USE_TOCF
  SID_Bcast(&tocf_params, sizeof(tocf_params_t), 0, SID.COMM_WORLD);
#endif

}
