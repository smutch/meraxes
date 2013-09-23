#include "meraxes.h"
#include "yaml.h"

#define STRING 101
#define INT    102
#define DOUBLE 103

static void inline store_param(
  char *value,
  int   used_tag[MAXTAGS],
  int   params_id[MAXTAGS],  
  void *params_addr[MAXTAGS],
  int  *tag_index)
{
  switch(params_id[*tag_index])
  {
    case DOUBLE:
      *((double *) params_addr[*tag_index]) = atof(value);
      break;
    case STRING:
      strcpy(params_addr[*tag_index], value);
      break;
    case INT:
      *((int *) params_addr[*tag_index]) = atoi(value);
      break;
  }
  used_tag[*tag_index] = 1;
  *tag_index=-1;
}


static void inline get_tag_index(
  int  *tag_index,       
  int   n_param,         
  char *value,           
  char  tag[MAXTAGS][50])
{
  *tag_index=-1;
  for(int i=0; i<n_param; i++)
  {
    if(strcmp(value, tag[i])==0)
    {
      *tag_index = i;
      break;
    }
  }
}


static void parse_param_file(
    char *fname,               
    char  tag[MAXTAGS][50],    
    int   n_param,             
    int   used_tag[MAXTAGS],   
    int   params_id[MAXTAGS],  
    void *params_addr[MAXTAGS])
{

  /*
   * Parse an input parameter file.
   */

  yaml_parser_t parser;
  yaml_document_t document;
  yaml_node_t *node;
  int node_index = 1;
  int tag_index = -2;

  // Create the Parser object. 
  yaml_parser_initialize(&parser);

  // Set a file input. 
  FILE *fd = fopen(fname, "rb");
  yaml_parser_set_input_file(&parser, fd);

  // Load the document
  yaml_parser_load(&parser, &document);

  while((node = yaml_document_get_node(&document, node_index))!=NULL)
  {
    switch(node->type)
    {
      case YAML_SCALAR_NODE:
        get_tag_index(&tag_index, n_param, (char *)node->data.scalar.value, tag);
        if(tag_index<0)
        {
          SID_log_error("in file %s:   Tag '%s' not allowed or multiple defined.", fname, node->data.scalar.value);
          yaml_document_delete(&document);
          yaml_parser_delete(&parser);
          fclose(fd);
          ABORT(EXIT_FAILURE);
        }
        node_index++;
        node = yaml_document_get_node(&document, node_index);
        if(node==NULL)
          break;
        else if(node->type == YAML_SCALAR_NODE)
        {
          store_param((char *)node->data.scalar.value, used_tag, params_id, params_addr, &tag_index);
        } else
          node_index--;
        break;

      case YAML_SEQUENCE_NODE:
        break;

      case YAML_NO_EVENT:
        break;

      case YAML_MAPPING_NODE:
        break;
    }
    node_index++;
  }

  // Delete the document
  yaml_document_delete(&document);

  // Delete the parser
  yaml_parser_delete(&parser);

  // Close the file
  fclose(fd);

}


void read_parameter_file(run_globals_struct *run_globals, char *fname)
{
  int i, n_param;
  int user_used_tag[MAXTAGS], defaults_used_tag[MAXTAGS], required_tag[MAXTAGS];
  char defaults_file[STRLEN];
  int    params_id[MAXTAGS];
  void   *params_addr[MAXTAGS];
  char   params_tag[MAXTAGS][50];

  run_params_struct *run_params = &(run_globals->params);

  // Initialise values and arrays
  n_param = 0;
  for (i = 0; i < MAXTAGS; i++) {
    user_used_tag[i] = 0;
    defaults_used_tag[i] = 0;
    required_tag[i] = 0;
  }

  if(SID.My_rank == 0)
    printf("\nreading parameter file:\n\n");

  strcpy(params_tag[n_param], "DefaultsFile");
  params_addr[n_param] = defaults_file;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "OutputDir");
  params_addr[n_param] = run_params->OutputDir;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "PhotometricTablesDir");
  params_addr[n_param] = run_params->PhotometricTablesDir;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "SSPModel");
  params_addr[n_param] = run_params->SSPModel;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "IMF");
  params_addr[n_param] = run_params->IMF;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "MagSystem");
  params_addr[n_param] = run_params->MagSystem;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "MagBands");
  params_addr[n_param] = run_params->MagBands;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "FileNameGalaxies");
  params_addr[n_param] = run_params->FileNameGalaxies;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "SimName");
  params_addr[n_param] = run_params->SimName;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "SimulationDir");
  params_addr[n_param] = run_params->SimulationDir;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "NEverySnap");
  params_addr[n_param] = &(run_params->NEverySnap);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "NScanSnap");
  params_addr[n_param] = &(run_params->NScanSnap);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "FileWithOutputSnaps");
  params_addr[n_param] = run_params->FileWithOutputSnaps;
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "FilesPerSnapshot");
  params_addr[n_param] = &(run_params->FilesPerSnapshot);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "TotalSimSnaps");
  params_addr[n_param] = &(run_params->TotalSimSnaps);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "FirstFile");
  params_addr[n_param] = &(run_params->FirstFile);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "LastFile");
  params_addr[n_param] = &(run_params->LastFile);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;
  
  strcpy(params_tag[n_param], "BoxSize");
  params_addr[n_param] = &(run_params->BoxSize);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "VolumeFactor");
  params_addr[n_param] = &(run_params->VolumeFactor);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "RandomSeed");
  params_addr[n_param] = &(run_params->RandomSeed);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "ThreshMajorMerger");
  params_addr[n_param] = &(run_params->ThreshMajorMerger);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "RecycleFraction");
  params_addr[n_param] = &(run_params->RecycleFraction);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "UnitVelocity_in_cm_per_s");
  params_addr[n_param] = &(run_globals->units.UnitVelocity_in_cm_per_s);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "UnitLength_in_cm");
  params_addr[n_param] = &(run_globals->units.UnitLength_in_cm);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "UnitMass_in_g");
  params_addr[n_param] = &(run_globals->units.UnitMass_in_g);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "Hubble_h");
  params_addr[n_param] = &(run_params->Hubble_h);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "BaryonFrac");
  params_addr[n_param] = &(run_params->BaryonFrac);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "OmegaM");
  params_addr[n_param] = &(run_params->OmegaM);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "OmegaK");
  params_addr[n_param] = &(run_params->OmegaK);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "OmegaLambda");
  params_addr[n_param] = &(run_params->OmegaLambda);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "PartMass");
  params_addr[n_param] = &(run_params->PartMass);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "MergerTimeFactor");
  params_addr[n_param] = &(run_params->MergerTimeFactor);
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  // Physics params

  strcpy(params_tag[n_param], "funcprop");
  params_addr[n_param] = &(run_params->physics).funcprop;
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "peak");
  params_addr[n_param] = &(run_params->physics).peak;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "peak_evo");
  params_addr[n_param] = &(run_params->physics).peak_evo;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "sigma");
  params_addr[n_param] = &(run_params->physics).sigma;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;
  
  strcpy(params_tag[n_param], "sigma_evo");
  params_addr[n_param] = &(run_params->physics).sigma_evo;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;
  
  strcpy(params_tag[n_param], "stellarfrac");
  params_addr[n_param] = &(run_params->physics).stellarfrac;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;
  
  strcpy(params_tag[n_param], "stellarfrac_evo");
  params_addr[n_param] = &(run_params->physics).stellarfrac_evo;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;
  
  strcpy(params_tag[n_param], "bhgrowthfactor");
  params_addr[n_param] = &(run_params->physics).bhgrowthfactor;
  required_tag[n_param] = 1;
  params_id[n_param++] = DOUBLE;

  strcpy(params_tag[n_param], "TOCF_Flag");
  params_addr[n_param] = &(run_params->TOCF_Flag);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;

  strcpy(params_tag[n_param], "TOCF_LogFileDir");
  params_addr[n_param] = &(run_params->TOCF_LogFileDir);
  required_tag[n_param] = 1;
  params_id[n_param++] = STRING;

  strcpy(params_tag[n_param], "TOCF_NThreads");
  params_addr[n_param] = &(run_params->TOCF_NThreads);
  required_tag[n_param] = 1;
  params_id[n_param++] = INT;


  
  // N.B. This part of the code is wasteful and should be updated!!! 

  // Parse the user parameter file first to get the default parameters file.
  parse_param_file(fname, params_tag, n_param, user_used_tag, params_id, params_addr);

  // Reset the user defined tags flags
  for (i = 0; i < MAXTAGS; i++)
    user_used_tag[i] = 0;

  // Now parse the default parameter file
  parse_param_file(defaults_file, params_tag, n_param, defaults_used_tag, params_id, params_addr);
  
  // Finally - parse the user file again to override any defaults.
  parse_param_file(fname, params_tag, n_param, user_used_tag, params_id, params_addr);

  for(i = 0; i < n_param; i++)
  {
    if((user_used_tag[i]==0) && (defaults_used_tag[i]==0) && (required_tag[i]==1))
    {
      SID_log_error("I miss a value for tag '%s' in parameter file '%s'.", params_tag[i], fname);
    }
  }

  if(SID.My_rank == 0)
  {
    for (i = 0; i < n_param; i++) {
      printf("%35s\t", params_tag[i]);
      
      switch (params_id[i])
        {
          case DOUBLE:
          printf("%g\n", *((double *) (params_addr[i])));
          break;
          case STRING:
          printf("%s\n", (char *) params_addr[i]);
          break;
          case INT:
          printf("%d\n", *((int *) (params_addr[i])));
          break;
        }

      if (i==(n_param-9))
        printf("\t\t%35s\n", "--- physics parameters ---");

    }
    printf("\n");
  }

  i = strlen(run_params->OutputDir);
  if(i > 0)
    if(run_params->OutputDir[i-1] != '/')
      strcat(run_params->OutputDir, "/");


}
