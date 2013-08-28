#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include "parameter_files/ANAL_PARAMS.H"
#include "parameter_files/INIT_PARAMS.H"

typedef struct gal_struct gal_struct;
struct gal_struct{
  float Pos[3];
  float Mvir;
  float StellarMass;
};

int read_galaxies(char *fname, int snapshot, gal_struct **gals)
{

  // Read in all galaxies within halos above a given mass threshold.

  herr_t status;
  hid_t file_id, group_id;
  hsize_t nfields;
  hsize_t start = 0;
  hsize_t nrecords =10;
  char group_name[128];
  char **field_names;
  char *requested_fields = {"Pos,Mvir,StellarMass"};
  int count = 0;

  // Open the file and group
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(group_name, "Snap%03d", snapshot);
  group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

  // Get the table info
  status = H5TBget_table_info(group_id, "Galaxies", &nfields, &nrecords);
  printf("Reading in %d galaxies... \n", (int)nrecords);

  // Setup the offsets and sizes and buffer for the read
  *gals = malloc(sizeof(gal_struct) * (size_t)nrecords);
  size_t gal_offsets[3]=
  {
    HOFFSET( gal_struct, Pos ),
    HOFFSET( gal_struct, Mvir ),
    HOFFSET( gal_struct, StellarMass )
  };
  size_t gal_sizes[3]=
  {
    sizeof((*gals)[0].Pos),
    sizeof((*gals)[0].Mvir),
    sizeof((*gals)[0].StellarMass)
  };
  size_t type_size = sizeof(gal_struct);

  // Read the galaxies
  status = H5TBread_fields_name(group_id, "Galaxies", requested_fields, start, nrecords, type_size, gal_offsets, gal_sizes, *gals);

  // Close the groups and file
  H5Gclose(group_id);
  H5Fclose(file_id);
  
  printf("...done\n");

  return (int)nrecords;
}

int generate_stellarmass_field(char *fname, int snapshot, float M_MIN, float *smfield)
{
  gal_struct *gals;
  float total_sm = 0.;
  int ngals;

  // Read in the galaxies
  printf("Creating stellar mass grid...\n");
  ngals = read_galaxies(fname, snapshot, &gals);
  if(ngals<=0)
    return 1;

  // Loop through and construct the stellar mass field
  // Note that Meraxes already outputs in h=0.705 units so no need to deal with
  // that here...
  for(int ii=0; ii<ngals; ii++)
    if ((gals[ii].Mvir*1.e10) >= M_MIN)
    {
      *(smfield + HII_R_FFT_INDEX(
            (int)(gals[ii].Pos[0]/BOX_LEN*HII_DIM),
            (int)(gals[ii].Pos[1]/BOX_LEN*HII_DIM),
            (int)(gals[ii].Pos[2]/BOX_LEN*HII_DIM))) 
        += (float)(gals[ii].StellarMass*1.e10);
      total_sm += gals[ii].StellarMass*1.e10;
    }

  printf("total stellar mass in grid = %.3e Msol\n", total_sm);
  printf("...done\n");

  // DEBUG
  // fprintf(stderr, "Writing smfield in test.hdf5...\n");
  // hid_t file_id = H5Fcreate("test.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  // hsize_t dims[] = {HII_TOT_NUM_PIXELS};
  // H5LTmake_dataset_float(file_id, "smfield", 1, dims, smfield);
  // dims[0] = ngals;
  // float temp[ngals];
  // for(int ii=0; ii<ngals; ii++)
  //   temp[ii] = gals[ii].Mvir*1.e10;
  // H5LTmake_dataset_float(file_id, "mvir", 1, dims, temp);
  // for(int ii=0; ii<ngals; ii++)
  //   temp[ii] = gals[ii].StellarMass*1.e10;
  // H5LTmake_dataset_float(file_id, "stellarmass", 1, dims, temp);
  // H5LTset_attribute_float(file_id, "/", "M_MIN", &M_MIN, 1);
  // H5Fclose(file_id);
  // fprintf(stderr, "...done\n");

  // Free the galaxies
  free(gals);

  return 0;
}
