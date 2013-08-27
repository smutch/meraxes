#include <hdf5.h>
#include <hdf5_hl.h>

typedef struct gal_struct gal_struct;
struct gal_struct{
  float Mvir;
  float Pos[3];
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
  char *requested_fields = {"Mvir,StellarMass"};
  int count = 0;

  // Open the file and group
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(group_name, "Snap%03d", snapshot);
  group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

  // Get the table info
  status = H5TBget_table_info(group_id, "Galaxies", &nfields, &nrecords);
  printf("Reading in %d galaxies... = %d\n", (int)nrecords);

  // Setup the offsets and sizes and buffer for the read
  *gals = malloc(sizeof(gal_struct) * (size_t)nrecords);
  size_t gal_offsets[3]=
  {
    HOFFSET( gal_struct, Mvir ),
    HOFFSET( gal_struct, Pos ),
    HOFFSET( gal_struct, StellarMass )
  };
  size_t gal_sizes[3]=
  {
    sizeof((*gals)[0].Mvir),
    sizeof((*gals)[0].Pos),
    sizeof((*gals)[0].StellarMass)
  };
  size_t type_size = sizeof(gal_struct);

  // Read the galaxies
  status = H5TBread_fields_name(group_id, "Galaxies", requested_fields, start, nrecords, type_size, gal_offsets, gal_sizes, *gals);

  // Close the groups and file
  H5Gclose(group_id);
  H5Fclose(file_id);

  return (int)nrecords;
}

int generate_stellarmass_field(char *fname, int snapshot, float M_MIN, float *smfield)
{
  gal_struct *gals;
  int ngals;

  // Read in the galaxies
  ngals = read_galaxies(fname, snapshot, &gals);
  if(ngals<=0)
    return 1;

  // Loop through and construct the stellar mass field
  for(int ii=0; ii<ngals; ii++)
    if ((gals[ii].Mvir*1.e10/hlittle) >= M_MIN)
      *(smfield + HII_R_FFT_INDEX(
            (int)(gals[ii].Pos[0]/BOX_LEN*HII_DIM/hlittle),
            (int)(gals[ii].Pos[1]/BOX_LEN*HII_DIM/hlittle),
            (int)(gals[ii].Pos[2]/BOX_LEN*HII_DIM/hlittle))) 
        += (float)(gals[ii].StellarMass*1.e10/hlittle);

  // Free the galaxies
  free(gals);

  return 0;
}
