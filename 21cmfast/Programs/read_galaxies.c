typedef struct gal_struct gal_struct;
struct gal_struct{
  float Mvir;
  float Pos[3];
  float StellarMass;
};

void read_galaxies(char *fname, int snapshot, float M_MIN, float **stellar_mass, float **pos)
{

  // Read in all galaxies within halos above a given mass threshold.

  herr_t status;
  hid_t file_id, group_id;
  hsize_t nfields;
  hsize_t start = 0;
  hsize_t nrecords =10;
  char group_name[128];
  char **field_names;
  char *requested_names = {"Mvir,StellarMass"};
  int count = 0;

  // Open the file and group
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  sprintf(group_name, "Snap%03d", snapshot);
  group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

  // Get the table info
  status = H5TBget_table_info(group_id, "Galaxies", &nfields, &nrecords);
  printf("Reading in %d galaxies... = %d\n", (int)nrecords);

  // Setup the offsets and sizes and buffer for the read
  gal_struct gals[(int)nrecords];
  size_t gal_offsets[3]=
  {
    HOFFSET( gal_struct, Mvir ),
    HOFFSET( gal_struct, Pos ),
    HOFFSET( gal_struct, StellarMass )
  };
  size_t gal_sizes[3]=
  {
    sizeof((*galaxy_field)[0].Mvir),
    sizeof((*galaxy_field)[0].Pos),
    sizeof((*galaxy_field)[0].StellarMass)
  };
  size_t type_size = sizeof(gal);

  // Read the galaxies
  status = H5TBread_fields_name(group_id, "Galaxies", requested_names, start, nrecords, type_size, gal_offsets, gal_sizes, gals);

  // Weed out all those galaxies in halos below our M_MIN value
  for(int ii=0,count=0; ii<nrecords; ii++)
  {
    gals[ii].Mvir *= 1.e10;
    gals[ii].StellarMass *= 1.e10;
    if(gals[ii].Mvir > M_MIN)
      count++;
  }

  *stellar_mass = malloc(sizeof(float) * (size_t)count);
  *pos = malloc(sizeof(float) * (size_t)count);
  for(int ii=0,count=0; ii<nrecords; ii++)
    if(gals[ii].Mvir > M_MIN)
    {
      (*stellar_mass)[count++] = gals[ii].StellarMass;
      (*pos)[count++] = gals[ii].Pos;
    }

  // Close the groups and file
  H5Gclose(group_id);
  H5Fclose(file_id);

  return;
}


