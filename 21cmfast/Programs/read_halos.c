#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "nbody.h"
#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/COSMOLOGY.H"

#ifdef DEBUG
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

/*
 * Read in the halo catalogues for a specified snapshot.
 */

typedef struct halo_catalog_struct halo_catalog_struct;
struct halo_catalog_struct{
  long long id_MBP;                    //!< ID of most bound particle in structure
  double    M_vir;                     //!< Bryan & Norman (ApJ 495, 80, 1998) virial mass [M_sol/h]
  int       n_particles;               //!< Number of particles in the structure
  float     position_COM[3];           //!< Centre-of-mass position      [Mpc/h]
  float     position_MBP[3];           //!< Most bound particle position [Mpc/h]
  float     velocity_COM[3];           //!< Centre-of-mass velocity      [km/s]
  float     velocity_MBP[3];           //!< Most bound particle velocity [km/s]
  float     R_vir;                     //!< Virial radius [Mpc/h]
  float     R_halo;                    //!< Distance of last halo particle from MBP [Mpc/h]
  float     R_max;                     //!< Radius of maximum circular velocity     [Mpc/h]
  float     V_max;                     //!< Maximum circular velocity               [km/s]
  float     sigma_v;                   //!< Total 3D velocity dispersion            [km/s]
  float     spin[3];                   //!< Specific angular momentum vector        [Mpc/h*km/s]
  float     q_triaxial;                //!< Triaxial shape parameter q=b/a
  float     s_triaxial;                //!< Triaxial shape parameter s=c/a
  float     shape_eigen_vectors[3][3]; //!< Normalized triaxial shape eigenvectors
  char      padding[8];                //!< Alignment padding
};

static void halo_catalog_filename(
  int   snapshot,  
  char *group_type,
  int   sub,       
  int  *i_layout,  
  char  fname[512])     
{

  bool flag_success = false;
  FILE *fin;

  // if we need to determine the filename structure...
  if (*i_layout==-1)
  {
    for (*i_layout=0; (*i_layout<4) && (flag_success==false); (*i_layout)++)
    {
      if (*i_layout==0)
        sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties/rockstar_fixed_%03d.catalog_%s_properties.%d", snapshot, group_type, snapshot, group_type, sub);
      else if (*i_layout==1)
        sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties/rockstar_fixed_%03d.catalog_%s_properties", snapshot, group_type, snapshot, group_type);
      else if (*i_layout==2)
        sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties.%d", snapshot, group_type, sub);
      else if (*i_layout==3)
        sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties", snapshot, group_type);
      
      if ((fin = fopen(fname, "rb"))!=NULL)
      {
        flag_success = true;
        fclose(fin);
        break;
      }
    }
  } 

  // ensure we have a valid i_layout value.
  if (*i_layout<0 && *i_layout>3)
  {
    fprintf(stderr, "cannot resolve catalogue filename.\n");
    ABORT(EXIT_FAILURE);
  }
   
  // provide the correct filename
  if (*i_layout==0)
    sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties/rockstar_fixed_%03d.catalog_%s_properties.%d", snapshot, group_type, snapshot, group_type, sub);
  else if (*i_layout==1)
    sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties/rockstar_fixed_%03d.catalog_%s_properties", snapshot, group_type, snapshot, group_type);
  else if (*i_layout==2)
    sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties.%d", snapshot, group_type, sub);
  else if (*i_layout==3)
    sprintf(fname, ROOT_PATH "/" SIM_NAME "/catalogs/rockstar_fixed_%03d.catalog_%s_properties", snapshot, group_type);

}


static void inline read_catalogs_header(
  FILE *fin,          
  int  *i_file,       
  int  *N_files,      
  int  *N_halos_file, 
  int  *N_halos_total )
{
  fread(i_file       , sizeof(int), 1, fin);
  fread(N_files      , sizeof(int), 1, fin);
  fread(N_halos_file , sizeof(int), 1, fin);
  fread(N_halos_total, sizeof(int), 1, fin);
}


int read_groups_field(
  int     snapshot,
  float  *halo_field,
  float   M_MIN)
{

  int                  flayout_switch=-1;
  FILE                *fin            = NULL;
  char                 fname[512];
  int                  dummy;
  int                  N_files;
  int                  N_halos_file;
  int                  N_halos_total;
  int                  halo_count = 0;
  halo_catalog_struct  catalog_halo;


  fprintf(stderr, "Reading halo groups...\n");

  // Read the first header to find out number of files and total number of halos
  halo_catalog_filename(snapshot, "groups", 0, &flayout_switch, fname);
  fin = fopen(fname, "rb");
  if (fin==NULL)
  {
    fprintf(stderr, "Failed to open file %s\n", fname);
    return -1;
  }
  read_catalogs_header(fin, &dummy, &N_files, &dummy, &N_halos_total);
  fclose(fin);

#ifdef DEBUG
  float *mass, (*pos)[3];
  mass = malloc(N_halos_total * sizeof(float));
  pos = malloc(N_halos_total * sizeof(float[3]));
#endif

  // Loop through each file and read in the halos
  for (int i_file=0; i_file<N_files; i_file++)
  {
    halo_catalog_filename(snapshot, "groups", i_file, &flayout_switch, fname);
    fin = fopen(fname, "rb");
    if (fin==NULL)
    {
      fprintf(stderr, "Failed to open file %s\n", fname);
      return -1;
    }
    read_catalogs_header(fin, &dummy, &dummy, &N_halos_file, &dummy);

    // Read in the halos
    for(int i_halo=0; i_halo<N_halos_file; i_halo++)
    {
      fread(&catalog_halo, sizeof(halo_catalog_struct), 1, fin);
      // fprintf(stderr, "Mass = %.2e :: R_vir(read) = %.2e :: Pos = [%.2f, %.2f, %.2f]\n", 
      //     catalog_halo.M_vir,
      //     catalog_halo.R_vir,
      //     catalog_halo.position_MBP[0],
      //     catalog_halo.position_MBP[1],
      //     catalog_halo.position_MBP[2]);
      // fprintf(stderr, "Rvir(calc) = %.2e\n", MtoR(catalog_halo.M_vir));
      // fprintf(stderr, "[x,y,z] = [%.2f, %.2f, %.2f]\n", 
      //     catalog_halo.position_MBP[0]/BOX_LEN*HII_DIM/hlittle,
      //     catalog_halo.position_MBP[1]/BOX_LEN*HII_DIM/hlittle,
      //     catalog_halo.position_MBP[2]/BOX_LEN*HII_DIM/hlittle);

      if ((catalog_halo.M_vir/hlittle) >= M_MIN)
      // if (catalog_halo.M_vir >= 172e8) // 20*Millennium particle mass
      // if (catalog_halo.M_vir >= 150e9) // 20*GiggleZ_main particle mass
        *(halo_field + HII_R_FFT_INDEX(
              (int)(catalog_halo.position_MBP[0]/BOX_LEN*HII_DIM/hlittle),
              (int)(catalog_halo.position_MBP[1]/BOX_LEN*HII_DIM/hlittle),
              (int)(catalog_halo.position_MBP[2]/BOX_LEN*HII_DIM/hlittle))) 
          += (float)(catalog_halo.M_vir/hlittle);
      
#ifdef DEBUG
      mass[halo_count] = catalog_halo.M_vir/hlittle;
      pos[halo_count][0] = catalog_halo.position_MBP[0]/BOX_LEN*HII_DIM/hlittle;
      pos[halo_count][1] = catalog_halo.position_MBP[1]/BOX_LEN*HII_DIM/hlittle;
      pos[halo_count][2] = catalog_halo.position_MBP[2]/BOX_LEN*HII_DIM/hlittle;
#endif

      halo_count++;
    }
    fclose(fin);
  }

  // Check to make sure we have read the expected number of halos
  if (halo_count != N_halos_total)
  {
    fprintf(stderr, "halo_count (%d) != N_halos_total (%d) in read!...\n", halo_count, N_halos_total);
    return -1;
  } else
    printf("Read %d halos\n", halo_count);

#ifdef DEBUG
  // Turn off error reporting in hdf5
  H5Eset_auto(H5P_DEFAULT, NULL, NULL);

  // Save the halo data
  char name[256];
  hid_t fout, group;
  if((fout = H5Fopen("../Boxes/halos.hdf5", H5F_ACC_RDWR, H5P_DEFAULT))<0)
    fout = H5Fcreate("../Boxes/halos.hdf5", H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(name, "/snap%04d", snapshot);
  if((group = H5Gcreate(fout, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT))<0)
    fprintf(stderr, "Group already exists in ../Boxes/halos.hdf5... Skipping...\n");
  else
  {
    hsize_t dims_mass[1];
    dims_mass[0] = halo_count;
    H5LTmake_dataset(group, "mass", 1, dims_mass, H5T_NATIVE_FLOAT, mass);
    hsize_t dims_pos[2];
    dims_pos[0] = halo_count;
    dims_pos[1] = 3;
    H5LTmake_dataset(group, "pos", 2, dims_pos, H5T_NATIVE_FLOAT, pos);
    H5Gclose(group);
  }
  free(pos);
  free(mass);
  H5Fclose(fout);
#endif

  fprintf(stderr, "...done!\n");
  return 0;
}

