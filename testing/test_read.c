#define _MAIN
#include "trees.h" 
#include "hdf5.h"
#include "hdf5_hl.h"

int main(int argc, char *argv[])
{
  
  // Initialise gbpSID
  SID_init(&argc, &argv, NULL);
  SID.fp_log = stdout;
  
  // Quick argument parsing test
  if(argc!=6){
    printf("Usage:\n\t%s sim total_sim_snaps n_every_snaps n_scan_snaps snapshot\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  int snapshot = atoi(argv[5]);

  Halo *halos;
  TreesHeader header = read_trees(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), snapshot, &halos);
  if (header.n_subgroups==0)
    SID_exit(ERROR_NONE);

  // Dump the results as an hdf5 file
  char   fname[STR_LEN];
  hid_t  file_id, dataspace_id;
  hsize_t dims[1] = {header.n_subgroups};
  hsize_t dims_2d[2] = {header.n_subgroups, 3};
  herr_t status;
  sprintf(fname, "output/snap%03d.hdf5", atoi(argv[5]));
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // fprintf(fout, "# n_groups : %d\n"        , header.n_groups);
  // fprintf(fout, "# n_subgroups : %d\n"     , header.n_subgroups);
  // fprintf(fout, "# n_halos_max : %d\n"     , header.n_halos_max);
  // fprintf(fout, "# n_trees_subgroup : %d\n", header.n_trees_subgroup);
  // fprintf(fout, "# n_trees_group : %d\n"   , header.n_trees_group);
  
  dataspace_id = H5Screate_simple(1, dims, NULL);

  int data_int[header.n_subgroups];
  long long data_ll[header.n_subgroups];
  double data_double[header.n_subgroups];
  float data_float[header.n_subgroups];
  float data_float3[header.n_subgroups][3];
  hid_t ds_type;

  // id
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].id;
  status = H5LTmake_dataset(file_id,"/id",1,dims,ds_type,data_int);

  // type
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].type;
  status = H5LTmake_dataset(file_id,"/type",1,dims,ds_type,data_int);

  // desc_id
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].desc_id;
  status = H5LTmake_dataset(file_id,"/desc_id",1,dims,ds_type,data_int);

  // tree_id
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].tree_id;
  status = H5LTmake_dataset(file_id,"/tree_id",1,dims,ds_type,data_int);

  // file_offset
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].file_offset;
  status = H5LTmake_dataset(file_id,"/file_offset",1,dims,ds_type,data_int);

  // file_index
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].file_index;
  status = H5LTmake_dataset(file_id,"/file_index",1,dims,ds_type,data_int);

  // tree_flags
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].tree_flags;
  status = H5LTmake_dataset(file_id,"/tree_flags",1,dims,ds_type,data_int);

  // n_subgroups
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].n_subgroups;
  status = H5LTmake_dataset(file_id,"/n_subgroups",1,dims,ds_type,data_int);

  // id_MBP
  ds_type = H5T_NATIVE_LLONG;
  for(int i=0; i<header.n_subgroups; i++)
    data_ll[i] = halos[i].id_MBP;
  status = H5LTmake_dataset(file_id,"/id_MBP",1,dims,ds_type,data_ll);

  // M_vir
  ds_type = H5T_NATIVE_DOUBLE;
  for(int i=0; i<header.n_subgroups; i++)
    data_double[i] = halos[i].M_vir;
  status = H5LTmake_dataset(file_id,"/M_vir",1,dims,ds_type,data_double);

  // n_particles
  ds_type = H5T_NATIVE_INT;
  for(int i=0; i<header.n_subgroups; i++)
    data_int[i] = halos[i].n_particles;
  status = H5LTmake_dataset(file_id,"/n_particles",1,dims,ds_type,data_int);

  // position_COM
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    for(int j=0; j<3; j++)
      data_float3[i][j] = halos[i].position_COM[j];
  status = H5LTmake_dataset(file_id,"/position_COM",2,dims_2d,ds_type,data_float3);

  // position_MBP
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    for(int j=0; j<3; j++)
      data_float3[i][j] = halos[i].position_MBP[j];
  status = H5LTmake_dataset(file_id,"/position_MBP",2,dims_2d,ds_type,data_float3);

  // velocity_COM
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    for(int j=0; j<3; j++)
      data_float3[i][j] = halos[i].velocity_COM[j];
  status = H5LTmake_dataset(file_id,"/velocity_COM",2,dims_2d,ds_type,data_float3);

  // velocity_MBP
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    for(int j=0; j<3; j++)
      data_float3[i][j] = halos[i].velocity_MBP[j];
  status = H5LTmake_dataset(file_id,"/velocity_MBP",2,dims_2d,ds_type,data_float3);

  // R_vir
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].R_vir;
  status = H5LTmake_dataset(file_id,"/R_vir",1,dims,ds_type,data_float);

  // R_halo
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].R_halo;
  status = H5LTmake_dataset(file_id,"/R_halo",1,dims,ds_type,data_float);

  // R_max
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].R_max;
  status = H5LTmake_dataset(file_id,"/R_max",1,dims,ds_type,data_float);

  // V_max
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].V_max;
  status = H5LTmake_dataset(file_id,"/V_max",1,dims,ds_type,data_float);

  // sigma_v
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].sigma_v;
  status = H5LTmake_dataset(file_id,"/sigma_v",1,dims,ds_type,data_float);

  // spin
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    for(int j=0; j<3; j++)
      data_float3[i][j] = halos[i].spin[j];
  status = H5LTmake_dataset(file_id,"/spin",2,dims_2d,ds_type,data_float3);

  // q_triaxial
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].q_triaxial;
  status = H5LTmake_dataset(file_id,"/q_triaxial",1,dims,ds_type,data_float);

  // s_triaxial
  ds_type = H5T_NATIVE_FLOAT;
  for(int i=0; i<header.n_subgroups; i++)
    data_float[i] = halos[i].s_triaxial;
  status = H5LTmake_dataset(file_id,"/s_triaxial",1,dims,ds_type,data_float);

  // clean up
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);


  // fprintf(fout, "index    type    id    desc_id    tree_id    file_offset    file_index    n_subgroups    M_vir    n_particles    V_max    R_vir    R_halo\n");

  // Halo *this_halo;
  // for (int i=0; i<header.n_subgroups; i++){
  //   this_halo = &(halos[i]);
  //   fprintf(fout, "%04d    %d    %d    %d    %d    %d    %d    %d    %.3e    %d    %.3e    %.3e    %.3e\n", 
  //       i, this_halo->type, this_halo->id,
  //       this_halo->desc_id, this_halo->tree_id, this_halo->file_offset,
  //       this_halo->file_index, this_halo->n_subgroups, this_halo->M_vir,
  //       this_halo->n_particles, this_halo->V_max, this_halo->R_vir,
  //       this_halo->R_halo);
  // }

  // fclose(fout);

  // Free allocated arrays
  free_trees(&halos);

  SID_exit(ERROR_NONE);

}

