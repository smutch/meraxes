#ifdef USE_TOCF

#include "meraxes.h"
#include <math.h>

static inline void read_identifier(FILE *fin, bool skip_flag)
{
  char identifier[32];

  fread(identifier, sizeof(identifier), 1, fin);
  if (skip_flag)
    SID_log_error("Skipping grid: %s...", identifier);
  else
    SID_log("Reading grid: %s...", SID_LOG_COMMENT, identifier);
}


int read_dm_grid(
  run_globals_t *run_globals,
  int            snapshot,
  int            i_grid,
  float         *grid_out)
{
  // N.B. We assume in this function that the grid has the fftw3 inplace
  // complex dft padding.

  char fname[512];
  FILE *fin;
  int n_cell[3];
  double box_size[3];
  int n_grids;
  int ma_scheme;
  int n_elem;
  float val;
  double mean              = 0.;
  double cell_volume       = 0.;
  double cell_volume_ratio = 0.;
  double resample_factor   = 1.;
  run_params_t *params  = &(run_globals->params);
  int HII_dim           = tocf_params.HII_dim;
  double *grid;

  // Construct the input filename
  sprintf(fname, "%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, snapshot);
  // ... and open
  fin = fopen(fname, "rb");
  if (!fin)
  {
    SID_log_error("Failed to open file: %s", fname);
    return(EXIT_FAILURE);
  }

  // Read the header
  fread(n_cell, sizeof(int), 3, fin);
  fread(box_size, sizeof(double), 3, fin);
  fread(&n_grids, sizeof(int), 1, fin);
  fread(&ma_scheme, sizeof(int), 1, fin);

  SID_log("Reading grid for snapshot %d", SID_LOG_OPEN, snapshot);
  SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, n_cell[0], n_cell[1], n_cell[2]);
  SID_log("box_size = [%.2f, %.2f, %.2f]", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);
  SID_log("ma_scheme = %d", SID_LOG_COMMENT, ma_scheme);

  if (n_grids != 4)
  {
    SID_log_error("n_grids != 4 as expected...");
    fclose(fin);
    return -1;
  }

  // Check if the grid in the file is higher resolution than we require
  if ((n_cell[0] != HII_dim) || (n_cell[1] != HII_dim) || (n_cell[2] != HII_dim))
  {
    resample_factor = (double)HII_dim / (double)n_cell[0];
    if (resample_factor > 1.0001)
    {
      SID_log_error("The dark matter density grid in this file has a resolution less than that required! Aborting!");
      fclose(fin);
      ABORT(EXIT_FAILURE);
    }
    SID_log("Using resample factor = %.3f", SID_LOG_COMMENT, resample_factor);
  }
  else
    resample_factor = 1;

  // Compute the total number of elements in each grid
  n_elem = n_cell[0] * n_cell[1] * n_cell[2];

  // Read the grids
  // Note that we are expecting them to be in a particular order here
  for (int ii = 0; ii < i_grid; ii++)
  {
    read_identifier(fin, true);
    fseek(fin, sizeof(float) * n_elem, SEEK_CUR);
  }
  read_identifier(fin, false);

  // malloc the grid
  grid = SID_calloc(sizeof(double) * HII_dim*HII_dim*HII_dim);

  if (i_grid == 0)  // density grid
  {
    // calculate the volume of a single cell
    cell_volume = pow(box_size[0] / (double)n_cell[0], 3);
    SID_log("cell_volume = %g", SID_LOG_COMMENT, cell_volume);

    // Read in the grid
    for (int i = 0; i < n_cell[0]; i++)
      for (int j = 0; j < n_cell[1]; j++)
        for (int k = 0; k < n_cell[2]; k++)
        {
          fread(&val, sizeof(float), 1, fin);
          mean += (double)val;
          *(grid + HII_R_FFT_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += (double)val;
        }

    SID_log("Sum of density values in high res grid = %g", SID_LOG_COMMENT, mean);
    SID_log("Total mass in high res grid = %g", SID_LOG_COMMENT, mean/cell_volume);

    double grid_sum = 0.0;
    for (int i = 0; i < HII_dim; i++)
      for (int j = 0; j < HII_dim; j++)
        for (int k = 0; k < HII_dim; k++)
          grid_sum += *(grid + HII_R_FFT_INDEX(i, j, k));
    SID_log("After reading grid sum = %g", SID_LOG_COMMENT, grid_sum);

    // In order to calculate the mean density we must actually take the mean of
    // the mass in each cell and then renormalise by the total volume.
    mean *= cell_volume / pow(box_size[0], 3);
    SID_log("Calculated mean = %.2e", SID_LOG_COMMENT, mean);

    // Loop through again and calculate the overdensity
    // i.e. (rho - rho_mean)/rho_mean
    grid_sum = 0.0;
    cell_volume_ratio = pow(box_size[0] / (double)HII_dim, 3) / cell_volume;
    for (int i = 0; i < HII_dim; i++)
      for (int j = 0; j < HII_dim; j++)
        for (int k = 0; k < HII_dim; k++)
        {
          *(grid + HII_R_FFT_INDEX(i, j, k)) = (*(grid + HII_R_FFT_INDEX(i, j, k)) / (cell_volume_ratio * mean)) - 1.;
          grid_sum += *(grid + HII_R_FFT_INDEX(i, j, k));
        }
    SID_log("Total delta sum of grid = %g", SID_LOG_COMMENT, grid_sum);
  }
  else // velocity component grid
  {
    // Read in the grid
    for (int i = 0; i < n_cell[0]; i++)
      for (int j = 0; j < n_cell[1]; j++)
        for (int k = 0; k < n_cell[2]; k++)
        {
          fread(&val, sizeof(float), 1, fin);
          *(grid + HII_R_FFT_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += (double)val;
        }
  }

  // copy the grid (double) to the output (float) and free
  for (int i = 0; i < HII_dim; i++)
    for (int j = 0; j < HII_dim; j++)
      for (int k = 0; k < HII_dim; k++)
        *(grid_out + HII_R_FFT_INDEX(i, j, k)) = (float)*(grid + HII_R_FFT_INDEX(i, j, k));

  SID_free(SID_FARG grid);

  float fgrid_sum = 0.0;
  for (int i = 0; i < HII_dim; i++)
    for (int j = 0; j < HII_dim; j++)
      for (int k = 0; k < HII_dim; k++)
        fgrid_sum += *(grid_out + HII_R_FFT_INDEX(i, j, k));
  SID_log("Total delta sum of FLOAT grid = %g", SID_LOG_COMMENT, fgrid_sum);

  SID_log("...done", SID_LOG_CLOSE);

  // Close the file
  fclose(fin);

  return 0;
}
#endif





// Original with play

//#ifdef USE_TOCF
//
//#include "meraxes.h"
//#include <math.h>
//
//static inline void read_identifier(FILE *fin, bool skip_flag)
//{
//  char identifier[32];
//
//  fread(identifier, sizeof(identifier), 1, fin);
//  if (skip_flag)
//    SID_log_error("Skipping grid: %s...", identifier);
//  else
//    SID_log("Reading grid: %s...", SID_LOG_COMMENT, identifier);
//}
//
//int read_dm_grid(
//  run_globals_t *run_globals,
//  int            snapshot,
//  int            i_grid,
//  float         *grid)
//{
//  // N.B. We assume in this function that the grid has the fftw3 inplace
//  // complex dft padding.
//
//  char fname[512];
//  FILE *fin;
//  int n_cell[3];
//  double box_size[3];
//  int n_grids;
//  int ma_scheme;
//  int n_elem;
//  float val;
//  float mean            = 0.;
//  float cell_volume     = 0.;
//  float resample_factor = 1.;
//  run_params_t *params  = &(run_globals->params);
//  int HII_dim           = tocf_params.HII_dim;
//
//  // Construct the input filename
//  sprintf(fname, "%s/grids/snapshot_%03d_dark_grid.dat", params->SimulationDir, snapshot);
//  // ... and open
//  fin = fopen(fname, "rb");
//  if (!fin)
//  {
//    SID_log_error("Failed to open file: %s", fname);
//    return(EXIT_FAILURE);
//  }
//
//  // Read the header
//  fread(n_cell, sizeof(int), 3, fin);
//  fread(box_size, sizeof(double), 3, fin);
//  fread(&n_grids, sizeof(int), 1, fin);
//  fread(&ma_scheme, sizeof(int), 1, fin);
//
//  SID_log("Reading grid for snapshot %d", SID_LOG_OPEN, snapshot);
//  SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, n_cell[0], n_cell[1], n_cell[2]);
//  SID_log("box_size = [%.2f, %.2f, %.2f]", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);
//  SID_log("ma_scheme = %d", SID_LOG_COMMENT, ma_scheme);
//
//  if (n_grids != 4)
//  {
//    SID_log_error("n_grids != 4 as expected...");
//    fclose(fin);
//    return -1;
//  }
//
//  // Check if the grid in the file is higher resolution than we require
//  if ((n_cell[0] != HII_dim) || (n_cell[1] != HII_dim) || (n_cell[2] != HII_dim))
//  {
//    resample_factor = (float)HII_dim / (float)n_cell[0];
//    if (resample_factor > 1.0001)
//    {
//      SID_log_error("The dark matter density grid in this file has a resolution less than that required! Aborting!");
//      fclose(fin);
//      ABORT(EXIT_FAILURE);
//    }
//    SID_log("Using resample factor = %.3f", SID_LOG_COMMENT, resample_factor);
//  }
//  else
//    resample_factor = 1;
//
//  // Compute the total number of elements in each grid
//  n_elem = n_cell[0] * n_cell[1] * n_cell[2];
//
//  // Read the grids
//  // Note that we are expecting them to be in a particular order here
//  for (int ii = 0; ii < i_grid; ii++)
//  {
//    read_identifier(fin, true);
//    fseek(fin, sizeof(float) * n_elem, SEEK_CUR);
//  }
//  read_identifier(fin, false);
//
//  // init the grid
//  for (int i = 0; i < HII_dim; i++)
//    for (int j = 0; j < HII_dim; j++)
//      for (int k = 0; k < HII_dim; k++)
//        *(grid + HII_R_FFT_INDEX(i, j, k)) = 0.;
//
//  if (i_grid == 0)  // density grid
//  {
//    // calculate the volume of a single cell
//    cell_volume = powf(box_size[0] / (float)n_cell[0], 3);
//    
//    SID_log("HIGH RES INFO:", SID_LOG_OPEN);   // PMG
//    SID_log("box_size = %g", SID_LOG_COMMENT, box_size[0]);   // PMG
//    SID_log("n_cell = %d", SID_LOG_COMMENT, n_cell[0]);   // PMG
//    SID_log("cell_volume = %g", SID_LOG_COMMENT, cell_volume);   // PMG
//    
//    
//    float sub_sum = 0.0;   // PMG
//    
//    // Read in the grid
//    for (int i = 0; i < n_cell[0]; i++)
//      for (int j = 0; j < n_cell[1]; j++)
//        for (int k = 0; k < n_cell[2]; k++)
//        {
//          fread(&val, sizeof(float), 1, fin);
//          val *= cell_volume; // now we have the mass in the cell
//          
//          if(i>=4 && i<=7 && j<=3 && k>=4 && k<=7)
//          {
//            sub_sum += val;   // PMG
////            SID_log("HR grid mass (%d, %d, %d) = %g", SID_LOG_COMMENT, i, j, k, val);   // PMG
//          }
//          
//          mean += val;
//          *(grid + HII_R_FFT_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += val;
//        }
//    
//    SID_log("sub_sum corre to (1, 0, 1) = %g", SID_LOG_COMMENT, sub_sum );   // PMG
//    
//    // This is just a counting check
////    for (int i = 0; i < 100; i++)   // PMG
////    SID_log("i = %3d \t (int)(i * resample_factor) = %d", SID_LOG_COMMENT, i, (int)(i * resample_factor));   // PMG
//        
//    SID_log("LR grid mass (1, 0, 1) = %g", SID_LOG_COMMENT, (*(grid + HII_R_FFT_INDEX(1, 0, 1))) );   // PMG
//    
//    
//    // This mean is the total mass in the high resolution box
//    SID_log("total mass HR = %g", SID_LOG_COMMENT, mean);   // PMG
//    
//    FILE *pmg_f1;   // PMG
//    char pmg_string1[100];   // PMG
//    sprintf(pmg_string1, "output-UVFB_OFF_HEF-5_Snaps-1-99_RERUN/grid_HR_mass_CHECK_snap%d.dat", snapshot);   // PMG
//    pmg_f1 = fopen(pmg_string1, "wt");   // PMG
//    fprintf(pmg_f1, "%g\n", mean);   // PMG
//    fclose(pmg_f1);   // PMG
//    
//    
//    
//    // In order to calculate the mean density we must actually take the mean of
//    // the mass in each cell and then renormalise by the total volume.
//    mean /= powf(box_size[0], 3);
//    //SID_log("Calculated mean = %.2e :: theory = %.2e", SID_LOG_COMMENT, mean*1.e10*hlittle*hlittle, OMm*RHOcrit);
//    
//    // This mean is the mean density of the high resolution box
//    SID_log("mean density HR = %g", SID_LOG_COMMENT, mean);   // PMG
//    SID_log("...done HR", SID_LOG_CLOSE);
//    
//    
//    // At this stage grid holds the total mass of each LR cell
//    SID_log("LR grid mass (1, 0, 1) = %g", SID_LOG_COMMENT, (*(grid + HII_R_FFT_INDEX(1, 0, 1))) );   // PMG
//    float LR_total_mass = 0.0;   // PMG
//    for (int i = 0; i < HII_dim; i++)   // PMG
//      for (int j = 0; j < HII_dim; j++)   // PMG
//        for (int k = 0; k < HII_dim; k++)   // PMG
//          LR_total_mass += (*(grid + HII_R_FFT_INDEX(i, j, k)));   // PMG
//    
//    
//    SID_log("LOW RES INFO:", SID_LOG_OPEN);   // PMG
//    SID_log("total mass LR = %g", SID_LOG_COMMENT, LR_total_mass);   // PMG
//    SID_log("...done LR", SID_LOG_CLOSE);
//
//    
//    
//    SID_log("HIGH RES INFO2:", SID_LOG_OPEN);   // PMG
//    SID_log("total mass HR = %g", SID_LOG_COMMENT, mean*powf(box_size[0], 3));   // PMG
//    SID_log("...done HR2", SID_LOG_CLOSE);
//    
//    
////    KEEP for later ! if(snapshot==9 || snapshot==17 || snapshot==26 || snapshot==35 || snapshot==42 || snapshot==50 || snapshot==62 || snapshot==77 || snapshot==99)
//    
//    
//    // Loop through again and calculate the overdensity
//    // i.e. (rho - rho_mean)/rho_mean
//    cell_volume = powf(box_size[0] / (float)HII_dim, 3);
//    for (int i = 0; i < HII_dim; i++)
//      for (int j = 0; j < HII_dim; j++)
//        for (int k = 0; k < HII_dim; k++)
//          *(grid + HII_R_FFT_INDEX(i, j, k)) = (*(grid + HII_R_FFT_INDEX(i, j, k)) / (cell_volume * mean)) - 1.;
//    
//    
//    
//    float total_lo_res_mass = 0.0;   // PMG
//    for (int i = 0; i < HII_dim; i++)   // PMG
//        for (int j = 0; j < HII_dim; j++)   // PMG
//            for (int k = 0; k < HII_dim; k++)   // PMG
//                total_lo_res_mass += (   (1.0 + *(grid + HII_R_FFT_INDEX(i, j, k)))*mean*cell_volume   );   // PMG
//    
//    SID_log("LOW RES INFO2:", SID_LOG_OPEN);   // PMG
//    SID_log("total mass LR = %g", SID_LOG_COMMENT, total_lo_res_mass);   // PMG
//    SID_log("...done LR2", SID_LOG_CLOSE);
//    
//    
//    
//    
//    // START - temp write file for checking
////    if(snapshot==9 || snapshot==17 || snapshot==26 || snapshot==35 || snapshot==42 || snapshot==50 || snapshot==62 || snapshot==77 || snapshot==99)
////    {
////        float total_lo_res_mass = 0.0;
////        
////        for (int i = 0; i < HII_dim; i++)
////            for (int j = 0; j < HII_dim; j++)
////                for (int k = 0; k < HII_dim; k++)
////                    total_lo_res_mass += (   (1.0 + *(grid + HII_R_FFT_INDEX(i, j, k)))*mean*cell_volume   );
////        
////        FILE *pmg_f2;
////        char pmg_string2[100];
////        
////        sprintf(pmg_string2, "output-UVFB_OFF_HEF-5_Snaps-1-99_RERUN/grid_LR_mass_CHECK_snap%d.dat", snapshot);
////        pmg_f2 = fopen(pmg_string2, "wt");
////        
////        fprintf(pmg_f2, "SNAPSHOT \t %d\n\n", snapshot);
////        fprintf(pmg_f2, "mean_density \t %g\n", mean);
////        fprintf(pmg_f2, "box_size \t %g\n", box_size[0]);
////        fprintf(pmg_f2, "n_cell_lo_res \t %d\n", HII_dim);
////        fprintf(pmg_f2, "cell_volume_lo_res \t %g\n", cell_volume);
////        fprintf(pmg_f2, "total_lo_res_mass \t %g\n", total_lo_res_mass);
////        
////        fclose(pmg_f2);
//        
//        
////        sprintf(pmg_string2, "output-UVFB_OFF_HEF-5_Snaps-1-99_RERUN/deltax_grid_CHECK_snap%d.dat", snapshot);
////    	pmg_f2 = fopen(pmg_string2, "wt");
////        
////        for (int i = 0; i < HII_dim; i++)
////          for (int j = 0; j < HII_dim; j++)
////            for (int k = 0; k < HII_dim - 1; k++)
////                fprintf(pmg_f2, "%g\n", (float)(*(grid + HII_R_FFT_INDEX(i, j, k)))  );
////        
////    	fclose(pmg_f2);
////    	
////    	sprintf(pmg_string2, "output-UVFB_OFF_HEF-5_Snaps-1-99_RERUN/DDeltax_grid_CHECK_snap%d.dat", snapshot);
////    	pmg_f2 = fopen(pmg_string2, "wt");
////        
////        for (int i = 0; i < HII_dim; i++)
////          for (int j = 0; j < HII_dim; j++)
////            for (int k = 0; k < HII_dim - 1; k++)
////                fprintf( pmg_f2, "%g\n", 1.0 + (float)(*(grid + HII_R_FFT_INDEX(i, j, k))) );
////        
////    	fclose(pmg_f2);
////    }
//    // END - temp write file for checking
//    
//    
//  }
//  else // velocity component grid
//  {
//    // Read in the grid
//    for (int i = 0; i < n_cell[0]; i++)
//      for (int j = 0; j < n_cell[1]; j++)
//        for (int k = 0; k < n_cell[2]; k++)
//        {
//          fread(&val, sizeof(float), 1, fin);
//          *(grid + HII_R_FFT_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += val;
//        }
//  }
//
//  SID_log("...done", SID_LOG_CLOSE);
//
//  // Close the file
//  fclose(fin);
//
//  return 0;
//}
//#endif
