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
    // N.B. We assume in this function that the grid has the fftw3 inplace complex dft padding.
    
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
    {
        resample_factor = 1;
    }
    
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
    
    // Malloc the grid
    grid = SID_calloc(sizeof(double) * HII_dim*HII_dim*HII_dim);
    
    // Initialise the grid
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
                *(grid + HII_R_FFT_INDEX(i, j, k)) = 0.0;
    

    if (i_grid == 0)  // density grid
    {
        // Calculate the volume of a single high resolution cell
        cell_volume = pow(box_size[0] / (double)n_cell[0], 3);
        
//        SID_log("HIGH RES INFO:", SID_LOG_OPEN);   // PMG
//        SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, n_cell[0], n_cell[1], n_cell[2]);   // PMG
//        SID_log("n_elem = %d", SID_LOG_COMMENT, n_elem);   // PMG
//        SID_log("box_size = [%.2f, %.2f, %.2f]", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);   // PMG
//        SID_log("cell_volume = %e", SID_LOG_COMMENT, cell_volume);   // PMG
//        
//        double my_HR_mass = 0.0;   // PMG
//        double my_min_density = 1.0e10;   // PMG
//        double my_max_density = -1.0e10;   // PMG
//        double my_mean_density = 0.0;   // PMG
        
        mean = 0.0;   /// PMG just to be sure to be sure
        
        // Read in the grid
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    fread(&val, sizeof(float), 1, fin);
                    mean += (double)val;
                    
//                    my_mean_density += (double)val;   // PMG
//                    
//                    if ((double)val<my_min_density) my_min_density = (double)val;   // PMG
//                    
//                    if ((double)val>my_max_density) my_max_density = (double)val;   // PMG
//                    
//                    my_HR_mass += (double)(val)*cell_volume;   // PMG mass = density * volume
//                    
//                    if (i == 32 && j == 32 && k == 32)   // PMG
//                    {
//                        SID_log("val(32,32,32) = %g", SID_LOG_COMMENT, val);   // PMG
//                        SID_log("mass(32,32,32) = %g", SID_LOG_COMMENT, (double)(val)*cell_volume);   // PMG
//                    }
                    
                    *(grid + HII_R_FFT_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += (double)val;
                }
        
//        my_mean_density /= n_elem;   // PMG
//        
//        SID_log("my_min_density = %g", SID_LOG_COMMENT, my_min_density);   // PMG
//        SID_log("      min mass = %g", SID_LOG_COMMENT, my_min_density*cell_volume);   // PMG
//        SID_log("my_max_density = %g", SID_LOG_COMMENT, my_max_density);   // PMG
//        SID_log("      max mass = %g", SID_LOG_COMMENT, my_max_density*cell_volume);   // PMG
//        
//        SID_log("mean density    = %g", SID_LOG_COMMENT, mean/n_elem);   // PMG
//        SID_log("my_mean_density = %g", SID_LOG_COMMENT, my_mean_density);   // PMG
//        SID_log("      mean mass = %g", SID_LOG_COMMENT, my_mean_density*cell_volume);   // PMG
//        
//        SID_log("Sum of density values in HR grid = %g", SID_LOG_COMMENT, mean);   // SM
//        SID_log("Total mass in HR grid = %g", SID_LOG_COMMENT, mean*cell_volume);   // PMG
//        SID_log("cf         my_HR_mass = %g", SID_LOG_COMMENT, my_HR_mass);   // PMG
        
        // Mean density from high res grid
        mean *= cell_volume / pow(box_size[0], 3);
//        SID_log("Calculated mean density from HR grid = %g", SID_LOG_COMMENT, mean);   // SM
//        SID_log("...done HR", SID_LOG_CLOSE);   // PMG
//        
//        SID_log("LOW RES INFO:", SID_LOG_OPEN);   // PMG
//        SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, HII_dim, HII_dim, HII_dim);   // PMG
//        SID_log("box_size = [%.2f, %.2f, %.2f]", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);   // PMG
        
        // At this point grid holds the summed densities in each LR cell
//        double grid_sum = 0.0;
//        for (int i = 0; i < HII_dim; i++)
//            for (int j = 0; j < HII_dim; j++)
//                for (int k = 0; k < HII_dim; k++)
//                    grid_sum += (double)( *(grid + HII_R_FFT_INDEX(i, j, k)) );
//        
//        SID_log("Sum of density values in LR grid = %g", SID_LOG_COMMENT, grid_sum);   // PMG
        
        // At this point grid holds the summed densities in each LR cell
        // Loop through again and calculate the overdensity
        // i.e. (rho - rho_mean)/rho_mean
//        double grid_sum = 0.0;
        cell_volume_ratio = pow(box_size[0] / (double)HII_dim, 3) / cell_volume;
        for (int i = 0; i < HII_dim; i++)
            for (int j = 0; j < HII_dim; j++)
                for (int k = 0; k < HII_dim; k++)
                {
                    *(grid + HII_R_FFT_INDEX(i, j, k)) = (*(grid + HII_R_FFT_INDEX(i, j, k)) / (cell_volume_ratio * mean)) - 1.;
//                    grid_sum += *(grid + HII_R_FFT_INDEX(i, j, k));
                }
        
//        SID_log("Total delta sum of LR DOUBLE grid = %g", SID_LOG_COMMENT, grid_sum);
//        float Total_mass_in_LR_grid = mean*pow(box_size[0]/(double)HII_dim,3)*(grid_sum + (float)(HII_dim*HII_dim*HII_dim));   // PMG
//        SID_log("Total mass in LR grid = %g", SID_LOG_COMMENT, Total_mass_in_LR_grid);   // PMG
//        SID_log("...done LR", SID_LOG_CLOSE);   // PMG
    }
    else // Velocity component grid
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
    
    // Copy the grid (double) to the output (float) and free
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
                *(grid_out + HII_R_FFT_INDEX(i, j, k)) = (float)*(grid + HII_R_FFT_INDEX(i, j, k));
    
    
    // PMG SECTION START
    
    float delta_LRF;
    float delta_LRF_min = 1.0e10;
    float delta_LRF_max = -1.0e10;
    float delta_LRF_sum = 0.0;
    float delta_LRF_mean;
    float delta_LRF_var = 0.0;
    float n_elem_LR = (float)(HII_dim*HII_dim*HII_dim);
    
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
            {
                delta_LRF = *(grid_out + HII_R_FFT_INDEX(i, j, k));
                
                delta_LRF_sum += delta_LRF;
                
                if (delta_LRF<delta_LRF_min) delta_LRF_min = delta_LRF;
                    
                if (delta_LRF>delta_LRF_max) delta_LRF_max = delta_LRF;
                
            }
    
    delta_LRF_mean = delta_LRF_sum/n_elem_LR;
    
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
            {
                delta_LRF = *(grid_out + HII_R_FFT_INDEX(i, j, k));
                delta_LRF_var += (delta_LRF - delta_LRF_mean)*(delta_LRF - delta_LRF_mean);
            }
    
    delta_LRF_var /= ((float)(n_elem_LR - 1.0));
    
    SID_log("delta_LRF_sum  = %g", SID_LOG_COMMENT, delta_LRF_sum);
    SID_log("delta_LRF_min  = %g", SID_LOG_COMMENT, delta_LRF_min);
    SID_log("delta_LRF_max  = %g", SID_LOG_COMMENT, delta_LRF_max);
    SID_log("delta_LRF_mean = %g", SID_LOG_COMMENT, delta_LRF_mean);
    SID_log("delta_LRF_var  = %g", SID_LOG_COMMENT, delta_LRF_var);
    
    // PMG SECTION END
    
    
    
    SID_free(SID_FARG grid);
    
    SID_log("...done", SID_LOG_CLOSE);
    
    // Close the file
    fclose(fin);
    
    return 0;
}
#endif

