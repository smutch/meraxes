#ifdef USE_TOCF

#include "meraxes.h"
#include <math.h>


/*
  ==============================================================================
  MAJOR CODE REVISION by Paul Geil (October 2014)
  ==============================================================================
  
  - A significant numerical accuracy bug was resolved by performing calculations
    using an array of doubles and then casting it as an array of floats (as
    required by 21cmfast)
    
  - Snapshots 53, 57, 61, 65 and 69 in TIAMAT dm grids contain anomalies [single
    highly oversense voxels at (0, 0, 0) and extensive slabs with zero density]
    which significantly affect their statistics. The density spikes have been
    remedied by resetting the density of the offending voxel to that of the
    average over its (6) neigbours. Anomalous zero density regions are reassigned
    to be the mean density of the (non-empty) grid.
    
  - The TIAMAT velocity grids have not been analysed for anomalies.
*/



static inline void read_identifier(FILE *fin, bool skip_flag)
{
    char identifier[32];
    
    fread(identifier, sizeof(identifier), 1, fin);
    
    if (skip_flag)
        SID_log_error("Skipping grid: %s...", identifier);
    else
        SID_log("Reading grid: %s...", SID_LOG_COMMENT, identifier);
}



static unsigned long long HR_INDEX(int i, int j, int k, int grid_dim)
{ 
    return (unsigned long long)(k + grid_dim*(j + grid_dim*i));
}



static long long find_HR_empty_count(double *grid_HR, int HR_dim_val)
{
    long long count_val = 0;
    
    for (int i = 0; i < HR_dim_val; i++)
        for (int j = 0; j < HR_dim_val; j++)
            for (int k = 0; k < HR_dim_val; k++)
            {
                if (*((double *)grid_HR + HR_INDEX(i, j, k, HR_dim_val)) == 0.0) count_val++;
                
            }
    
    return (long long)(count_val);
}

static double find_HR_non_empty_mean(double *grid_HR, int HR_dim_val, long long N)
{
    double val;
    double mean_val = 0.0;
    
    for (int i = 0; i < HR_dim_val; i++)
        for (int j = 0; j < HR_dim_val; j++)
            for (int k = 0; k < HR_dim_val; k++)
            {
                val = *((double *)grid_HR + HR_INDEX(i, j, k, HR_dim_val));
                
                if (val > 0.0) mean_val += val;
                
            }
    
    mean_val /= (double)(N);
    
    return mean_val;
}






int read_dm_grid(
    run_globals_t *run_globals,
    int            snapshot,
    int            i_grid,
    float         *grid_out)
{
    // N.B. We assume in this function that the grid has the fftw3 inplace complex dft padding.
    
    char   fname[512];
    FILE   *fin;
    int    n_cell[3];
    double box_size[3];
    int    n_grids;
    int    ma_scheme;
    int    n_elem;
    float  val;
    double mean              = 0.;
    double cell_volume       = 0.;
    double cell_volume_ratio = 0.;
    double resample_factor   = 1.;
    int    HII_dim           = tocf_params.HII_dim;
    double *grid;
    double *grid_HR;
    int    HR_dim;
    
    run_params_t *params = &(run_globals->params);
    
    
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
    
    SID_log("Reading grid for snapshot %d", SID_LOG_OPEN | SID_LOG_TIMER, snapshot);
    SID_log("n_cell = [%d, %d, %d]", SID_LOG_COMMENT, n_cell[0], n_cell[1], n_cell[2]);
    SID_log("box_size = [%.2f, %.2f, %.2f]", SID_LOG_COMMENT, box_size[0], box_size[1], box_size[2]);
    SID_log("ma_scheme = %d", SID_LOG_COMMENT, ma_scheme);
    
    // Assuming the grid is cubic!
    HR_dim = n_cell[0];
    
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
    grid    = SID_calloc(sizeof(double) * HII_dim*HII_dim*HII_dim);
    grid_HR = SID_calloc(sizeof(double) * n_elem);
    
    // Initialise the grids (just in case!)
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
                *(grid + HII_R_INDEX(i, j, k)) = 0.0;
    
    for (int i = 0; i < n_cell[0]; i++)
        for (int j = 0; j < n_cell[1]; j++)
            for (int k = 0; k < n_cell[2]; k++)
                *(grid_HR + HR_INDEX(i, j, k, HR_dim)) = 0.0;
    
    
    if (i_grid == 0)  // Density grid
    {
        // Read in the grid and assign grid_HR values
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    fread(&val, sizeof(float), 1, fin);
                    *(grid_HR + HR_INDEX(i, j, k, HR_dim)) = (double)val;
                }
        
        
        // ************ START OF QUICK FIXES FOR TIAMAT !!!
        
        long long empty_count;
        long long non_empty_count;
        double    non_empty_mean;
        
        // For these problematic snapshots only
        
        if (snapshot==53 || snapshot==57 || snapshot==61 || snapshot==65 || snapshot==69)
        {
            // From previous analysis, we know that the offending maximum spike voxel for these snapshots is (0, 0, 0)
            // Now reset its value to that of the average over its neigbours
            
            SID_log("Revaluating (0,0,0) dm density voxel:", SID_LOG_OPEN);
            
            *(grid_HR + HR_INDEX(0, 0, 0, HR_dim)) = ( *(grid_HR + HR_INDEX(1, 0, 0, HR_dim)) +
                                                        *(grid_HR + HR_INDEX(0, 1, 0, HR_dim)) +
                                                        *(grid_HR + HR_INDEX(0, 0, 1, HR_dim)) +
                                                        *(grid_HR + HR_INDEX(HR_dim - 1, 0, 0, HR_dim)) +
                                                        *(grid_HR + HR_INDEX(0, HR_dim - 1, 0, HR_dim)) +
                                                        *(grid_HR + HR_INDEX(0, 0, HR_dim - 1, HR_dim))) / 6.0;
            SID_log("...done", SID_LOG_CLOSE);
            
            // From previous analysis, we know that there are anomalous empty regions
            // Now reassign these voxels with the mean of the (non-empty) grid
            
            SID_log("Reassigning empty dm density voxels:", SID_LOG_OPEN);
            
            empty_count = find_HR_empty_count(grid_HR, HR_dim);
            non_empty_count = n_elem - empty_count;
            non_empty_mean = find_HR_non_empty_mean(grid_HR, HR_dim, non_empty_count);
            
            for (int i = 0; i < HR_dim; i++)
                for (int j = 0; j < HR_dim; j++)
                    for (int k = 0; k < HR_dim; k++)
                    {
                        if (*(grid_HR + HR_INDEX(i, j, k, HR_dim)) == 0.0) *(grid_HR + HR_INDEX(i, j, k, HR_dim)) = non_empty_mean;
                        
                    }
            
            SID_log("...done", SID_LOG_CLOSE);
        }
        
        // ************ END OF QUICK FIXES FOR TIAMAT !!!
        
        
        // Regrid
        mean = 0.0;
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    mean += *(grid_HR + HR_INDEX(i, j, k, HR_dim));
                    *(grid + HII_R_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += *(grid_HR + HR_INDEX(i, j, k, HR_dim));
                }
        
        // Calculate the volume of a single high resolution cell
        cell_volume = pow(box_size[0] / (double)n_cell[0], 3);
        
        // Mean density from high res grid
        mean *= cell_volume / pow(box_size[0], 3);
        
        // At this point grid holds the summed densities in each LR cell
        // Loop through again and calculate the overdensity
        // i.e. (rho - rho_mean)/rho_mean
        cell_volume_ratio = pow(box_size[0] / (double)HII_dim, 3) / cell_volume;
        for (int i = 0; i < HII_dim; i++)
            for (int j = 0; j < HII_dim; j++)
                for (int k = 0; k < HII_dim; k++)
                {
                    *(grid + HII_R_INDEX(i, j, k)) = (*(grid + HII_R_INDEX(i, j, k)) / (cell_volume_ratio * mean)) - 1.;
                }
        
    }
    else // Velocity component grid
    {
        // Read in the grid
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    fread(&val, sizeof(float), 1, fin);
                    *(grid + HII_R_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += (double)val;
                }
    }
    
    // Copy the grid (double) to the output (float) and free
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
                *(grid_out + HII_R_FFT_INDEX(i, j, k)) = (float)*(grid + HII_R_INDEX(i, j, k));
    
    
    SID_free(SID_FARG grid);
    SID_free(SID_FARG grid_HR);
    
    SID_log("...done", SID_LOG_CLOSE);
    
    // Close the file
    fclose(fin);
    
    return 0;
}
#endif
