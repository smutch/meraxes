#ifdef USE_TOCF

#include "meraxes.h"
#include <math.h>


/*
  ==============================================================================
  MAJOR CODE REVISION by Paul Geil (Octobner 2014)
  ==============================================================================
  
  - A significant numerical accuracy bug was resolved by performing calculations
    using an array of doubles and then casting it as an array of floats (as
    required by 21cmfast)
    
  - Snapshots 53, 57, 61, 65 and 69 in TIAMAT dm grids contain anomalies [single
    highly oversense voxels at (0, 0, 0) and extensive slabs with zero density]
    which significantly affect their statistics. The density spikes have been
    remedied by resetting the density of the offending voxel to that of the
    average over its (6) neigbours. Anomalous zero density regions are left
    untreated.
    
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

static long long find_LR_empty_count(double *grid, int LR_dim_val)
{
    long long count_val = 0;
    
    for (int i = 0; i < LR_dim_val; i++)
        for (int j = 0; j < LR_dim_val; j++)
            for (int k = 0; k < LR_dim_val; k++)
            {
                if (*((double *)grid + HII_R_INDEX(i, j, k)) == -1.0) count_val++;
                
            }
    
    return (long long)(count_val);
}

static void write_list_of_HR_empties(double *grid_HR, int HR_dim_val, int snap, long long empty_count, char *out_dir)
{
    FILE *f1;
    char file1[128];
    
    sprintf(file1, "%s/TIAMAT_dm_density_snap%d_%d_EMPTIES.dat", out_dir, snap, empty_count);
    f1 = fopen(file1, "wt");
    
    for (int i = 0; i < HR_dim_val; i++)
        for (int j = 0; j < HR_dim_val; j++)
            for (int k = 0; k < HR_dim_val; k++)
            {
                if (*(grid_HR + HR_INDEX(i, j, k, HR_dim_val)) == 0.0) fprintf(f1, "%d\t%d\t%d\n", i, j, k);
                
            }
    
    fclose(f1); 
}

static void write_list_of_LR_empties(double *grid, int LR_dim_val, int snap, long long empty_count, char *out_dir)
{
    FILE *f1;
    char file1[128];
    
    sprintf(file1, "%s/TIAMAT_LR_overdensity_snap%d_%d_EMPTIES.dat", out_dir, snap, empty_count);
    f1 = fopen(file1, "wt");
    
    for (int i = 0; i < LR_dim_val; i++)
        for (int j = 0; j < LR_dim_val; j++)
            for (int k = 0; k < LR_dim_val; k++)
            {
                if (*(grid + HII_R_INDEX(i, j, k)) == -1.0) fprintf(f1, "%d\t%d\t%d\n", i, j, k);
                
            }
    
    fclose(f1); 
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
    
    SID_log("Reading grid for snapshot %d", SID_LOG_OPEN, snapshot);
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
        
        
        if (snapshot==53 || snapshot==57 || snapshot==61 || snapshot==65 || snapshot==69)
        {
            FILE     *f1_pmg;
            char      file1_pmg[128];
            long long empty_count_HR;
            long long empty_count_LR;
            int       i, j, k;
            
            empty_count_HR = find_HR_empty_count(grid_HR, HR_dim);
            SID_log("empty_count_HR = %d", SID_LOG_COMMENT, empty_count_HR)
            
            if(empty_count_HR)
            {
                write_list_of_HR_empties(grid_HR, HR_dim, snapshot, empty_count_HR, run_globals->params.OutputDir);
                
                if (snapshot==53)
                {
                    i = 580;
                    j = 580;
                    k = 260;
                }
                
                if (snapshot==57)
                {
                    i = 550;
                    j = 130;
                    k = 520;
                }
                
                if (snapshot==61)
                {
                    i = 710;
                    j = 325;
                    k = 770;
                }
                
                if (snapshot==65)
                {
                    i = 1020;
                    j = 620;
                    k = 620;
                }
                
                if (snapshot==69)
                {
                    i = 515;
                    j = 770;
                    k = 260;
                }
                
                write_HR_i_slice(grid_HR, HR_dim, snapshot, i, run_globals->params.OutputDir);
                write_HR_j_slice(grid_HR, HR_dim, snapshot, j, run_globals->params.OutputDir);
                write_HR_k_slice(grid_HR, HR_dim, snapshot, k, run_globals->params.OutputDir);
            }
             
        }
        
        
        
//        FILE *f1_pmg;
//        char file1_pmg[128];
//        long long empty_count_HR;
//        long long empty_count_LR;
        
//        if (snapshot==53)
//        {
//            empty_count_HR = find_HR_empty_count(grid_HR, HR_dim);
//            SID_log("empty_count_HR = %d", SID_LOG_COMMENT, empty_count_HR);
//            
//            
//            if(empty_count_HR)
//            {
//                write_list_of_HR_empties(grid_HR, HR_dim, snapshot, empty_count_HR, run_globals->params.OutputDir);
//                
//                int i = 580;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int j = 0; j < HR_dim; j++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int j = 580;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int k = 260;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int j = 0; j < HR_dim; j++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//            }
//            
//        }
        
//        if (snapshot==57)
//        {
//            empty_count_HR = find_HR_empty_count(grid_HR, HR_dim);
//            SID_log("empty_count_HR = %d", SID_LOG_COMMENT, empty_count_HR);
//            
//            if(empty_count_HR)
//            {
//                write_list_of_HR_empties(grid_HR, HR_dim, snapshot, empty_count_HR, run_globals->params.OutputDir);
//                
//                int i = 550;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int j = 0; j < HR_dim; j++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int j = 130;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int k = 520;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int j = 0; j < HR_dim; j++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//            }
//            
//        }
        
//        if (snapshot==61)
//        {
//            empty_count_HR = find_HR_empty_count(grid_HR, HR_dim);
//            SID_log("empty_count_HR = %d", SID_LOG_COMMENT, empty_count_HR);
//            
//            if(empty_count_HR)
//            {
//                write_list_of_HR_empties(grid_HR, HR_dim, snapshot, empty_count_HR, run_globals->params.OutputDir);
//                
//                int i = 710;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int j = 0; j < HR_dim; j++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int j = 325;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int k = 770;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int j = 0; j < HR_dim; j++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//            }
//            
//        }
        
//        if (snapshot==65)
//        {
//            empty_count_HR = find_HR_empty_count(grid_HR, HR_dim);
//            SID_log("empty_count_HR = %d", SID_LOG_COMMENT, empty_count_HR);
//            
//            if(empty_count_HR)
//            {
//                write_list_of_HR_empties(grid_HR, HR_dim, snapshot, empty_count_HR, run_globals->params.OutputDir);
//                
//                int i = 1020;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int j = 0; j < HR_dim; j++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int j = 620;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int k = 620;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int j = 0; j < HR_dim; j++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//            }
//            
//        }
        
//        if (snapshot==69)
//        {
//            empty_count_HR = find_HR_empty_count(grid_HR, HR_dim);
//            SID_log("empty_count_HR = %d", SID_LOG_COMMENT, empty_count_HR);
//            
//            if(empty_count_HR)
//            {
//                write_list_of_HR_empties(grid_HR, HR_dim, snapshot, empty_count_HR, run_globals->params.OutputDir);
//                
//                int i = 515;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int j = 0; j < HR_dim; j++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int j = 770;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int k = 0; k < HR_dim; k++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//                
//                int k = 260;
//                sprintf(file1_pmg, "%s/TIAMAT_dm_density_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int j = 0; j < HR_dim; j++)
//                {
//                    for (int i = 0; i < HR_dim; i++)
//                    {
//                        fprintf(f1_pmg, "%g\t", *(grid_HR + HR_INDEX(i, j, k, HR_dim)));
//                    }
//                    fprintf(f1_pmg, "\n");
//                }
//                fclose(f1_pmg);
//            }
//            
//        }
        
        
        
        // QUICK FIX FOR TIAMAT !!!
        //
        // From previous analysis, we know that the offending maximum spike voxel for these snapshots is (0, 0, 0)
        // Now reset its value to that of the average over its neigbours
        //
        // If this is a problem snapshot then do the averaging
        if (snapshot==53 || snapshot==57 || snapshot==61 || snapshot==65 || snapshot==69)
        {
            SID_log("Revaluating problem dm density voxel:", SID_LOG_OPEN);
            
            *(grid_HR + HR_INDEX(0, 0, 0, HR_dim)) = ( *(grid_HR + HR_INDEX(1, 0, 0, HR_dim)) +
                                                       *(grid_HR + HR_INDEX(0, 1, 0, HR_dim)) +
                                                       *(grid_HR + HR_INDEX(0, 0, 1, HR_dim)) +
                                                       *(grid_HR + HR_INDEX(HR_dim - 1, 0, 0, HR_dim)) +
                                                       *(grid_HR + HR_INDEX(0, HR_dim - 1, 0, HR_dim)) +
                                                       *(grid_HR + HR_INDEX(0, 0, HR_dim - 1, HR_dim))) / 6.0;
            SID_log("...done", SID_LOG_CLOSE);
        }
        
        
        
        
        
        
        
        
        
        
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
        
        
        if (snapshot==53)
        {
            empty_count_LR = find_LR_empty_count(grid, HII_dim);
            SID_log("empty_count_LR = %d", SID_LOG_COMMENT, empty_count_LR);
            
            if(empty_count_LR)
            {
//                sprintf(file1_pmg, "%s/TIAMAT_LR_overdensity_snap%d_%d_EMPTIES.dat", run_globals->params.OutputDir, snapshot, empty_count_LR);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int i = 0; i < HII_dim; i++)
//                    for (int j = 0; j < HII_dim; j++)
//                        for (int k = 0; k < HII_dim; k++)
//                        {
//                            if (*(grid + HII_R_INDEX(i, j, k)) == -1.0) fprintf(f1_pmg, "%d\t%d\t%d\n", i, j, k);
//                            
//                        }
//                fclose(f1_pmg);
                
                write_list_of_LR_empties(grid, HII_dim, snapshot, empty_count_LR, run_globals->params.OutputDir);
                
                int i = 145;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int j = 0; j < HII_dim; j++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int j = 145;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int k = 65;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int j = 0; j < HII_dim; j++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
            }
            
        }
        
        if (snapshot==57)
        {
            empty_count_LR = find_LR_empty_count(grid, HII_dim);
            SID_log("empty_count_LR = %d", SID_LOG_COMMENT, empty_count_LR);
            
            if(empty_count_LR)
            {
//                sprintf(file1_pmg, "%s/TIAMAT_LR_overdensity_snap%d_%d_EMPTIES.dat", run_globals->params.OutputDir, snapshot, empty_count_LR);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int i = 0; i < HII_dim; i++)
//                    for (int j = 0; j < HII_dim; j++)
//                        for (int k = 0; k < HII_dim; k++)
//                        {
//                            if (*(grid + HII_R_INDEX(i, j, k)) == -1.0) fprintf(f1_pmg, "%d\t%d\t%d\n", i, j, k);
//                            
//                        }
//                fclose(f1_pmg);
                
                write_list_of_LR_empties(grid, HII_dim, snapshot, empty_count_LR, run_globals->params.OutputDir);
                
                int i = 139;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int j = 0; j < HII_dim; j++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int j = 34;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int k = 132;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int j = 0; j < HII_dim; j++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
            }
        }
        
        if (snapshot==61)
        {
            empty_count_LR = find_LR_empty_count(grid, HII_dim);
            SID_log("empty_count_LR = %d", SID_LOG_COMMENT, empty_count_LR);
            
            if(empty_count_LR)
            {
//                sprintf(file1_pmg, "%s/TIAMAT_LR_overdensity_snap%d_%d_EMPTIES.dat", run_globals->params.OutputDir, snapshot, empty_count_LR);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int i = 0; i < HII_dim; i++)
//                    for (int j = 0; j < HII_dim; j++)
//                        for (int k = 0; k < HII_dim; k++)
//                        {
//                            if (*(grid + HII_R_INDEX(i, j, k)) == -1.0) fprintf(f1_pmg, "%d\t%d\t%d\n", i, j, k);
//                            
//                        }
//                fclose(f1_pmg);
                
                write_list_of_LR_empties(grid, HII_dim, snapshot, empty_count_LR, run_globals->params.OutputDir);
                
                int i = 179;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int j = 0; j < HII_dim; j++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int j = 83;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int k = 194;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int j = 0; j < HII_dim; j++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
            }
            
        }
        
        if (snapshot==65)
        {
            empty_count_LR = find_LR_empty_count(grid, HII_dim);
            SID_log("empty_count_LR = %d", SID_LOG_COMMENT, empty_count_LR);
            
            if(empty_count_LR)
            {
//                sprintf(file1_pmg, "%s/TIAMAT_LR_overdensity_snap%d_%d_EMPTIES.dat", run_globals->params.OutputDir, snapshot, empty_count_LR);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int i = 0; i < HII_dim; i++)
//                    for (int j = 0; j < HII_dim; j++)
//                        for (int k = 0; k < HII_dim; k++)
//                        {
//                            if (*(grid + HII_R_INDEX(i, j, k)) == -1.0) fprintf(f1_pmg, "%d\t%d\t%d\n", i, j, k);
//                            
//                        }
//                fclose(f1_pmg);
                
                write_list_of_LR_empties(grid, HII_dim, snapshot, empty_count_LR, run_globals->params.OutputDir);
                
                int i = 255;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int j = 0; j < HII_dim; j++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int j = 157;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int k = 157;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int j = 0; j < HII_dim; j++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
            }
            
        }
        
        if (snapshot==69)
        {
            empty_count_LR = find_LR_empty_count(grid, HII_dim);
            SID_log("empty_count_LR = %d", SID_LOG_COMMENT, empty_count_LR);
            
            if(empty_count_LR)
            {
//                sprintf(file1_pmg, "%s/TIAMAT_LR_overdensity_snap%d_%d_EMPTIES.dat", run_globals->params.OutputDir, snapshot, empty_count_LR);
//                f1_pmg = fopen(file1_pmg, "wt");
//                for (int i = 0; i < HII_dim; i++)
//                    for (int j = 0; j < HII_dim; j++)
//                        for (int k = 0; k < HII_dim; k++)
//                        {
//                            if (*(grid + HII_R_INDEX(i, j, k)) == -1.0) fprintf(f1_pmg, "%d\t%d\t%d\n", i, j, k);
//                            
//                        }
//                fclose(f1_pmg);
                
                write_list_of_LR_empties(grid, HII_dim, snapshot, empty_count_LR, run_globals->params.OutputDir);
                
                int i = 130;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_i%d_snap%d.dat", run_globals->params.OutputDir, i, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int j = 0; j < HII_dim; j++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int j = 194;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_j%d_snap%d.dat", run_globals->params.OutputDir, j, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int k = 0; k < HII_dim; k++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
                
                int k = 66;
                sprintf(file1_pmg, "%s/TIAMAT_LR_density_contrast_grid_slice_k%d_snap%d.dat", run_globals->params.OutputDir, k, snapshot);
                f1_pmg = fopen(file1_pmg, "wt");
                for (int j = 0; j < HII_dim; j++)
                {
                    for (int i = 0; i < HII_dim; i++)
                    {
                        fprintf(f1_pmg, "%g\t", *(grid + HII_R_INDEX(i, j, k)) + 1.0);
                    }
                    fprintf(f1_pmg, "\n");
                }
                fclose(f1_pmg);
            }
            
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
