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
    
    double *grid_HR;   // PMG: for HR stats
    unsigned long long HR_INDEX;   // PMG: for HR stats
    unsigned long long HR_INDEX_neighbour1;   // PMG: for anomaly fix
    unsigned long long HR_INDEX_neighbour2;   // PMG: for anomaly fix
    unsigned long long HR_INDEX_neighbour3;   // PMG: for anomaly fix
    unsigned long long HR_INDEX_neighbour4;   // PMG: for anomaly fix
    unsigned long long HR_INDEX_neighbour5;   // PMG: for anomaly fix
    unsigned long long HR_INDEX_neighbour6;   // PMG: for anomaly fix
    FILE *f1_pmg;
    char file1_pmg[128];
    
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
    grid    = SID_calloc(sizeof(double) * HII_dim*HII_dim*HII_dim);
    grid_HR = SID_calloc(sizeof(double) * n_elem);   // PMG: for HR stats
    
    // Initialise the grid
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
                *(grid + HII_R_INDEX(i, j, k)) = 0.0;
    
    for (int i = 0; i < n_cell[0]; i++)   // PMG: for HR stats
        for (int j = 0; j < n_cell[1]; j++)
            for (int k = 0; k < n_cell[2]; k++)
            {
                HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));   // PMG: for HR stats
                *(grid_HR + HR_INDEX) = 0.0;   // PMG: for HR stats
            }
    
    
    if (i_grid == 0)  // Density grid
    {
        // Read in the grid and assign grid_HR values
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    fread(&val, sizeof(float), 1, fin);
                    HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));   // PMG: for anomaly fix
                    *(grid_HR + HR_INDEX) = (double)val;   // PMG: for anomaly fix
                }
        
        // PMG: for anomaly fix
        // We now know that the offending maximum spike voxel for these snapshots is (0, 0, 0)
        // Now reset its value to that of the average over its neigbours
        // For these snapshots there are also anomalous empty regions - these are left untreated
        
        // If this is a problem snapshot then do the averaging
        if (snapshot==53 || snapshot==57 || snapshot==61 || snapshot==65 || snapshot==69)
        {
            SID_log("Revaluating problem voxel:", SID_LOG_OPEN);
            
            int i = 0;
            int j = 0;
            int k = 0;
            
            HR_INDEX            = (unsigned long long)(0 + n_cell[2]*(0 + n_cell[1]*0));
            HR_INDEX_neighbour1 = (unsigned long long)(0 + n_cell[2]*(0 + n_cell[1]*1));
            HR_INDEX_neighbour2 = (unsigned long long)(0 + n_cell[2]*(1 + n_cell[1]*0));
            HR_INDEX_neighbour3 = (unsigned long long)(1 + n_cell[2]*(0 + n_cell[1]*0));
            HR_INDEX_neighbour4 = (unsigned long long)(0 + n_cell[2]*(0 + n_cell[1]*(n_cell[0]-1)));
            HR_INDEX_neighbour5 = (unsigned long long)(0 + n_cell[2]*((n_cell[1]-1) + n_cell[1]*0));
            HR_INDEX_neighbour6 = (unsigned long long)((n_cell[2]-1) + n_cell[2]*(0 + n_cell[1]*0));
            
            *(grid_HR + HR_INDEX) = ( *(grid_HR + HR_INDEX_neighbour1) +
                                      *(grid_HR + HR_INDEX_neighbour2) +
                                      *(grid_HR + HR_INDEX_neighbour3) +
                                      *(grid_HR + HR_INDEX_neighbour4) +
                                      *(grid_HR + HR_INDEX_neighbour5) +
                                      *(grid_HR + HR_INDEX_neighbour6)) / 6.0;
            
            SID_log("...done", SID_LOG_CLOSE);
        }
        
        // PMG: for anomaly fix
        // Regrid
        mean = 0.0;
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
                    mean += *(grid_HR + HR_INDEX);
                    *(grid + HII_R_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += *(grid_HR + HR_INDEX);
                }
        
        // Calculate the volume of a single high resolution cell
        cell_volume = pow(box_size[0] / (double)n_cell[0], 3);
        
        // Mean density from high res grid
        mean *= cell_volume / pow(box_size[0], 3);
        
        // PMG HR stats added START
        // At this point grid_HR holds the densities in each HR cell
        // Loop through again and calculate the overdensity
        // i.e. (rho - rho_mean)/rho_mean
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
                    *(grid_HR + HR_INDEX) = (*(grid_HR + HR_INDEX))/mean - 1.0;
                }
        // PMG HR stats added END
        
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
    
    
    // PMG SECTION START
    // -----------------------------------------------------------------
    // -----------------------------------------------------------------
    
    // We find that the dm grids of Tiamat snapshots 53, 57, 61, 65, 69 have single outlier voxels
    // We want to find these and reassign their overdensity value with the average of their neighbours
    
    // HRD basic stats (min, max, mean, variance)
    // --------------------------------------
    
    double delta_HRD;
    double delta_HRD_min = 1.0e10;
    double delta_HRD_max = -1.0e10;
    double delta_HRD_sum = 0.0;
    double delta_HRD_mean;
    double delta_HRD_var = 0.0;
    double n_elem_HR = (double)(n_elem);
    
    for (int i = 0; i < n_cell[0]; i++)
        for (int j = 0; j < n_cell[1]; j++)
            for (int k = 0; k < n_cell[2]; k++)
            {
                HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
                delta_HRD = *(grid_HR + HR_INDEX);
                delta_HRD_sum += delta_HRD;
                
                if (delta_HRD<delta_HRD_min) delta_HRD_min = delta_HRD;
                    
                if (delta_HRD>delta_HRD_max) delta_HRD_max = delta_HRD;
                
            }
    
    delta_HRD_mean = delta_HRD_sum/n_elem_HR;
    
    for (int i = 0; i < n_cell[0]; i++)
        for (int j = 0; j < n_cell[1]; j++)
            for (int k = 0; k < n_cell[2]; k++)
            {
                HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
                delta_HRD = *(grid_HR + HR_INDEX);
                delta_HRD_var += (delta_HRD - delta_HRD_mean)*(delta_HRD - delta_HRD_mean);
            }
    
    delta_HRD_var /= ((double)(n_elem_HR - 1.0));
    
    SID_log("delta_HRD_sum  = %g", SID_LOG_COMMENT, delta_HRD_sum);
    SID_log("delta_HRD_min  = %g", SID_LOG_COMMENT, delta_HRD_min);
    SID_log("delta_HRD_max  = %g", SID_LOG_COMMENT, delta_HRD_max);
    SID_log("delta_HRD_mean = %g", SID_LOG_COMMENT, delta_HRD_mean);
    SID_log("delta_HRD_var  = %g", SID_LOG_COMMENT, delta_HRD_var);
    
    // Write overdensity stats to file
    sprintf(file1_pmg, "%s/delta_HR_grid_stats_snap%d.dat", run_globals->params.OutputDir, snapshot);
    f1_pmg = fopen(file1_pmg, "wt");
    fprintf(f1_pmg, "min\t%g\n", delta_HRD_min);
    fprintf(f1_pmg, "max\t%g\n", delta_HRD_max);
    fprintf(f1_pmg, "mean\t%g\n", delta_HRD_mean);
    fprintf(f1_pmg, "var\t%g\n", delta_HRD_var);
    fclose(f1_pmg);
    
    
    // HRD histogram
    // --------------------------------------
    
    unsigned long bin_count_HRD;
    int    n_bins_HRD = 100;
    double delta_HRD_bin_width = (delta_HRD_max + 1.0)/(double)(n_bins_HRD);
    double delta_HRD_bin_min, delta_HRD_bin_mid, delta_HRD_bin_max;
    
//    SID_log("n_bins_HRD           = %d", SID_LOG_COMMENT, n_bins_HRD);
//    SID_log("delta_HRD_bin_width  = %g", SID_LOG_COMMENT, delta_HRD_bin_width);
//    SID_log("Bins for HR grid:", SID_LOG_OPEN);
    
    sprintf(file1_pmg, "%s/delta_HR_grid_hist_snap%d.dat", run_globals->params.OutputDir, snapshot);
    f1_pmg = fopen(file1_pmg, "wt");
    
    for (int bin = 0; bin <= n_bins_HRD; bin++)
    {
        delta_HRD_bin_min = -1.0 + (double)(bin)*delta_HRD_bin_width;
        delta_HRD_bin_mid = delta_HRD_bin_min + delta_HRD_bin_width/2.0;
        delta_HRD_bin_max = delta_HRD_bin_min + delta_HRD_bin_width;
        
        bin_count_HRD = 0;
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
                    delta_HRD = *(grid_HR + HR_INDEX);
                    
                    if (delta_HRD >= delta_HRD_bin_min && delta_HRD < delta_HRD_bin_max)
                        bin_count_HRD++;
                    
                }
        
        fprintf(f1_pmg, "%f\t%d\n", delta_HRD_bin_mid, bin_count_HRD);
//        SID_log("%3d\t%f\t%f\t%f\t%d", SID_LOG_COMMENT, bin, delta_HRD_bin_min, delta_HRD_bin_mid, delta_HRD_bin_max, bin_count_HRD);
    }
    fclose(f1_pmg);
    
//    SID_log("...done", SID_LOG_CLOSE);
    
    
    // Find bad voxel (maximum) if this is a problem snapshot
    // --------------------------------------
    
//    int i_max = -1;
//    int j_max = -1;
//    int k_max = -1;
//    delta_HRD_max = -1.0e10;
//    
//    if (snapshot==53 || snapshot==57 || snapshot==61 || snapshot==65 || snapshot==69)
//    {
//        SID_log("Finding problem voxel (maximum):", SID_LOG_OPEN);
//        
//        for (int i = 0; i < n_cell[0]; i++)
//            for (int j = 0; j < n_cell[1]; j++)
//                for (int k = 0; k < n_cell[2]; k++)
//                {
//                    HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
//                    delta_HRD = *(grid_HR + HR_INDEX);
//                    
//                    if (delta_HRD>delta_HRD_max)
//                    {
//                        delta_HRD_max = delta_HRD;
//                        i_max = i;
//                        j_max = j;
//                        k_max = k;
//                    }
//                    
//                }
//        
//        SID_log("delta(%d, %d, %d) = %g", SID_LOG_COMMENT, i_max, j_max, k_max, delta_HRD_max);
//        SID_log("...done", SID_LOG_CLOSE);
//        
//        sprintf(file1_pmg, "%s/TIAMAT_dmgrid_anomaly_max_snap%d.dat", run_globals->params.OutputDir, snapshot);
//        f1_pmg = fopen(file1_pmg, "wt");
//        fprintf(f1_pmg, "%d\t%d\t%d\n", i_max, j_max, k_max);
//        fclose(f1_pmg);
//    }
    
    
    // Find bad voxel/s (-1) if this is a problem snapshot
    // --------------------------------------
    
//    if (snapshot==53 || snapshot==57 || snapshot==61 || snapshot==65 || snapshot==69)
//    {
//        SID_log("Finding problem voxel/s (-1):", SID_LOG_OPEN);
//        
//        sprintf(file1_pmg, "%s/TIAMAT_dmgrid_anomaly_min_snap%d.dat", run_globals->params.OutputDir, snapshot);
//        f1_pmg = fopen(file1_pmg, "wt");
//        
//        for (int i = 0; i < n_cell[0]; i++)
//            for (int j = 0; j < n_cell[1]; j++)
//                for (int k = 0; k < n_cell[2]; k++)
//                {
//                    HR_INDEX = (unsigned long long)(k + n_cell[2]*(j + n_cell[1]*i));
//                    delta_HRD = *(grid_HR + HR_INDEX);
//                    
//                    if (delta_HRD==-1.0)
//                    {
//                        SID_log("delta(%d, %d, %d) = -1", SID_LOG_COMMENT, i, j, k);
//                        fprintf(f1_pmg, "%d\t%d\t%d\n", i, j, k);
//                    }
//                    
//                }
//        
//        SID_log("...done", SID_LOG_CLOSE);
//        
//        fclose(f1_pmg);
//    }
    
    
    
    
    // LRF basic stats (min, max, mean, variance)
    // --------------------------------------
    
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
    
    // Write overdensity stats to file
    sprintf(file1_pmg, "%s/delta_LR_grid_stats_snap%d.dat", run_globals->params.OutputDir, snapshot);
    f1_pmg = fopen(file1_pmg, "wt");
    fprintf(f1_pmg, "min\t%g\n", delta_LRF_min);
    fprintf(f1_pmg, "max\t%g\n", delta_LRF_max);
    fprintf(f1_pmg, "mean\t%g\n", delta_LRF_mean);
    fprintf(f1_pmg, "var\t%g\n", delta_LRF_var);
    fclose(f1_pmg);
    
    
    // LRF histogram
    // --------------------------------------
    
    unsigned long bin_count_LRF;
    int   n_bins_LRF = 100;
    float delta_LRF_bin_width = (delta_LRF_max + 1.0)/(float)(n_bins_LRF);
    float delta_LRF_bin_min, delta_LRF_bin_mid, delta_LRF_bin_max;
    
    sprintf(file1_pmg, "%s/delta_LR_grid_hist_snap%d.dat", run_globals->params.OutputDir, snapshot);
    f1_pmg = fopen(file1_pmg, "wt");
    
    for (int bin = 0; bin <= n_bins_LRF; bin++)
    {
        delta_LRF_bin_min = -1.0 + (float)(bin)*delta_LRF_bin_width;
        delta_LRF_bin_mid = delta_LRF_bin_min + delta_LRF_bin_width/2.0;
        delta_LRF_bin_max = delta_LRF_bin_min + delta_LRF_bin_width;
        
        bin_count_LRF = 0;
        for (int i = 0; i < HII_dim; i++)
            for (int j = 0; j < HII_dim; j++)
                for (int k = 0; k < HII_dim; k++)
                {
                    delta_LRF = *(grid_out + HII_R_FFT_INDEX(i, j, k));
                    
                    if (delta_LRF >= delta_LRF_bin_min && delta_LRF < delta_LRF_bin_max)
                        bin_count_LRF++;
                    
                }
        
        fprintf(f1_pmg, "%f\t%d\n", delta_LRF_bin_mid, bin_count_LRF);
    }
    fclose(f1_pmg);
    
    // PMG SECTION END
    // -----------------------------------------------------------------
    // -----------------------------------------------------------------
    
    
    SID_free(SID_FARG grid);
    SID_free(SID_FARG grid_HR);   // PMG: for HR stats
    
    SID_log("...done", SID_LOG_CLOSE);
    
    // Close the file
    fclose(fin);
    
    return 0;
}
#endif
