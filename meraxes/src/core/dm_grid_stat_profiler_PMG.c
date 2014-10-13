#include "meraxes.h"
#include <fftw3.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>


void call_dm_grid_stats_PMG(run_globals_t *run_globals, int snapshot, int nout_gals)
{
    
    int total_n_out_gals = 0;
    
    tocf_grids_t *grids = &(run_globals->tocf_grids);
    
    SID_log("Getting ready to call dm_grid_stats_PMG...", SID_LOG_OPEN);
    
    // Check to see if there are actually any galaxies at this snapshot
    SID_Allreduce(&nout_gals, &total_n_out_gals, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
    if (total_n_out_gals == 0)
    {
        SID_log("No galaxies in the simulation - skipping...", SID_LOG_CLOSE);
        return;
    }
    
    SID_log("...done", SID_LOG_CLOSE);
    
    if (SID.My_rank == 0)
    {
        SID_log("Calling dm_grid_stats_PMG...", SID_LOG_OPEN);
        
        // Read in the dark matter density grid
        dm_grid_stats_PMG(run_globals, snapshot, 0);   // i_grid must be 0 !
        
        SID_log("...done", SID_LOG_CLOSE);
    }
    
}


static inline void read_identifier_PMG(FILE *fin, bool skip_flag)
{
    char identifier[32];
    
    fread(identifier, sizeof(identifier), 1, fin);
    
    if (skip_flag)
        SID_log_error("Skipping grid: %s...", identifier);
    else
        SID_log("Reading grid: %s...", SID_LOG_COMMENT, identifier);
}


int dm_grid_stats_PMG(
    run_globals_t *run_globals,
    int            snapshot,
    int            i_grid)
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
    float  *grid_float;
    
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
    
    SID_log("Reading grid for snapshot %d ------ BUG FIXED TEST STAT WRITER", SID_LOG_OPEN, snapshot);
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
        read_identifier_PMG(fin, true);
        fseek(fin, sizeof(float) * n_elem, SEEK_CUR);
    }
    read_identifier_PMG(fin, false);
    
    // Malloc the grids
    grid       = SID_calloc(sizeof(double) * HII_dim*HII_dim*HII_dim);
    grid_float = SID_calloc(sizeof(float) * HII_dim*HII_dim*HII_dim);
    
    // Initialise the grid
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
            {
                *(grid + HII_R_INDEX(i, j, k)) = 0.0;
                *(grid_float + HII_R_INDEX(i, j, k)) = 0.0;
            }
    
    
    if (i_grid == 0)  // Density grid
    {
        // Calculate the volume of a single high resolution cell
        cell_volume = pow(box_size[0] / (double)n_cell[0], 3);
        
        // Read in the grid
        mean = 0.0;
        for (int i = 0; i < n_cell[0]; i++)
            for (int j = 0; j < n_cell[1]; j++)
                for (int k = 0; k < n_cell[2]; k++)
                {
                    fread(&val, sizeof(float), 1, fin);
                    mean += (double)val;
                    *(grid + HII_R_INDEX((int)(i * resample_factor), (int)(j * resample_factor), (int)(k * resample_factor))) += (double)val;
                }
        
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
    
    
    // Copy the grid (double) to the output (float) and free
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
                *(grid_float + HII_R_INDEX(i, j, k)) = (float)*(grid + HII_R_INDEX(i, j, k));
    
    
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
                delta_LRF = *(grid_float + HII_R_INDEX(i, j, k));
                
                delta_LRF_sum += delta_LRF;
                
                if (delta_LRF<delta_LRF_min) delta_LRF_min = delta_LRF;
                    
                if (delta_LRF>delta_LRF_max) delta_LRF_max = delta_LRF;
                
            }
    
    delta_LRF_mean = delta_LRF_sum/n_elem_LR;
    
    for (int i = 0; i < HII_dim; i++)
        for (int j = 0; j < HII_dim; j++)
            for (int k = 0; k < HII_dim; k++)
            {
                delta_LRF = *(grid_float + HII_R_INDEX(i, j, k));
                delta_LRF_var += (delta_LRF - delta_LRF_mean)*(delta_LRF - delta_LRF_mean);
            }
    
    delta_LRF_var /= ((float)(n_elem_LR - 1.0));
    
    SID_log("delta_LRF_sum  = %g", SID_LOG_COMMENT, delta_LRF_sum);
    SID_log("delta_LRF_min  = %g", SID_LOG_COMMENT, delta_LRF_min);
    SID_log("delta_LRF_max  = %g", SID_LOG_COMMENT, delta_LRF_max);
    SID_log("delta_LRF_mean = %g", SID_LOG_COMMENT, delta_LRF_mean);
    SID_log("delta_LRF_var  = %g", SID_LOG_COMMENT, delta_LRF_var);
    
    // Write overdensity stats to file -- temp
    FILE *f1_pmg;
    char file1_pmg[128];
    sprintf(file1_pmg, "output_TEST/delta_grid_stats_snap%d_TEST.dat", snapshot);
    f1_pmg = fopen(file1_pmg, "wt");
    fprintf(f1_pmg, "min\t%g\n", delta_LRF_min);
    fprintf(f1_pmg, "max\t%g\n", delta_LRF_max);
    fprintf(f1_pmg, "mean\t%g\n", delta_LRF_mean);
    fprintf(f1_pmg, "var\t%g\n", delta_LRF_var);
    fclose(f1_pmg);
    
    // PMG SECTION END
    
    
    SID_free(SID_FARG grid);
    SID_free(SID_FARG grid_float);
    
    SID_log("...done", SID_LOG_CLOSE);
    
    // Close the file
    fclose(fin);
    
    return 0;
}




