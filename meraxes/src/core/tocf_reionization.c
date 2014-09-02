#ifdef USE_TOCF

#include "meraxes.h"
#include <fftw3.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <assert.h>

void set_HII_eff_factor(run_globals_t *run_globals)
{
  // Use the params passed to Meraxes via the input file to set the HII
  // ionizing efficiency factor.

  physics_params_t *params = &(run_globals->params.physics);

  tocf_params.HII_eff_factor *= (params->ReionNionPhotPerBary / 4000.0) *
                               (params->ReionEscapeFrac / 0.15) * (1.0 / (1.0 + params->ReionMeanNRec));
  
  tocf_params.HII_eff_factor = 450000;//25000;   // PMG FORCED to make a reasonable ionisation field at snapshot 18
  
}


// NOTE: I am following the indexing conventions off 21cmFAST here.  I have no
// idea what the unsigned long long memory offset stuff is all about
// with the fftwf_complex arrays that have been cast as floats.  Why not still
// just use common indexing in square brackets?
// i.e. ((float *)deltax)[k+(2*(HII_dim/2+1))*(j+HII_dim*i)]

void call_find_HII_bubbles(run_globals_t *run_globals, int snapshot, int nout_gals)
{
  // Thin wrapper round find_HII_bubbles

  int total_n_out_gals = 0;

  tocf_grids_t *grids = &(run_globals->tocf_grids);

  SID_log("Getting ready to call find_HII_bubbles...", SID_LOG_OPEN);

  // Check to see if there are actually any galaxies at this snapshot
  SID_Allreduce(&nout_gals, &total_n_out_gals, 1, SID_INT, SID_SUM, SID.COMM_WORLD);
  if (total_n_out_gals == 0)
  {
    SID_log("No galaxies in the simulation - skipping...", SID_LOG_CLOSE);
    return;
  }

  // Construct the stellar mass grid
  construct_stellar_grids(run_globals);
    
  SID_log("...done", SID_LOG_CLOSE);

  if (SID.My_rank == 0)
  {
    // Read in the dark matter density grid
    read_dm_grid(run_globals, snapshot, 0, (float*)(grids->deltax));

    // Make copies of the stellar and deltax grids before sending them to
    // 21cmfast.  This is because the floating precision fft--ifft introduces
    // rounding error that we don't want to store...
    memcpy(grids->stars_copy , grids->stars , sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    memcpy(grids->sfr_copy   , grids->sfr   , sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    memcpy(grids->deltax_copy, grids->deltax, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);

    SID_log("Calling find_HII_bubbles...", SID_LOG_OPEN);
    // TODO: Fix if snapshot==0
    grids->global_xH = find_HII_bubbles(run_globals->ZZ[snapshot], run_globals->ZZ[snapshot - 1],
        grids->xH,
        grids->stars,
        grids->stars_filtered,
        grids->deltax,
        grids->deltax_filtered,
        grids->sfr,                 // grids->sfr or NULL
        grids->sfr_filtered,        // grids->sfr_filtered or NULL
        grids->z_at_ionization,
        grids->J_21_at_ionization,
        grids->J_21,
        grids->mfp,
        grids->N_rec,
        grids->N_rec_filtered
        );
    
    SID_log("grids->global_xH = %g", SID_LOG_COMMENT, grids->global_xH);
    
    // copy the original (non fourier transformed) grids back into place
    memcpy(grids->stars , grids->stars_copy , sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    memcpy(grids->sfr   , grids->sfr_copy   , sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    memcpy(grids->deltax, grids->deltax_copy, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
  }

  // send the global_xH value to all cores
  SID_Bcast(&(grids->global_xH), sizeof(float), 0, SID.COMM_WORLD);

  SID_log("...done", SID_LOG_CLOSE);
}


void malloc_reionization_grids(run_globals_t *run_globals)
{
  tocf_grids_t *grids = &(run_globals->tocf_grids);

  grids->xH                 = NULL;
  grids->stars              = NULL;
  grids->stars_filtered     = NULL;
  grids->deltax             = NULL;
  grids->deltax_filtered    = NULL;
  grids->sfr                = NULL;
  grids->sfr_filtered       = NULL;
  grids->z_at_ionization    = NULL;
  grids->J_21_at_ionization = NULL;
  grids->J_21               = NULL;
  grids->mfp                = NULL;
  grids->N_rec              = NULL;
  grids->N_rec_filtered     = NULL;

  grids->global_xH = 1.0;

  if (run_globals->params.TOCF_Flag)
  {
    if (SID.My_rank == 0)
    {
      // Only rank 0 will call 21cmFAST and so it is the only rank that needs to carry around these grids
      grids->xH              = (float*)fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
      grids->stars_filtered  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->stars_copy      = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->deltax          = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->deltax_filtered = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->deltax_copy     = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->sfr_filtered    = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->sfr_copy        = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    }
    grids->stars           = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    grids->sfr             = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);

    if (tocf_params.uvb_feedback)
    {
      if (SID.My_rank == 0)
      {
        // Only rank 0 will call 21cmFAST and so it is the only rank that needs to carry around these grids
        grids->z_at_ionization    = (float*)fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
        grids->J_21_at_ionization = (float*)fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
        grids->J_21               = (float*)fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
      }
      grids->Mvir_crit          = (float*)fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
    }

    if (tocf_params.compute_mfp && (SID.My_rank == 0))
    {
      // Only rank 0 will call 21cmFAST and so it is the only rank that needs to carry around these grids
      grids->mfp            = (float*)fftwf_malloc(sizeof(float) * HII_TOT_NUM_PIXELS);
      grids->N_rec          = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      grids->N_rec_filtered = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    }

    SID_log("Initialising grids...", SID_LOG_COMMENT);

    if (SID.My_rank == 0)
    {
      for (int ii = 0; ii < HII_TOT_NUM_PIXELS; ii++)
        grids->xH[ii] = 1.0;
    }

    if (tocf_params.uvb_feedback)
    {
      if (SID.My_rank == 0)
      {
        memset(grids->J_21_at_ionization, 0., sizeof(float) * HII_TOT_NUM_PIXELS);
        memset(grids->J_21              , 0., sizeof(float) * HII_TOT_NUM_PIXELS);
        for (int ii = 0; ii < HII_TOT_NUM_PIXELS; ii++)
          grids->z_at_ionization[ii] = -1;
      }

      memset(grids->Mvir_crit         , 0 , sizeof(float) * HII_TOT_NUM_PIXELS);
    }

    if (SID.My_rank == 0)
    {
      memset(grids->stars_filtered , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      memset(grids->stars_copy     , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      memset(grids->deltax         , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      memset(grids->deltax_filtered, 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      memset(grids->deltax_copy    , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      memset(grids->sfr_filtered   , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      memset(grids->sfr_copy       , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);

      if (tocf_params.compute_mfp)
      {
        memset(grids->mfp           , 0., sizeof(float) * HII_TOT_NUM_PIXELS);
        memset(grids->N_rec         , 0 , sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
        memset(grids->N_rec_filtered, 0 , sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
      }
    }

    memset(grids->stars          , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);
    memset(grids->sfr            , 0, sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS);


    SID_log(" ...done", SID_LOG_CLOSE);
  }
}


void free_reionization_grids(run_globals_t *run_globals)
{
  SID_log("Freeing reionization grids...", SID_LOG_OPEN);

  tocf_grids_t *grids = &(run_globals->tocf_grids);

  if (SID.My_rank == 0)
  {
    if (tocf_params.compute_mfp)
    {
      fftwf_free(grids->N_rec_filtered);
      fftwf_free(grids->N_rec);
      fftwf_free(grids->mfp);
    }
    if (tocf_params.uvb_feedback)
    {
      fftwf_free(grids->J_21);
      fftwf_free(grids->J_21_at_ionization);
      fftwf_free(grids->z_at_ionization);
    }
    fftwf_free(grids->sfr_filtered);
    fftwf_free(grids->sfr_copy);
    fftwf_free(grids->deltax_filtered);
    fftwf_free(grids->deltax_copy);
    fftwf_free(grids->deltax);
    fftwf_free(grids->stars_filtered);
    fftwf_free(grids->stars_copy);
    fftwf_free(grids->xH);
  }

  if (tocf_params.uvb_feedback)
    fftwf_free(grids->Mvir_crit);

  fftwf_free(grids->stars);
  fftwf_free(grids->sfr);

  SID_log(" ...done", SID_LOG_CLOSE);
}


int find_cell(float pos, double box_size)
{
  int HII_dim = tocf_params.HII_dim;

  int cell = (int)floor(pos / box_size * (double)HII_dim);

  cell = (cell < 0) ? 0 : cell;
  cell = (cell >= HII_dim) ? HII_dim-1 : cell;

  return cell;
}


void construct_stellar_grids(run_globals_t *run_globals)
{
  galaxy_t *gal;
  int i, j, k;
  double box_size      = (double)(run_globals->params.BoxSize);
  double Hubble_h      = run_globals->params.Hubble_h;
  float *stellar_grid  = (float*)(run_globals->tocf_grids.stars);
  float *sfr_grid      = (float*)(run_globals->tocf_grids.sfr);
  int HII_dim          = tocf_params.HII_dim;
  run_units_t *units   = &(run_globals->units);

  SID_log("Constructing stellar mass and sfr grids...", SID_LOG_OPEN);

  // init the grid
  for (int ii = 0; ii < HII_TOT_FFT_NUM_PIXELS; ii++)
  {
    *(stellar_grid + ii) = 0.0;
    *(sfr_grid + ii)     = 0.0;
  }

  // Loop through each valid galaxy and add its stellar mass to the appropriate cell
  gal = run_globals->FirstGal;
  while (gal != NULL)
  {
    // TODO: Note that I am including ghosts here.  We will need to check the
    // validity of this.  By definition, if they are ghosts then their host
    // halo hasn't been identified at this time step and hence they haven't
    // been evolved.  Their properties (Sfr, StellarMass, etc.) will all have
    // been set when they were last identified.
    if (gal->Type < 3)
    {
      // TODO: for type 2 galaxies these positions will be set from the last
      // time they were identified.  If their host halo has moved significantly
      // since then, these positions won't reflect that and the satellites will
      // be spatially disconnected from their hosts.  We will need to fix this
      // at some point.

      // TODO: Get Greg to fix these positions to obey PBC!!
      if (gal->Pos[0] >= box_size)
        gal->Pos[0] -= box_size;
      else if (gal->Pos[0] < 0.0)
        gal->Pos[0] += box_size;
      i = find_cell(gal->Pos[0], box_size);

      if (gal->Pos[1] >= box_size)
        gal->Pos[1] -= box_size;
      else if (gal->Pos[1] < 0.0)
        gal->Pos[1] += box_size;
      j = find_cell(gal->Pos[1], box_size);

      if (gal->Pos[2] >= box_size)
        gal->Pos[2] -= box_size;
      else if (gal->Pos[2] < 0.0)
        gal->Pos[2] += box_size;
      k = find_cell(gal->Pos[2], box_size);

      assert((i >= 0) && (i < HII_dim));
      assert((j >= 0) && (j < HII_dim));
      assert((k >= 0) && (k < HII_dim));

      *(stellar_grid + HII_R_FFT_INDEX(i, j, k)) += gal->StellarMass;
      *(sfr_grid + HII_R_FFT_INDEX(i, j, k))     += gal->Sfr;
    }
    gal = gal->Next;
  }

  // Collect all grid cell values onto rank 0 which will actually call 21cmFAST
  SID_Allreduce(SID_IN_PLACE, stellar_grid, HII_TOT_FFT_NUM_PIXELS, SID_FLOAT, SID_SUM, SID.COMM_WORLD);
  SID_Allreduce(SID_IN_PLACE, sfr_grid, HII_TOT_FFT_NUM_PIXELS, SID_FLOAT, SID_SUM, SID.COMM_WORLD);

  if (SID.My_rank == 0)
  {
    // Do one final pass to put the grid in the correct (real) units (Msol or Msol/s)
    for (int i = 0; i < HII_dim; i++)
      for (int j = 0; j < HII_dim; j++)
        for (int k = 0; k < HII_dim; k++)
        {
          if (*(stellar_grid + HII_R_FFT_INDEX(i, j, k)) > 0)
            *(stellar_grid + HII_R_FFT_INDEX(i, j, k)) *= (1.e10 / Hubble_h);
          if (*(sfr_grid + HII_R_FFT_INDEX(i, j, k)) > 0)
            *(sfr_grid + HII_R_FFT_INDEX(i, j, k)) *= (units->UnitMass_in_g / units->UnitTime_in_s / SOLAR_MASS);

          // Check for under/overflow
          if (*(stellar_grid + HII_R_FFT_INDEX(i, j, k)) < 0)
            *(stellar_grid + HII_R_FFT_INDEX(i, j, k)) = 0;
          if (*(sfr_grid + HII_R_FFT_INDEX(i, j, k)) < 0)
            *(sfr_grid + HII_R_FFT_INDEX(i, j, k)) = 0;
        }
  }

  SID_log("done", SID_LOG_CLOSE);
}


void save_tocf_grids(run_globals_t *run_globals, hid_t parent_group_id, int snapshot)
{    
    if (SID.My_rank == 0)
    {
        tocf_grids_t *grids = &(run_globals->tocf_grids);
        hsize_t dims        = HII_TOT_NUM_PIXELS;
        int   HII_dim       = tocf_params.HII_dim;
        float *grid;
        float *ps;
        float *delta_T;
        int   ps_nbins;
        float average_deltaT;
        hid_t group_id;
        
        int   size_dist_nbins;
        float *size_dist = NULL;
        
        
        // Save tocf grids
        // ----------------------------------------------------------------------------------------------------
        
        SID_log("Saving tocf grids...", SID_LOG_OPEN);
        
        group_id = H5Gcreate(parent_group_id, "Grids", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // float grids
        H5LTmake_dataset_float(group_id, "xH", 1, &dims, grids->xH);
        
        if (tocf_params.uvb_feedback)
        {
            H5LTmake_dataset_float(group_id, "J_21", 1, &dims, grids->J_21);
            H5LTmake_dataset_float(group_id, "J_21_at_ionization", 1, &dims, grids->J_21_at_ionization);
            H5LTmake_dataset_float(group_id, "z_at_ionization", 1, &dims, grids->z_at_ionization);
            H5LTmake_dataset_float(group_id, "Mvir_crit", 1, &dims, grids->Mvir_crit);
        }
        
        if (tocf_params.compute_mfp)
            H5LTmake_dataset_float(group_id, "MFP", 1, &dims, grids->mfp);
        
        H5LTset_attribute_float(group_id, "xH", "global_xH", &(grids->global_xH), 1);
        
        // fftw padded grids
        grid = (float*)SID_calloc(HII_TOT_NUM_PIXELS * sizeof(float));
        
        // for (int ii = 0; ii < HII_dim; ii++)
        //   for (int jj = 0; jj < HII_dim; jj++)
        //     for (int kk = 0; kk < HII_dim; kk++)
        //       grid[HII_R_INDEX(ii, jj, kk)] = *((float*)(grids->stars) + HII_R_FFT_INDEX(ii, jj, kk));
        // H5LTmake_dataset_float(group_id, "StellarMass", 1, &dims, grid);
        
        // memset((void*)grid, 0, sizeof(float) * HII_TOT_NUM_PIXELS);
        // for (int ii = 0; ii < HII_dim; ii++)
        //   for (int jj = 0; jj < HII_dim; jj++)
        //     for (int kk = 0; kk < HII_dim; kk++)
        //       grid[HII_R_INDEX(ii, jj, kk)] = *((float*)(grids->sfr) + HII_R_FFT_INDEX(ii, jj, kk)) * SEC_PER_YEAR;
        // H5LTmake_dataset_float(group_id, "Sfr", 1, &dims, grid);
        
        memset((void*)grid, 0, sizeof(float) * HII_TOT_NUM_PIXELS);
        for (int ii = 0; ii < HII_dim; ii++)
            for (int jj = 0; jj < HII_dim; jj++)
                for (int kk = 0; kk < HII_dim; kk++)
                    grid[HII_R_INDEX(ii, jj, kk)] = *((float*)(grids->deltax) + HII_R_FFT_INDEX(ii, jj, kk));
        H5LTmake_dataset_float(group_id, "deltax", 1, &dims, grid);
        
        if (tocf_params.compute_mfp)
        {
            memset((void*)grid, 0, sizeof(float) * HII_TOT_NUM_PIXELS);
            for (int ii = 0; ii < HII_dim; ii++)
                for (int jj = 0; jj < HII_dim; jj++)
                    for (int kk = 0; kk < HII_dim; kk++)
                        grid[HII_R_INDEX(ii, jj, kk)] = *((float*)(grids->N_rec) + HII_R_FFT_INDEX(ii, jj, kk));
            H5LTmake_dataset_float(group_id, "N_rec", 1, &dims, grid);
        }
        
        
        // Run delta_T_ps
        // ----------------------------------------------------------------------------------------------------
        
        SID_log("Calculating delta_T box and power spectrum...", SID_LOG_OPEN);
        
        memset((void*)grid, 0, sizeof(float) * HII_TOT_NUM_PIXELS);
        delta_T_ps(
            run_globals->ZZ[snapshot],
            tocf_params.numcores,
            grids->xH,
            (float*)(grids->deltax),
            NULL,
            NULL,
            &average_deltaT,
            grid,
            &ps,
            &ps_nbins);
        
        H5LTmake_dataset_float(group_id, "delta_T", 1, &dims, grid);
        
        dims = ps_nbins * 3;
        H5LTmake_dataset_float(parent_group_id , "PowerSpectrum", 1               , &dims          , ps);
        H5LTset_attribute_int(parent_group_id  , "PowerSpectrum", "nbins"         , &ps_nbins      , 1);
        H5LTset_attribute_float(parent_group_id, "PowerSpectrum", "average_deltaT", &average_deltaT, 1);
        
        // START PMG ADDED
        // ====================================================
        // Remove this when python plotting is sorted?
        
        // Create ps file
        
        FILE *f1;
        char filename[128];
        int i, j, k, i_cut, j_cut, k_cut;
        
        // Create delta_T slice file
        
        i_cut = (int)(HII_dim/2);
        j_cut = (int)(HII_dim/2);
        k_cut = (int)(HII_dim/2);
        
        sprintf(filename, "output/delta_T_i%d_%d.dat", i_cut, snapshot);
        f1 = fopen(filename, "wt");
        for(k=0;k<HII_dim;k++)
        {   
            for(j=0;j<HII_dim;j++)
            {
                fprintf(f1, "%g\t", grid[HII_R_INDEX(i_cut,j,k)]);
            }
            fprintf(f1, "\n");
        }
        fclose(f1);
        SID_log("Creating %s", SID_LOG_COMMENT, filename);
        
        sprintf(filename, "output/delta_T_j%d_%d.dat", j_cut, snapshot);
        f1 = fopen(filename, "wt");
        for(k=0;k<HII_dim;k++)
        {   
            for(i=0;i<HII_dim;i++)
            {
                fprintf(f1, "%g\t", grid[HII_R_INDEX(i,j_cut,k)]);
            }
            fprintf(f1, "\n");
        }
        fclose(f1);
        SID_log("Creating %s", SID_LOG_COMMENT, filename);
        
        sprintf(filename, "output/delta_T_k%d_%d.dat", k_cut, snapshot);
        f1 = fopen(filename, "wt");
        for(j=0;j<HII_dim;j++)
        {   
            for(i=0;i<HII_dim;i++)
            {
                fprintf(f1, "%g\t", grid[HII_R_INDEX(i,j,k_cut)]);
            }
            fprintf(f1, "\n");
        }
        fclose(f1);
        SID_log("Creating %s", SID_LOG_COMMENT, filename);
        
        // Create ps file
        sprintf(filename, "output/ps_%d.dat", snapshot);
        f1 = fopen(filename, "wt");
        for(i=1;i<ps_nbins;i++) fprintf(f1, "%g\t%g\t%g\n", ps[0+3*i], ps[1+3*i], ps[2+3*i]);
        fclose(f1);
        SID_log("Creating %s", SID_LOG_COMMENT, filename);
        
        // END PMG ADDED
        // ====================================================
        
        free(delta_T);
        free(ps);
        
        SID_log("...done", SID_LOG_CLOSE);   // delta_T
        
        
        // Run call_gen_size_dist if neutral fraction is (0, 0.8]
        // ----------------------------------------------------------------------------------------------------
        
        if((run_globals->tocf_grids.global_xH <= 0.95) && (run_globals->tocf_grids.global_xH != 0.0))
        {   
            SID_log("Calling call_gen_size_dist for snapshot %d...", SID_LOG_OPEN, snapshot);
            
            if (tocf_params.region_type_flag) SID_log("Region type: neutral", SID_LOG_COMMENT);
            else SID_log("Region type: ionised", SID_LOG_COMMENT);
            
            call_gen_size_dist(run_globals, snapshot, tocf_params.region_type_flag, &size_dist, &size_dist_nbins);
            
            SID_log("size_dist_nbins = %d", SID_LOG_COMMENT, size_dist_nbins);
            
            // ====================================================
            // Remove this when python plotting is sorted?
            
            // Create size_dist file
            if (tocf_params.region_type_flag) sprintf(filename, "output/size_dist_neutral_%d.dat", snapshot);
            else                              sprintf(filename, "output/size_dist_ionised_%d.dat", snapshot);
            
            f1 = fopen(filename, "wt");
            for(i=0;i<size_dist_nbins;i++) fprintf(f1, "%g\t%g\n", size_dist[0+2*i], size_dist[1+2*i]);
            fclose(f1);
            SID_log("Creating %s", SID_LOG_COMMENT, filename);
            
            // ====================================================
            
            // Write call_gen_size_dist output to hdf5 file
            dims = size_dist_nbins * 2;
            H5LTmake_dataset_float(parent_group_id, "RegionSizeDist", 1,             &dims,                         size_dist);
            H5LTset_attribute_int(parent_group_id , "RegionSizeDist", "nbins",       &size_dist_nbins,              1);
            H5LTset_attribute_int(parent_group_id , "RegionSizeDist", "region_type", &tocf_params.region_type_flag, 1);
            
            free(size_dist);   // allocated within gen_size_dist
            
            SID_log("...done", SID_LOG_CLOSE);   // gen_size_dist
        }
        else
        {
            SID_log("NOT calling call_gen_size_dist for snapshot %d (out of global_xH range)", SID_LOG_COMMENT, snapshot);
        }
        
//        SID_log("Calling call_gen_size_dist for snapshot %d...", SID_LOG_OPEN, snapshot);
//        
//        if (tocf_params.region_type_flag) SID_log("Region type: neutral", SID_LOG_COMMENT);
//        else SID_log("Region type: ionised", SID_LOG_COMMENT);
//        
//        call_gen_size_dist(run_globals, snapshot, tocf_params.region_type_flag, &size_dist, &size_dist_nbins);
//        
//        // If ionisation field is single phase then size_dist won't have been written to
//        // and will still be NULL (so don't try to write to file!)
//        if (size_dist != NULL)
//        {
//            SID_log("size_dist_nbins = %d", SID_LOG_COMMENT, size_dist_nbins);
//        
//            //for (i=0; i<size_dist_nbins; i++) SID_log("%3d\t%1.3f\t%g", SID_LOG_COMMENT, i, size_dist[0+2*i], size_dist[1+2*i]);
//        
//            // Create ps file
//            
//            // START PMG ADDED
//            // ====================================================
//            // Remove this when python plotting is sorted?
//            
//            // Create size_dist file
//            if (tocf_params.region_type_flag) sprintf(filename, "output/size_dist_neutral_%d.dat", snapshot);
//            else                              sprintf(filename, "output/size_dist_ionised_%d.dat", snapshot);
//            
//            f1 = fopen(filename, "wt");
//            for(i=0;i<size_dist_nbins;i++) fprintf(f1, "%g\t%g\n", size_dist[0+2*i], size_dist[1+2*i]);
//            fclose(f1);
//            SID_log("Creating %s", SID_LOG_COMMENT, filename);
//            
//            // END PMG ADDED
//            // ====================================================
//            
//            // Write call_gen_size_dist output to hdf5 file
//            dims = size_dist_nbins * 2;
//            H5LTmake_dataset_float(parent_group_id, "RegionSizeDist", 1,             &dims,                         size_dist);
//            H5LTset_attribute_int(parent_group_id , "RegionSizeDist", "nbins",       &size_dist_nbins,              1);
//            H5LTset_attribute_int(parent_group_id , "RegionSizeDist", "region_type", &tocf_params.region_type_flag, 1);
//            
//            free(size_dist);   // allocated within gen_size_dist
//            
//            SID_log("...done", SID_LOG_CLOSE);   // gen_size_dist
//        }
        
        // ----------------------------------------------------------------------------------------------------
            
        SID_free(SID_FARG grid);
        H5Gclose(group_id);
        
        SID_log("...done", SID_LOG_CLOSE);   // Saving tocf grids
    }
    
}


void check_if_reionization_complete(run_globals_t *run_globals)
{
    SID_log("Checking if reionization complete... (global_xH = %.2f)", SID_LOG_COMMENT, run_globals->tocf_grids.global_xH);
    
    // If the global_xH value is less than 1%, stop doing reionisation
    if (run_globals->tocf_grids.global_xH < 0.01)
        run_globals->params.TOCF_Flag = 0;
}


void call_gen_size_dist(run_globals_t *run_globals, int snapshot, int region_type_flag, float **size_dist, int *size_dist_nbins)
{
    // Thin wrapper around gen_size_dist
    
    float nf;
    tocf_grids_t *grids = &(run_globals->tocf_grids);
    
    // Check if the ionisation field is fully neutral or ionised
    nf = run_globals->tocf_grids.global_xH;
    
    if ((nf==0) || (nf==1))
    {
        SID_log("call_gen_size_dist: The ionisation field is only a single phase... skipping", SID_LOG_CLOSE);
        return;
    }
    
    gen_size_dist(run_globals->ZZ[snapshot], region_type_flag, grids->xH, size_dist, size_dist_nbins);
    
    //SID_log("address of size_dist_nbins (within call_gen_size_dist) = %p", SID_LOG_COMMENT, size_dist_nbins);
    SID_log("value of size_dist_nbins (within call_gen_size_dist) = %d", SID_LOG_COMMENT, *size_dist_nbins);
    
    SID_log("...done", SID_LOG_CLOSE);
}

#endif
