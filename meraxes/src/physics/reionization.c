#include "meraxes.h"
#include <math.h>

#ifdef USE_TOCF

void calculate_Mvir_crit(run_globals_t *run_globals, double redshift)
{
  // Calculate the critical Mvir value in each grid cell (ala Sobacchi & Mesinger 2013b)
  
  int            HII_dim        = tocf_params.HII_dim;
  float          Mvir_atomic;
  float          cell_Mvir_crit;

  float          m_0_sm         = tocf_params.m_0_sm;
  float          a_sm           = tocf_params.a_sm;
  float          b_sm           = tocf_params.b_sm;
  float          c_sm           = tocf_params.c_sm;
  float          d_sm           = tocf_params.d_sm;

  float         *Mvir_crit      = run_globals->tocf_grids.Mvir_crit;
  float         *J_21_at_ion    = run_globals->tocf_grids.J_21_at_ionization;
  float         *z_at_ion       = run_globals->tocf_grids.z_at_ionization;

  // init
  memset(Mvir_crit, 0, sizeof(float)*HII_TOT_NUM_PIXELS);
  Mvir_atomic = convert_Tvir_to_Mvir(redshift, tocf_params.ion_tvir_min);

  // Loop through each cell and calculate the value of Mvir_crit
  for(int ii=0; ii<HII_dim; ii++)
  {
    for(int jj=0; jj<HII_dim; jj++)
    {
      for(int kk=0; kk<HII_dim; kk++)
      {
        // Initialise critical mass to atomic cooling mass
        cell_Mvir_crit = Mvir_atomic;
        
        // If this cell was ionized in the past then calculate the critical
        // mass using the UVB feedback prescription of Sobacchi & Mesinger
        // 2013b
        if(z_at_ion[HII_R_INDEX(ii,jj,kk)] > redshift)
        {
          cell_Mvir_crit = m_0_sm*pow((1.0+redshift)/10.0, a_sm) * pow(J_21_at_ion[HII_R_INDEX(ii,jj,kk)], b_sm)*
            pow((1.0-pow((1.0+redshift)/(1.0+z_at_ion[HII_R_INDEX(ii,jj,kk)]), c_sm)), d_sm);

          // DEBUG
          // SID_log("Cell was ionized in past -> Mvir_atomic=%.2e, cell_Mvir_crit=%.2e", SID_LOG_COMMENT, Mvir_atomic, cell_Mvir_crit);
          // SID_log("\tJ_21_at_ion=%.2e, z_at_ion=%.2e, redshift=%.2e, m_0_sm=%.2e, a_sm=%.2e, b_sm=%.2e, c_sm=%.2e, d_sm=%.2e", SID_LOG_COMMENT, 
          //     J_21_at_ion[HII_R_INDEX(ii,jj,kk)], z_at_ion[HII_R_INDEX(ii,jj,kk)], redshift, m_0_sm, a_sm, b_sm, c_sm, d_sm);
        }

        // Save the critical mass to the grid
        Mvir_crit[HII_R_INDEX(ii,jj,kk)] = (Mvir_atomic > cell_Mvir_crit) ? Mvir_atomic : cell_Mvir_crit;
      }
    }
  }
}

bool check_reionization_cooling(run_globals_t *run_globals, halo_t *halo)
{

  if(tocf_params.uvb_feedback)
  {
    bool    flag;
    float   Mvir;
    double  box_size    = run_globals->params.BoxSize;
    float  *M_crit_grid = run_globals->tocf_grids.Mvir_crit;

    // Find which cell this halo lies in
    int i = find_cell((halo->Pos)[0], box_size);
    int j = find_cell((halo->Pos)[1], box_size);
    int k = find_cell((halo->Pos)[2], box_size);

    // If the halo virial mass is below the critical for this cell then set the
    // cooling flag to false, else set it to true
    Mvir = halo->Mvir*1.e10/run_globals->params.Hubble_h;
    flag = (Mvir < M_crit_grid[HII_R_INDEX(i,j,k)]) ? false : true;

    // DEBUG
    // SID_log("Mvir=%.2e, M_crit_grid=%.2e, cooling_flag=%d", SID_LOG_COMMENT, Mvir, M_crit_grid[HII_R_INDEX(i,j,k)], flag);

    return flag;
  } else
  {
    return true;
  }

}

#endif
