import numpy as np
from tqdm import tqdm
from dragons import munge

import emcee

def smf_uncert(sm, volume, edges, redshift, n_samples=20):

    # sigma--redshift relation taken from Behroozi et al. 2013 (used by Lu et
    # al. 2013)
    sigma = 0.07 + 0.04*redshift

    min_phi = np.ones(edges.shape[0]-1)*1e10
    max_phi = np.zeros(edges.shape[0]-1)

    for ii in tqdm(xrange(n_samples), desc="Generating uncerts"):
        pert_sm = sm + np.random.normal(0, sigma, sm.shape[0])
        mf = munge.mass_function(pert_sm, volume, edges, poisson_uncert=True)

        phi = mf[:,1]-mf[:,2]
        sel = phi<min_phi
        min_phi[sel] = phi[sel]

        phi = mf[:,1]+mf[:,2]
        sel = phi>max_phi
        max_phi[sel] = phi[sel]

    return min_phi, max_phi
