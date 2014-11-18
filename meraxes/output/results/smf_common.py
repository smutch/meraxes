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

        phi = mf[:, 1]-mf[:, 2]
        sel = phi < min_phi
        min_phi[sel] = phi[sel]

        phi = mf[:, 1]+mf[:, 2]
        sel = phi > max_phi
        max_phi[sel] = phi[sel]

    return min_phi, max_phi


def fit_schechter(m, phi, min_phi, max_phi,
                  nwalkers=100, nsamples=1000, nburn=200,
                  theta_init=[10.68, 1.24e-4, -1.74],
                  offset_frac=0.01,
                  mstar_lim=[9.0, 12.0],
                  phi_star_lim=[6e-6, 6e-2],
                  alpha_lim=[-3, -0.5]):

    """ Use emcee to fit a Schechter function to a galaxy SMF.
    N.B. m_star is really log10(m_star)

    Returns marginalised best fit and associated 16 & 18pc confidence intervals
    on parameters.
    """

    def schechter(m_star, phi_star, alpha, m):
        """Schechter function"""

        frac = 10.**m / 10.**m_star

        return np.log(10.0) * 10.**phi_star * frac**(alpha+1) *\
            np.exp(-frac)

    def lnprior(theta):
        """Flat priors within specified limits"""

        m_star, phi_star, alpha = theta
        phi_star = 10.**phi_star
        if mstar_lim[0] < m_star < mstar_lim[1]\
           and phi_star_lim[0] < phi_star < phi_star_lim[1]\
           and alpha_lim[0] < alpha < alpha_lim[1]:
            return 0.0

        return -np.inf

    def lnprob(theta, m, phi, min_phi, max_phi):
        """Calculate the ln(prob) value."""

        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        m_star, phi_star, alpha = theta
        model = schechter(m_star, phi_star, alpha, m)
        uncert = np.zeros(phi.shape)

        # deal with uneven uncerts
        sel = model <= phi
        uncert[sel] = (phi-min_phi)[sel]
        uncert[~sel] = (max_phi-phi)[~sel]

        return lp - 0.5*np.sum((phi-model)**2 / uncert**2)

    # Initialise the starting positions for our walkers.
    # Note that we are sampling log10(phi_star).
    theta_init[1] = np.log10(theta_init[1])
    theta_init = np.array(theta_init)
    init_pos = [theta_init + offset_frac*theta_init*np.random.randn(3) for i in
                xrange(nwalkers)]

    # Run emcee
    sampler = emcee.EnsembleSampler(nwalkers, 3, lnprob,
                                    args=(m, phi, min_phi, max_phi))
    sampler.run_mcmc(init_pos, nsamples)

    # Get rid of our burn in, convert back to linear phi_star and then generate
    # marginalised 16 & 84% confidence intervals for each parameter.
    chain = sampler.chain[:, nburn:, :].reshape((-1, 3))
    chain[:, 1] = 10.**chain[:, 1]
    m_star, phi_star, alpha = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                  zip(*np.percentile(chain, [16, 50, 84],
                                                     axis=0)))

    # Return the marginalised best fit param values with confidence intervals
    return m_star, phi_star, alpha
