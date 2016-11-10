#!/usr/bin/env python

import numpy as np
from dragons import meraxes, munge
from tqdm import tqdm
import emcee

__author__ = "Simon Mutch"
__date__ = "2015-04-14"


def schechter(m_star, phi_star, alpha, m):
    frac = 10.**m / 10.**m_star
    return np.log(10.0) * 10.**phi_star * frac**(alpha+1) * np.exp(-frac)


def _fit_schechter(m, phi, delta_phi,
                   nwalkers=100, nsamples=1000, nburn=200,
                   theta_init=[10.93, -3.91, -1.74],
                   offset_frac=0.05,
                   m_star_lim=[6.0, 14.0],
                   phi_star_lim=[-6, 1],
                   alpha_lim=[-5, 0],
                   debug_fname=None):

    """ Use emcee to fit a Schechter function to a galaxy SMF.
    N.B. m_star is really log10(m_star)

    Returns marginalised best fit and associated 16 & 18pc confidence
    intervals on parameters.
    """

    def lnprior(theta):
        """Flat priors within specified limits"""
        m_star, phi_star, alpha = theta
        phi_star = 10.**phi_star
        if m_star_lim[0] < m_star < m_star_lim[1]\
           and phi_star_lim[0] < phi_star < phi_star_lim[1]\
           and alpha_lim[0] < alpha < alpha_lim[1]:
            return 0.0

        return -np.inf

    def lnlike(theta, m, p, dp):
        m_star, phi_star, alpha = theta
        model = schechter(m_star, phi_star, alpha, m)
        return -0.5*np.sum((p-model)**2/(dp**2))

    def lnprob(theta, m, p, dp):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, m, p, dp)

    # Initialise the starting positions for our walkers.
    theta_init = np.array(theta_init)
    pos = [theta_init + offset_frac*np.random.randn(3) for i in
           range(nwalkers)]

    # Run emcee
    sampler = emcee.EnsembleSampler(nwalkers, 3, lnprob, a=5.0,
                                    args=(m, phi, delta_phi))
    sampler.run_mcmc(pos, nsamples)

    # Get rid of our burn in, convert back to linear phi_star and then generate
    # marginalised 16 & 84% confidence intervals for each parameter.
    chain = sampler.chain[:, nburn:, :].reshape((-1, 3))

    # DEBUG
    if debug_fname is not None:
        import matplotlib.pyplot as plt
        from triangle import triangle
        fig, axs = plt.subplots(3, 1)
        for ii in xrange(0, nwalkers, 10):
            for ia, ax in enumerate(axs):
                ax.plot(sampler.chain[ii, :, ia], lw=1, color='k')
        fig.savefig("%s-chains.png" % debug_fname)
        fig = triangle.corner(chain, labels=[r"$\log_{10}(M_*)$",
                                             r"$\log_{10}(\phi_*)$",
                                             r"$\alpha$"],
                              truths=theta_init)
        fig.savefig("%s-posteriors.png" % debug_fname)

    m_star, phi_star, alpha = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                  zip(*np.percentile(chain, [16, 50, 84],
                                                     axis=0)))

    # Return the marginalised best fit param values with confidence intervals
    return m_star, phi_star, alpha


class Smf(object):
    """The base class for a smf belonging to a particular redshift or model."""

    def __init__(self, fname, redshift, snapshot=None):
        self.fname = fname
        if snapshot is None:
            self.snapshot, self.redshift =\
                meraxes.check_for_redshift(fname, redshift, tol=0.1)
        else:
            self.snapshot = snapshot
            self.redshift = meraxes.grab_redshift(fname, snapshot)
        self.simprops = meraxes.read_input_params(fname, quiet=True)
        self.volume = self.simprops["Volume"]

        self.min_phi = None
        self.max_phi = None
        self.median_phi = None
        self.masses = None
        self.mf = None

    def set_masses(self, masses):
        self.masses = np.log10(masses*1e10)

    def read_masses(self):
        sm = meraxes.io.read_gals(self.fname, snapshot=self.snapshot,
                                  props=["StellarMass",],
                                  quiet=True)["StellarMass"]
        # drop zeros
        sel = sm > 0
        self.masses = np.log10(sm[sel]*1e10)

    def generate(self, bins="knuth", limits=None):
        if self.masses is None:
            self.read_masses()

        if type(limits) is list:
            if len(limits) != 2:
                raise IndexError("len(limits) must equal 2")
            if limits[0] is None:
                limits[0] = self.masses.min()
            if limits[1] is None:
                limits[1] = self.masses.max()

        self.mf, self.edges =\
            munge.mass_function(self.masses, self.volume, bins,
                                return_edges=True, poisson_uncert=True,
                                range=limits)

    def gen_uncert(self, n_samples=100):
        # sigma--redshift relation taken from Behroozi et al. 2013 (used by Lu
        # et al. 2013)
        sigma = 0.07 + 0.04*self.redshift

        n_gals = np.zeros([n_samples, self.edges.shape[0]-1])
        bin_width = self.edges[1] - self.edges[0]

        for ii in tqdm(xrange(n_samples), desc="Generating uncerts"):
            pert_sm = self.masses + np.random.normal(0, sigma,
                                                     self.masses.shape[0])
            n_gals[ii, :], _ = np.histogram(pert_sm, bins=self.edges)

        # add poisson uncert to each bin count and then find the 68th
        # percentile
        self.max_phi = np.percentile(n_gals+np.sqrt(n_gals), 68.0, axis=0)\
            / (self.volume * bin_width)

        # subtract poisson uncert from each bin count and then find the 32nd
        # percentile
        self.min_phi = np.percentile(n_gals-np.sqrt(n_gals), 32.0, axis=0)\
            / (self.volume * bin_width)

        # find the median
        self.median_phi = np.median(n_gals.astype(float), axis=0)/(self.volume
                                                                   * bin_width)

    def plot(self, ax, **kwargs):
        if self.mf is None:
            self.generate()

        if self.min_phi is not None:
            phi, min_phi, max_phi = [np.log10(v) for v in (self.median_phi,
                                                           self.min_phi,
                                                           self.max_phi)]
            min_phi[np.isinf(min_phi)] = -7

            self.line, = ax.plot(self.mf[:, 0], phi, **kwargs)
            self.patch = ax.fill_between(self.mf[:, 0], min_phi, max_phi,
                                         color=self.line.get_color(),
                                         alpha=0.3)

        else:
            self.line, = ax.plot(self.mf[:, 0], np.log10(self.mf[:, 1]),
                                 **kwargs)

    def set_axlabel(self, ax):
        ax.set_xlabel(r"$\log_{10}(M_* / {\rm M_{\odot}})$")
        ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")

    def fit_schechter(self, limits=None):
        mass = self.mf[:, 0].copy()
        min_phi = self.min_phi.copy()
        max_phi = self.max_phi.copy()

        if limits is not None:
            sel = (mass >= limits[0]) & (mass < limits[1])
            mass = mass[sel]
            min_phi = min_phi[sel]
            max_phi = max_phi[sel]

        delta_phi = 0.5 * (max_phi - min_phi)
        phi = min_phi + delta_phi
        params = _fit_schechter(mass, phi, delta_phi)
        return params
