#!/usr/bin/env python

"""Plot the z=0 SMF.

Usage: smf-z0.py <fname>

"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from ssimpl import plotutils, meraxes, munge
import pandas as pd
from astropy import log
import os

__author__ = "Simon Mutch"
__date__   = "2014-05-16"

__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))


def plot_smf_z0(gals, simprops, ax):

    print "Plotting z=0 SMF..."

    # read the observed smf
    obs = pd.read_table(os.path.join(__script_dir__,
                                     "../../utils/obs_datasets/smf/Baldry2008-total.txt"),
                        delim_whitespace=True,
                        header = None,
                        skiprows = 10,
                        usecols = range(1,7),
                        names = ["sm", "n_gal", "phi", "poisson", "min_phi", "max_phi"])

    # convert obs to same hubble value
    obs.sm += np.log10(0.7/simprops["Hubble_h"])
    for col in ["phi", "min_phi", "max_phi"]:
        obs[col] /= (0.7/simprops["Hubble_h"])**3

    # sm = stellar mass
    sm = np.log10(gals.StellarMass[gals.StellarMass > 0] * 1.0e10)

    n_dropped = gals.shape[0] - sm.shape[0]
    if n_dropped > 0:
        log.warn("Dropped %d galaxies (%.1f%% of total) with stellar mass <= 0" %
                (n_dropped, float(n_dropped)/gals.shape[0]*100))

    smf = munge.mass_function(sm, simprops["Volume"], "knuth")

    # plot the model
    ax.plot(smf[:,0], np.log10(smf[:,1]), label="Meraxes")

    # # plot the hmf
    # hm = np.log10(gals.Mvir[gals.Type == 0] * 1.0e10)
    # hmf = munge.mass_function(hm, simprops["Volume"], 50)
    # ax.plot(np.log10(0.05*0.17)+hmf[:,0], np.log10(hmf[:,1]), label="HMF (Mvir*0.05*0.17)", ls='--', color='k')

    # plot the observations
    ax.errorbar(obs.sm, np.log10(obs.phi), yerr=[np.log10(obs.phi) - np.log10(obs.min_phi),
                                                 np.log10(obs.max_phi) - np.log10(obs.phi)],
            label="Baldry et al. 2008",
            ls="none",)

    ax.axis([7,12.5, -4, -1])

    ax.set_xlabel(r"$\log_{10}(M_* / {\rm M_{\odot}})$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")

    ax.text(0.05, 0.05,
           "z=0\nh={:.2f}\nChabrier IMF".format(simprops["Hubble_h"]),
           verticalalignment="bottom",
           horizontalalignment="left",
           transform=ax.transAxes)


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    snap, redshift = meraxes.io.check_for_redshift(fname, 0.0, tol=0.1)

    props = ("StellarMass", "Mvir", "Type")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot_smf_z0(gals, simprops, ax)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="upper right")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname), "plots/smf-z0.pdf")
    plt.savefig(output_fname)
