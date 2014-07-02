#!/usr/bin/env python

"""Plot the z=5 SMF.

Usage: smf-z5.py <fname>

"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from ssimpl import plotutils, meraxes, munge
import pandas as pd
from astropy import log
import os

__author__ = "Simon Mutch"
__date__   = "2014-05-01"

__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))


def plot_smf_z5(gals, simprops, ax):

    print "Plotting z=5 SMF..."

    # generate the model smf
    # sm = stellar mass
    sm = np.log10(gals.StellarMass[gals.StellarMass > 0] * 1.0e10)

    n_dropped = gals.shape[0] - sm.shape[0]
    if n_dropped > 0:
        log.warn("Dropped %d galaxies (%.1f%% of total) with stellar mass <= 0" %
                (n_dropped, float(n_dropped)/gals.shape[0]*100))

    smf = munge.mass_function(sm, simprops["Volume"], "knuth")

    # plot the model
    ax.plot(smf[:,0], np.log10(smf[:,1]), label="Meraxes")

    # read the observed smf
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/smf/Gonzalez11_z5_smf.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 3,
        names = ["sm", "log_phi", "m_err", "p_err"])

    # convert obs to same hubble value
    obs.sm += np.log10(0.7/simprops["Hubble_h"])
    for col in ["log_phi", "m_err", "p_err"]:
        obs[col] -= 3.0*np.log10(0.7/simprops["Hubble_h"])

    # plot the observations
    ax.errorbar(obs.sm, obs.log_phi, yerr=[obs.m_err, obs.p_err],
                label="Gonzalez et al. 2011", ls="none",
                lw=4, capsize=0)

    # do it all again for the next set of observations
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/smf/Katsianis_z5_smf.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 3,
        names = ["sm", "log_phi", "err"])

    # convert obs to same hubble value and IMF
    # -0.26 dex conversion from Salpeter to Chabrier
    obs.sm += np.log10(0.702/simprops["Hubble_h"]) - 0.26
    for col in ["log_phi", "err"]:
        obs[col] -= 3.0*np.log10(0.702/simprops["Hubble_h"])

    # plot the observations
    ax.errorbar(obs.sm, obs.log_phi, yerr=obs.err,
                label="Katsianis et al. 2014", ls="none",
                lw=4, capsize=0)

    # add some text
    ax.text(0.05,0.05, "z=5\nh={:.2f}\nChabrier IMF".format(simprops["Hubble_h"]),
           horizontalalignment="left",
           verticalalignment="bottom",
           transform=ax.transAxes)

    ax.set_xlim([7,11])

    ax.set_xlabel(r"$\log_{10}(M_* / {\rm M_{\odot}})$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    snap, redshift = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

    props = ("StellarMass", "Mvir", "Type")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot_smf_z5(gals, simprops, ax)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="upper right")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname), "plots/smf-z5.pdf")
    plt.savefig(output_fname)
