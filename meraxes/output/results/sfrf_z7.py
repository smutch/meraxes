#!/usr/bin/env python

"""Plot the z=7 SFRF.

Usage: sfrf_z7.py <fname> [Hubble_h]

Arguments:
    Hubble_h   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from dragons import plotutils, meraxes, munge
import pandas as pd
from astropy import log
import os

__author__ = "Simon Mutch"
__date__   = "2014-07-04"

__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))


def plot(gals, simprops, ax, h):

    print "Plotting z=7 SFRF..."

    markers = ['o', 's', 'D', 'p', 'v', '^'][::-1]

    # generate the model sfrf
    sfr = np.log10(gals.Sfr[(gals.Sfr > 0) & (gals.GhostFlag == 0)])

    n_dropped = gals.shape[0] - sfr.shape[0]
    if n_dropped > 0:
        log.warn("Dropped %d galaxies (%.1f%% of total) with sfr <= 0" %
                (n_dropped, float(n_dropped)/gals.shape[0]*100))

    sfrf = munge.mass_function(sfr, simprops["Volume"], "knuth",
                               poisson_uncert=True)

    # plot the model
    l, = ax.plot(sfrf[:,0], np.log10(sfrf[:,1]), lw=4, label="Meraxes")
    yerr = [np.log10(sfrf[:,1]-sfrf[:,2]), np.log10(sfrf[:,1]+sfrf[:,2])]
    yerr[0][np.isinf(yerr[0])] = -7
    ax.fill_between(sfrf[:,0], yerr[0], yerr[1], color=l.get_color(), alpha=0.3)

    # read the observed sfrf
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/sfrf/Katsianis_z7_sfrf-smit.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 4,
        names = ["log_sfr", "log_phi", "err"])

    # convert obs to same hubble value
    obs.log_sfr += 2.0*np.log10(0.7/h)
    for col in ["log_phi", "err"]:
        obs[col] -= 3.0*np.log10(0.7/h)

    # plot the observations
    ax.errorbar(obs.log_sfr, obs.log_phi, yerr=obs.err,
                label="Katsianis+ 2014 (from Smit+ 2012)", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # do it all again for the next set of observations
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/sfrf/Katsianis_z7_sfrf-bouwens.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 4,
        names = ["log_sfr", "log_phi", "err"])

    # convert obs to same hubble value
    obs.log_sfr += 2.0*np.log10(0.7/h)
    for col in ["log_phi", "err"]:
        obs[col] -= 3.0*np.log10(0.7/h)

    # plot the observations
    ax.errorbar(obs.log_sfr, obs.log_phi, yerr=obs.err,
                label="Katsianis+ 2014 (from Bouwens+ 2014)", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # and again
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/sfrf/Duncan14_SFRF_SED_z7.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 8,
        names = ["log_sfr", "phi", "merr", "perr"])
    obs.merr[obs.merr >= obs.phi] = obs.phi - 1e-10

    # convert obs to same hubble value and IMF
    obs.log_sfr += 2.0*np.log10(0.702/h)
    obs.log_sfr += 0.25  # IMF correction Chabrier -> Salpeter
    for col in ["phi", "merr", "perr"]:
        obs[col] /= (0.7**3/h**3)

    # plot the observations
    ax.errorbar(obs.log_sfr, np.log10(obs.phi),
                yerr=[np.log10(obs.phi / (obs.phi - obs.merr)),
                      np.log10(1.0 + (obs.perr / obs.phi))],
                label="Duncan et al. 2014 (SEDs)", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # and again
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/sfrf/Duncan14_SFRF_UV_z7.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 8,
        names = ["log_sfr", "phi", "merr", "perr"])
    obs.merr[obs.merr >= obs.phi] = obs.phi - 1e-10

    # convert obs to same hubble value and IMF
    obs.log_sfr += 2.0*np.log10(0.702/h)
    obs.log_sfr += 0.25  # IMF correction Chabrier -> Salpeter
    for col in ["phi", "merr", "perr"]:
        obs[col] /= (0.7**3/h**3)

    # plot the observations
    ax.errorbar(obs.log_sfr, np.log10(obs.phi),
                yerr=[np.log10(obs.phi / (obs.phi - obs.merr)),
                      np.log10(1.0 + (obs.perr / obs.phi))],
                label="Duncan et al. 2014 (UV)", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # add some text
    ax.text(0.95,0.95, "z=7\nh={:.2f}".format(h),
           horizontalalignment="right",
           verticalalignment="top",
           transform=ax.transAxes)

    ax.set_xlim([-1,3])
    ax.set_ylim([-6,-1])

    ax.set_xlabel(r"$\log_{10}({\rm SFR} / ({\rm M_{\odot}} {\rm yr^{-1}}))$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    if args['Hubble_h'] is False:
        h = 0.7
    else:
        h = float(args['Hubble_h'])
    meraxes.set_little_h(h)

    snap, redshift = meraxes.io.check_for_redshift(fname, 7.0, tol=0.1)

    props = ("Sfr", "Type", "GhostFlag")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot(gals, simprops, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="lower left", fontsize="small")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname), "plots/sfrf-z7.pdf")
    plt.savefig(output_fname)
