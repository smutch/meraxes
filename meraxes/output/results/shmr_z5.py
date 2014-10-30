#!/usr/bin/env python

"""Plot the z=5 SHMR.

Usage: shmr_z5.py <fname> [Hubble_h]

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
__date__   = "2014-07-08"

__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))


def plot(gals, simprops, ax, h):

    print "Plotting z=5 SHMR..."

    sample_frac = 0.25
    axis_lim = (8,12,-6,0)

    # generate the model stellar--halo mass fractions
    sel = (gals.StellarMass > 0) & (gals.Type == 0)
    mvir = gals.FOFMvir[sel]

    stars = np.zeros(gals.shape[0], 'f')
    sel_valid = gals.CentralGal >= 0
    for c, sm in gals[sel_valid][["CentralGal", "StellarMass"]]:
        stars[c] += sm
    stars = stars[sel]
    shm = stars / mvir
    mvir = np.log10(mvir * 1.e10)

    # randomly sample
    sel = np.random.choice(np.arange(shm.shape[0]), int(sample_frac*shm.shape[0]))
    shm = shm[sel]
    mvir = mvir[sel]

    # plot the data as a hexbin plot
    im = ax.hexbin(mvir, np.log10(shm), bins="log", cmap=plt.cm.gist_earth_r, zorder=-1)
    cbar = plt.colorbar(im, ax=ax, label="Meraxes galaxies")
    cbar.solids.set_edgecolor("face")
    log.info("done hexbin")

    # Calculate the bin edges
    dm, bins = munge.scotts_bin_width(mvir, True)

    # find which bin each galaxy is in and work out the median shm in each one
    ind = np.digitize(mvir, bins) - 1
    med_shm = [np.median(shm[ind==ii]) for ii in xrange(len(bins)-1)]

    # plot the median line
    x_vals = (bins[1] - bins[0])/2.0 + bins[:-1]
    ax.plot(x_vals, np.log10(med_shm), label="Meraxes median")

    # read the simulation data of Tescari et al. 2014
    # hm = halo mass
    # sm = stellar mass
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/shmr/Tescari_z5_sim_shmr.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 5,
        names = ["hm", "sm"])

    # convert values to same hubble value and IMF
    obs["sm"] /= h
    obs["hm"] /= h

    # plot the values
    pc = ax.scatter(np.log10(obs.hm*1.e10), np.log10(obs.sm/obs.hm),
               label="Tescari et al. 2014",
               marker='o', s=7,
               color="0.5",
               alpha=0.3)
    pc.set_rasterized(True)

    hm = np.log10(obs.hm*1.e10)
    dm, bins = munge.scotts_bin_width(hm, True)
    ind = np.digitize(hm, bins) - 1
    shm = obs.sm/obs.hm
    med_shm = [np.median(shm[ind==ii]) for ii in xrange(len(bins)-1)]
    x_vals = (bins[1] - bins[0])/2.0 + bins[:-1]
    ax.plot(x_vals, np.log10(med_shm), color="0.1", label="Tescari et al. (2014) median")

    # # do it all again for the next set of observations
    # obs = pd.read_table(os.path.join(__script_dir__,
    #     "../../utils/obs_datasets/smf/Katsianis_z5_smf.txt"),
    #     delim_whitespace=True,
    #     header = None,
    #     skiprows = 3,
    #     names = ["sm", "log_phi", "err"])

    # # convert obs to same hubble value and IMF
    # # -0.16 dex conversion from Salpeter to Chabrier
    # obs.sm += np.log10(0.702/h)
    # for col in ["log_phi", "err"]:
    #     obs[col] -= 3.0*np.log10(0.702/h)

    # # plot the observations
    # ax.errorbar(obs.sm, obs.log_phi, yerr=obs.err,
    #             label="Katsianis et al. 2014", ls="None",
    #             lw=4, capsize=0)

    # add some text
    ax.text(0.95,0.05, "z=5\nh={:.2f}\nSalpeter IMF".format(h),
           horizontalalignment="right",
           verticalalignment="bottom",
           transform=ax.transAxes)

    ax.axis(axis_lim)

    ax.set_xlabel(r"$\log_{10}(M_{\rm vir} / {\rm M_{\odot}})$")
    ax.set_ylabel(r"$\log_{10}(M_* / M_{\rm vir})$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    if args['Hubble_h'] is False:
        h = 0.702
    else:
        h = float(args['Hubble_h'])

    snap, redshift = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

    props = ("StellarMass", "FOFMvir", "Type", "CentralGal")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True, h=h)
    gals = gals.view(np.recarray)

    plotutils.init_style(theme="white_bg")

    fig, ax = plt.subplots(1,1)
    ax.grid('off')
    plot(gals, simprops, ax, h)
    ax.tick_params(axis="both", which='both', color=plt.rcParams["axes.edgecolor"], length=3)
    ax.legend(loc="upper left", fontsize="small")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname), "plots/shmr-z5.pdf")
    plt.savefig(output_fname)
