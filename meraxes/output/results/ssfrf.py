#!/usr/bin/env python

"""Plot the specific star formation rate function for a given redshift.

Usage: ssfrf_z5.py <fname> <redshift> [Hubble_h]

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

    print "Plotting z=5 ssfrf..."

    # generate the model smf
    # sm = stellar mass
    old_count = gals.shape[0]
    ssfr = gals.Sfr / (gals.StellarMass*1.e10)
    ssfr = np.compress((ssfr > 0) & (gals.GhostFlag==0), ssfr)
    new_count = ssfr.shape[0]


    n_dropped = old_count - new_count
    if n_dropped > 0:
        log.warn("Dropped %d galaxies (%.1f%% of total) with sfr <= 0" %
                (n_dropped, float(n_dropped)/new_count*100))

    lower_lim = -12.0
    ssfr = np.log10(ssfr)
    ssfr = np.compress(ssfr > lower_lim, ssfr)
    ssfrf = munge.mass_function(ssfr, simprops["Volume"], "scotts")

    # plot the model
    ax.plot(ssfrf[:,0], np.log10(ssfrf[:,1]), label="Meraxes")

    ax.set_xlim([lower_lim, -7])
    ax.set_ylim([-4, -0.5])

    ax.set_xlabel(r"$\log_{10}({\rm sSFR} / ({\rm yr^{-1}}))$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    if args['Hubble_h'] is False:
        h = 0.702
    else:
        h = float(args['Hubble_h'])
    redshift = float(args['<redshift>'])

    snap, redshift = meraxes.io.check_for_redshift(fname, redshift, tol=0.1)

    props = ("Sfr", "StellarMass", "GhostFlag")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True, h=h)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot(gals, simprops, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="lower left", fontsize="small")

    # add some text
    ax.text(0.95,0.95, "z=5\nh={:.2f}".format(h),
           horizontalalignment="right",
           verticalalignment="top",
           transform=ax.transAxes)

    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/ssfrf-z{:d}.pdf".format(int(redshift)))
    plt.savefig(output_fname)
