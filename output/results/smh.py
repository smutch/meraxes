#!/usr/bin/env python

"""Plot the global stellar mass density history.

Usage: smh.py <fname> [Hubble_h]

Arguments:
    Hubble_h   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from dragons import meraxes
from astropy.utils.console import ProgressBar
import os

__author__ = "Simon Mutch"
__date__ = "2014-08-13"

__script_dir__ = os.path.dirname(os.path.realpath(__file__))


def plot(fname, ax, h):

    props = ("StellarMass", "GhostFlag")
    simprops = meraxes.io.read_input_params(fname, h)
    snaps, redshifts, lt_times = meraxes.io.read_snaplist(fname, h)
    volume = simprops["Volume"]

    # Stellar mass limit
    sm_lim = simprops["PartMass"]*100*simprops["BaryonFrac"]*0.1

    smd = np.zeros(snaps.shape)

    # generate the model smh
    print "Generating smh..."
    with ProgressBar(snaps.shape[0]) as bar:
        for snap in snaps:
            try:
                gals = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                                            h=h, quiet=True)
            except IndexError:
                continue
            gals = np.compress((gals['GhostFlag'] == 0) &
                               (gals['StellarMass'] >= sm_lim), gals)

            smd[snap] = gals['StellarMass'].sum()
            bar.update()

    smd /= volume

    print "Plotting star formation history..."

    # plot the model
    l, = ax.plot(redshifts, np.log10(smd), label="Meraxes")

    # add some text
    ax.text(0.95, 0.95, "h={:.2f}\nSalpeter IMF\n".format(h)+
            r"log$_{10}(M/M_{\odot}$)"+" > {:.2e}".format(sm_lim*1e10)+
            r"M$_{\odot}$",
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([5, 25])
    ax.set_ylim([-7, -2])

    ax.set_xlabel(r"redshift")
    ax.set_ylabel(r"$\log_{10}\psi ({\rm 10^{10}M_{\odot}\, Mpc^{-3}})$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    if args['Hubble_h'] is False:
        h = 0.702
    else:
        h = float(args['Hubble_h'])

    fig, ax = plt.subplots(1, 1)
    plot(fname, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="lower left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname), "plots/smh.pdf")
    plt.savefig(output_fname)
