#!/usr/bin/env python

"""Plot the global star formation rate history.

Usage: sfh.py <fname> [Hubble_h]

Arguments:
    Hubble_h   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from ssimpl import meraxes
from astropy.utils.console import ProgressBar
import os

__author__ = "Simon Mutch"
__date__ = "2014-08-13"

__script_dir__ = os.path.dirname(os.path.realpath(__file__))


def plot(fname, ax, h):

    props = ("Sfr", "GhostFlag")
    simprops = meraxes.io.read_input_params(fname, h)
    snaps, redshifts, lt_times = meraxes.io.read_snaplist(fname, h)
    volume = simprops["Volume"]

    sfrd = np.zeros(snaps.shape)

    # generate the model sfh
    print "Generating SFH..."
    with ProgressBar(snaps.shape[0]) as bar:
        for snap in snaps:
            try:
                gals = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                                            h=h, quiet=True)
            except IndexError:
                continue
            gals = gals.view(np.recarray)
            gals = np.compress(gals.GhostFlag == 0, gals)

            sfrd[snap] = gals.Sfr.sum() / volume
            bar.update()

    print "Plotting star formation history..."

    # plot the model
    l, = ax.plot(redshifts, np.log10(sfrd), label="Meraxes")

    # add some text
    ax.text(0.95, 0.95, "h={:.2f}\nSalpeter IMF\n".format(h),
            # r"log$_{10}$(SFR)"+" > {:.2f}".format(sfr_limit)
            # +r"M$_{\odot}$/yr",
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([5, 25])
    ax.set_ylim([-4, -1])

    ax.set_xlabel(r"redshift")
    ax.set_ylabel(r"$\log_{10}\psi ({\rm M_{\odot}yr^{-1}Mpc^{-3}})$")


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
    output_fname = os.path.join(os.path.dirname(fname), "plots/sfh.pdf")
    plt.savefig(output_fname)
