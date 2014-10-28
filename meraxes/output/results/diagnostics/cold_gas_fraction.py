#!/usr/bin/env python

"""Plot the cold gas fraction distribution.

Usage: cold_gas_fraction.py <fname> <redshift> [--h=<Hubble_h>]

Arguments:
    redshift   redshift to plot [default: last available]

Options:
    --h=<Hubble_h>   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from dragons import plotutils, meraxes, munge
import os

__author__ = "Simon Mutch"
__date__ = "2014-10-07"

args = docopt(__doc__)


def plot(gals, simprops, redshift, ax, h):

    print "Plotting z={:g} cold gas fractions...".format(redshift)

    # generate the results
    min_sm = (10.0**7.5)/1e10
    fracs = np.log10(np.compress((gals.StellarMass > min_sm) & (gals.ColdGas > 0),
                                 gals.ColdGas / gals.StellarMass))

    mf = munge.mass_function(fracs, simprops["Volume"], 50,
                              poisson_uncert=True)

    # plot the model
    l, = ax.plot(mf[:,0], np.log10(mf[:,1]), ls="-", lw=4, label="Meraxes")
    yerr = [np.log10(mf[:,1]-mf[:,2]), np.log10(mf[:,1]+mf[:,2])]
    yerr[0][np.isinf(yerr[0])] = -7
    ax.fill_between(mf[:,0], yerr[0], yerr[1], color=l.get_color(), alpha=0.3)

    # add some text
    ax.text(0.95,0.95, "z={:d}\nh={:.2f}\nSalpeter IMF\n".format(int(redshift), h)+
            r"$M_*>$"+"{:.1e}".format(min_sm*1e10)+r"${\rm M_{\odot}}$",
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([-5,3])
    ax.set_ylim([-6,0])

    ax.set_xlabel(r"$\log_{10}$(gas fraction)")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    h = float(args['--h'])

    snap, redshift = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

    props = ("StellarMass", "ColdGas", "Type", "Sfr")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True, h=h)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot(gals, simprops, redshift, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="upper left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/cold_gas_fraction-z{:d}.png".format(int(redshift)))
    plt.savefig(output_fname)
