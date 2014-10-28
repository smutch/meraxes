#!/usr/bin/env python

"""Plot the used cold gas fractions.

Usage: cold_gas_fraction.py <fname> <redshift> [--h=<Hubble_h>]

Arguments:
    redshift   redshift to plot [default: last available]

Options:
    --h=<Hubble_h>   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from dragons import meraxes, munge
import os

__author__ = "Simon Mutch"
__date__ = "2014-10-07"

args = docopt(__doc__)

def plot(gals, simprops, redshift, ax, h):

    print "Plotting z={:g} gas maximum fractions used distrib...".format(redshift)

    # generate the results
    sel = gals.MaxReheatFrac > 0
    gals = np.compress(sel, gals)

    # bin by stellar mass
    nbins = 20
    edges = np.linspace(6,11.5,nbins+1)
    width = edges[1]-edges[0]
    centers = edges[:-1] + 0.5*width
    ind = np.digitize(np.log10(gals.StellarMass*1e10), edges)-1

    # plot the model
    vals = [gals.MaxReheatFrac[sel] for sel in [ind==ii for ii in xrange(0, ind.max()+1)]]
    ax.boxplot(vals, positions=centers[:ind.max()+1], widths=width*0.75)

    # add some text
    ax.text(0.95,0.95, "z={:d}\nh={:.2f}\nSalpeter IMF\n".format(int(redshift), h),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([5.5,12])
    # ax.set_ylim([-6,0])

    plt.setp(ax.get_xticklabels(), rotation=45)

    ax.set_xlabel(r"$\log_{10}(M_{*}/\rm M_{\odot})$")
    ax.set_ylabel(r"used cold gas fraction")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    redshift = float(args['<redshift>'])
    h = float(args['--h'])

    snap, redshift = meraxes.io.check_for_redshift(fname, redshift, tol=0.1)

    props = ("StellarMass", "ColdGas", "DiskScaleLength", "Type",
             "MaxReheatFrac", "MaxEjectedFrac")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                                          sim_props=True, h=h,
                                          h_conv={'MaxReheatFrac' : lambda x,
                                                  h: x})
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot(gals, simprops, redshift, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="upper left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/gas_maxfracs-z{:d}.png".format(int(redshift)))
    plt.savefig(output_fname)
