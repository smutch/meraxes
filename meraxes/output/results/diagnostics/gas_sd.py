#!/usr/bin/env python

"""Plot the gas surface density distribution.

Usage: cold_gas_fraction.py <fname> <redshift> [--h=<Hubble_h>]

Arguments:
    redshift   redshift to plot [default: last available]

Options:
    --h=<Hubble_h>   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from ssimpl import meraxes, munge
import os

__author__ = "Simon Mutch"
__date__ = "2014-10-07"

args = docopt(__doc__)


def plot_sd(gals, simprops, redshift, ax, h):

    print "Plotting z={:g} gas surface density distrib...".format(redshift)

    # generate the results
    sel = gals.ColdGas > 0
    sd = np.log10(np.compress(sel, gals.ColdGas*1e10 / (np.pi * (3.0*gals.DiskScaleLength*1000.)**2)))

    mf = munge.mass_function(sd, simprops["Volume"], 50, poisson_uncert=True)

    # plot the model
    l, = ax.plot(mf[:,0], np.log10(mf[:,1]), ls="-", lw=4, label="Meraxes")
    yerr = [np.log10(mf[:,1]-mf[:,2]), np.log10(mf[:,1]+mf[:,2])]
    yerr[0][np.isinf(yerr[0])] = -10
    ax.fill_between(mf[:,0], yerr[0], yerr[1], color=l.get_color(), alpha=0.3)

    # add some text
    ax.text(0.95,0.95, "z={:d}\nh={:.2f}\nSalpeter IMF\n".format(int(redshift), h),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([0,12])
    ax.set_ylim([-6,0])

    ax.set_xlabel(r"$\log_{10}(\Sigma_{\rm gas}/[{\rm M_{\odot}/kpc^2}])$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


def plot_ds(gals, simprops, redshift, ax, h):

    print "Plotting z={:g} disk size distrib...".format(redshift)

    for gtype in range(3):
        # generate the results
        sel = (gals.ColdGas > 0) & (gals.Type == gtype)
        ds = np.log10(np.compress(sel, gals.DiskScaleLength*3*1000.))

        mf = munge.mass_function(ds, simprops["Volume"], 50, poisson_uncert=True)

        # plot the model
        l, = ax.plot(mf[:,0], np.log10(mf[:,1]), ls="-", lw=4, label="Type {:d}".format(gtype))
        yerr = [np.log10(mf[:,1]-mf[:,2]), np.log10(mf[:,1]+mf[:,2])]
        yerr[0][np.isinf(yerr[0])] = -10
        ax.fill_between(mf[:,0], yerr[0], yerr[1], color=l.get_color(), alpha=0.3)

    # add some text
    ax.text(0.95,0.95, "z={:d}\nh={:.2f}\nSalpeter IMF\n".format(int(redshift), h),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([-2.5,2])
    ax.set_ylim([-6,1])

    ax.set_xlabel(r"$\log_{10}(r_{\rm disk}/{\rm kpc})$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    h = float(args['--h'])

    snap, redshift = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

    props = ("StellarMass", "ColdGas", "DiskScaleLength", "Type")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True, h=h)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot_sd(gals, simprops, redshift, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="upper left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/gas_sd-z{:d}.png".format(int(redshift)))
    plt.savefig(output_fname)

    fig, ax = plt.subplots(1,1)
    plot_ds(gals, simprops, redshift, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="upper left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/disk_size-z{:d}.png".format(int(redshift)))
    plt.savefig(output_fname)

