#!/usr/bin/env python

"""Plot the supernova feedback parameters as a function of different
properties.

Usage: sn_feedback_params.py <fname> <redshift> [--h=<Hubble_h>]

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

def calc_eff_sn_params(gals, simprops):
    """calculate the SN params"""

    SnReheatParam = simprops['SnReheatEff']\
            * (0.5 + (gals.Vmax/simprops['SnReheatNorm'])**(-simprops['SnReheatScaling']))
    SnReheatParam[SnReheatParam > simprops['SnReheatLimit']] = simprops['SnReheatLimit']

    SnEjectionParam = simprops['SnEjectionEff']\
            * (0.5 + (gals.Vmax/simprops['SnEjectionNorm'])**(-simprops['SnEjectionScaling']))
    SnEjectionParam[SnEjectionParam > 1.0] = 1.0

    return SnReheatParam, SnEjectionParam


def plot(gals, simprops, redshift, axs, h):

    print "Plotting z={:g} SN feedback params...".format(redshift)

    # bin by stellar mass
    nbins = 20
    edges = np.linspace(7,11,nbins+1)
    width = edges[1]-edges[0]
    centers = edges[:-1] + 0.5*width
    ind = np.digitize(np.log10(gals.StellarMass*1e10), edges)-1

    # calc SN params
    SnReheatParam, SnEjectionParam = calc_eff_sn_params(gals, simprops)

    # plot the model
    vals = [SnReheatParam[sel] for sel in [ind==ii for ii in xrange(0, ind.max()+1)]]
    axs[0].boxplot(vals, positions=centers[:ind.max()+1], widths=width*0.75)

    # add some text
    axs[0].text(0.95,0.95, "z={:d}\nh={:.2f}\nSalpeter IMF\n".format(int(redshift), h),
            horizontalalignment="right",
            verticalalignment="top",
            transform=axs[0].transAxes)

    axs[0].set_ylim([0,simprops['SnReheatLimit']+2])

    plt.setp(axs[0].get_xticklabels(), visible=False)

    axs[0].set_ylabel(r"SnReheatEff")

    vals = [SnEjectionParam[sel] for sel in [ind==ii for ii in xrange(0, ind.max()+1)]]
    axs[1].boxplot(vals, positions=centers[:ind.max()+1], widths=width*0.75)

    axs[1].set_ylim([0,1.2])

    plt.setp(axs[1].get_xticklabels(), rotation=45)

    axs[1].set_xlabel(r"$\log_{10}(M_{*}/\rm M_{\odot})$")
    axs[1].set_ylabel(r"SnEjectionEff")


def plot_vmax(gals, simprops, redshift, ax, h):

    print "Plotting z={:g} Vmax as func of M*...".format(redshift)

    # bin by stellar mass
    nbins = 20
    edges = np.linspace(7,11,nbins+1)
    width = edges[1]-edges[0]
    centers = edges[:-1] + 0.5*width
    ind = np.digitize(np.log10(gals.StellarMass*1e10), edges)-1

    # plot the model
    vals = [gals.Vmax[sel] for sel in [ind==ii for ii in xrange(0, ind.max()+1)]]
    ax.boxplot(vals, positions=centers[:ind.max()+1], widths=width*0.75)

    # add some text
    ax.text(0.05,0.95, "z={:d}\nh={:.2f}\nSalpeter IMF\n".format(int(redshift), h),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes)

    # ax.set_ylim([0,])
    plt.setp(ax.get_xticklabels(), rotation=45)

    ax.set_ylabel(r"$V_{\rm max}\ [{\rm km/s}]$")
    ax.set_xlabel(r"$\log_{10}(M_{*}/\rm M_{\odot})$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    redshift = float(args['<redshift>'])
    h = float(args['--h'])

    snap, redshift = meraxes.io.check_for_redshift(fname, redshift, tol=0.1)

    props = ("StellarMass", "Vmax",)
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                                          sim_props=True, h=h)
    gals = gals.view(np.recarray)

    fig, axs = plt.subplots(2, 1, sharex=True)
    plt.subplots_adjust(hspace=0.02)
    plot(gals, simprops, redshift, axs, h)
    # ax.yaxis.set_tick_params(which='both', color='w')
    # axs[1].legend(loc="upper left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/sn_feedback_params-z{:d}.png".format(int(redshift)))
    plt.savefig(output_fname)

    fig, ax = plt.subplots(1, 1)
    plot_vmax(gals, simprops, redshift, ax, h)
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname),
                                "plots/Vmax_vs_stellarmass-z{:d}.png".format(int(redshift)))
    plt.savefig(output_fname)
