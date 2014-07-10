#!/usr/bin/env python

"""Plot all available results for a particular redshift.

Usage: allresults.py <fname> <redshift> [--hubble=<Hubble_h> --ext=<extension>]

Arguments:
    fname      Input Meraxes file
    redshift   Redshift to plot

Options:
    --hubble=<Hubble_h>  Hubble constant [default: 0.702]
    --ext=<extension>    Ouput file extension [default: pdf]
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from ssimpl import meraxes, plotutils
from docopt import docopt
import sfrf_z5, smf_z5, shmr_z5

__author__ = "Simon Mutch"
__date__   = "2014-07-08"

__available_redshifts__ = [0, 5]

def savefig(fig, fname_in, fname_out, ext):
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname_in),
                                "plots/"+fname_out+"."+ext)
    print "Saving to:",output_fname
    plt.savefig(output_fname)


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    h = float(args['--hubble'])
    redshift = float(args['<redshift>'])
    ext = args['--ext']

    # Check that a valid redshift has been requested
    assert redshift in __available_redshifts__,\
            "Available redshifts are {:s}".format(__available_redshifts__)

    # init the plot style
    plotutils.init_style()

    if redshift == 5:
        snap, redshift = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

        props = ("Sfr", "StellarMass", "Type", "FOFMvir", "CentralGal")
        gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                sim_props=True, h=h)
        gals = gals.view(np.recarray)

        # SMF
        fig, ax = plt.subplots(1,1)
        smf_z5.plot(gals, simprops, ax, h)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="upper right", fontsize="small")
        savefig(fig, fname, "smf-z5", ext)

        # SFRF
        fig, ax = plt.subplots(1,1)
        sfrf_z5.plot(gals, simprops, ax, h)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="lower left", fontsize="small")
        savefig(fig, fname, "sfrf-z5", ext)

        # SHMR
        rc = plt.rcParams.copy()
        plotutils.init_style(theme="white_bg")
        fig, ax = plt.subplots(1,1)
        ax.grid('off')
        shmr_z5.plot(gals, simprops, ax, h)
        ax.tick_params(axis="both", which='both', color=plt.rcParams["axes.edgecolor"], length=3)
        ax.legend(loc="upper left", fontsize="small")
        fig.tight_layout()
        savefig(fig, fname, "shmr-z5", ext)
        plt.rcParams = rc
