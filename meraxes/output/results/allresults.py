#!/usr/bin/env python

"""Plot all available results for a particular redshift.

Usage: allresults.py <fname> [--z=<redshift> --hubble=<Hubble_h> --ext=<extension>]

Arguments:
    fname      Input Meraxes file

Options:
    --z=<redshift>       Redshift to plot (unspecified = all avail.)
    --hubble=<Hubble_h>  Hubble constant (default: sim defined)
    --ext=<extension>    Ouput file extension [default: pdf]
"""

import os
import numpy as np

from savefig import monkey_patch
monkey_patch()
import matplotlib.pyplot as plt

from dragons import meraxes, plotutils
from docopt import docopt

from smf import Smf
import sfrf_z5, smf_z5, shmr_z5
import sfrf_z6, smf_z6
import sfrf_z7, smf_z7
import sfh

__author__ = "Simon Mutch"
__date__   = "2014-07-08"

__available_redshifts__ = [5,6,7]


def savefig(fig, fname_in, fname_out, ext, simprops):
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname_in),
                                "plots/"+fname_out+"."+ext)
    print "Saving to:",output_fname

    # ensure output directory exists
    dirname = os.path.dirname(output_fname)
    if not os.path.exists(dirname):
        os.mkdir(dirname)

    plt.savefig(output_fname, extra_info=simprops)


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    ext = args['--ext']
    redshift = args['--z']

    h = args['--hubble']
    if h is None:
        h = meraxes.read_input_params(fname)["Hubble_h"]
    else:
        h = float(h)

    meraxes.set_little_h(h)

    if redshift:
        redshift = float(redshift)
        # Check that a valid redshift has been requested
        assert redshift in __available_redshifts__,\
                "Available redshifts are {:s}".format(__available_redshifts__)

    # init the plot style
    plotutils.init_style()
    plt.rcParams["savefig.dpi"] = 120

    props = ("Sfr", "StellarMass", "Type", "FOFMvir", "CentralGal",
             "GhostFlag")

    if not redshift or (redshift == 5):
        snap, _ = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

        gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                sim_props=True)
        gals = gals.view(np.recarray)
        simprops["h"] = h

        # SMF
        fig, ax = plt.subplots(1,1)
        model = Smf(fname, -1, snapshot=snap)
        smf_z5.plot(model, ax)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="upper right", fontsize="small")
        savefig(fig, fname, "smf-z5", ext, simprops)

        # SFRF
        fig, ax = plt.subplots(1,1)
        sfrf_z5.plot(gals, simprops, ax, h)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="lower left", fontsize="small")
        savefig(fig, fname, "sfrf-z5", ext, simprops)

        # SHMR
        # rc = plt.rcParams.copy()
        # plotutils.init_style(theme="white_bg")
        # fig, ax = plt.subplots(1,1)
        # ax.grid('off')
        # shmr_z5.plot(gals, simprops, ax, h)
        # ax.tick_params(axis="both", which='both', color=plt.rcParams["axes.edgecolor"], length=3)
        # ax.legend(loc="upper left", fontsize="small")
        # fig.tight_layout()
        # savefig(fig, fname, "shmr-z5", ext, simprops)
        # plt.rcParams = rc


    if not redshift or (redshift == 6):
        snap, _ = meraxes.io.check_for_redshift(fname, 6.0, tol=0.1)

        gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                sim_props=True, h=h)
        gals = gals.view(np.recarray)
        simprops["h"] = h

        # SMF
        fig, ax = plt.subplots(1,1)
        model = Smf(fname, -1, snapshot=snap)
        smf_z6.plot(model, ax)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="upper right", fontsize="small")
        savefig(fig, fname, "smf-z6", ext, simprops)

        # SFRF
        fig, ax = plt.subplots(1,1)
        sfrf_z6.plot(gals, simprops, ax, h)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="lower left", fontsize="small")
        savefig(fig, fname, "sfrf-z6", ext, simprops)


    if not redshift or (redshift == 7):
        snap, _ = meraxes.io.check_for_redshift(fname, 7.0, tol=0.1)

        gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
                sim_props=True, h=h)
        gals = gals.view(np.recarray)
        simprops["h"] = h

        # SMF
        fig, ax = plt.subplots(1,1)
        model = Smf(fname, -1, snapshot=snap)
        smf_z7.plot(model, ax)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="upper right", fontsize="small")
        savefig(fig, fname, "smf-z7", ext, simprops)

        # SFRF
        fig, ax = plt.subplots(1,1)
        sfrf_z7.plot(gals, simprops, ax, h)
        ax.yaxis.set_tick_params(which='both', color='w')
        ax.legend(loc="lower left", fontsize="small")
        savefig(fig, fname, "sfrf-z7", ext, simprops)

    # if not redshift:
    #     # star formation history
    #     fig, ax = plt.subplots(1,1)
    #     sfh.plot(fname, ax, h)
    #     ax.legend(loc="lower left", fontsize="small")
    #     savefig(fig, fname, "sfh", ext, simprops)
