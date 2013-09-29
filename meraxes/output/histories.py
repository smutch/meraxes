#!/usr/bin/env python

"""Plot histories...

Usage: histories.py <master_file> <id_list> ... [--output=<dir_path> --format=<ext>]

Options:
    --output=<dir_path>   Target directory for output figures [default: ./plots]
    --format=<ext>        Image format (e.g. png, pdf) [default: png]

"""

import os
import sys
__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.join(__script_dir__, "../utils"))

import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from docopt import docopt
from astropy import log
from astropy.utils.console import ProgressBar
import random

from ssimpl.meraxes.io import read_gals

__author__ = "Simon Mutch"
__date__   = "2013/06/25"

def histories(master_file, id_list, snapshot=None, output_dir='./plots/',
               fig_format='png'):

    plots.set_custom_rcParams(plt.rcParams)

    # If the target directory doesn't exist then create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get a list of available snapshots
    with h5.File(master_file, "r") as fin:
        groups = np.array(fin.keys())
        sel = np.array([v.find("Snap") for v in groups])>=0
        snaplist = np.array([g[-3:] for g in groups[sel]], 'i4')

        # Grab the datatype for a galaxy
        gal_dtype = fin['Snap%03d/Galaxies'%(snaplist[0])].dtype

        if id_list[0] == -1:
            # Grab the size of the lowest redshift snapshot
            last_snap_size = fin['Snap%03d/Galaxies'%(snaplist[-1])].size
            id_list = random.sample(fin['Snap%03d/Galaxies'%(snaplist[-1])]['ID'], 10)

    # Allocate necessary arrays
    zlist = np.zeros(snaplist.size)
    gals = np.zeros((len(id_list), snaplist.size), dtype=gal_dtype)
    gals[:]['Type'] = -1

    # Read in the galaxies
    print "Constructing histories..."
    with ProgressBar(snaplist.size) as bar:
        for ii, snap in enumerate(snaplist[::-1]):
            snap_gals, sim_props, descendant_inds = read_gals(master_file,
                                                              snapshot=snap,
                                                              sim_props=True,
                                                              descendant_inds=True,
                                                              verbose=False)
            zlist[ii] = sim_props['Redshift']
            if snap_gals.size<1:
                continue
            for jj, g in enumerate(gals):
                sel = [snap_gals['ID']==id_list[jj],]
                if np.any(sel):
                    g[ii] = snap_gals[sel]
                else:
                    g[ii]['Type'] = -1
            bar.update()

    # Plot the halo mass accretion histories
    print "Plotting..."
    fig, ax = plt.subplots(1,1)
    plots.histories.mvir(gals, zlist, sim_props, ax, ylim=(8,15), xlim=None, bins=20, ls='-', lw=1, color='k', alpha=0.5)
    # plt.legend(loc='lower left', numpoints=1, fontsize='small', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mvir_histories."+fig_format))

    # import IPython; IPython.embed()

if __name__ == '__main__':
    args = docopt(__doc__)
    if args['<id_list>'][0]=='random':
        np.random.seed(789379)
        args['<id_list>'] = np.array([-1,],'i4')
    else:
        args['<id_list>'] = np.asarray(args['<id_list>'], 'i4')

    histories(args['<master_file>'], args['<id_list>'], output_dir=args['--output'],
               fig_format=args['--format'])
