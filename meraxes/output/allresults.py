#!/usr/bin/env python

"""Plot all results...

Usage: allresults.py <master_file> [--output=<dir_path> --format=<ext>]

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

from samtools.io import read_gals
from samtools.plots import *

__author__ = "Simon Mutch"
__date__   = "2013/06/19"

def allresults(master_file, snapshot=None, output_dir='./plots/',
               fig_format='png'):

    set_custom_rcParams(plt.rcParams)

    # If the target directory doesn't exist then create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read in the galaxies
    gals, sim_props = read_gals(master_file, snapshot=snapshot, sim_props=True)
    redshift = sim_props['Redshift']
    hubble_h = sim_props['Hubble_h']

    # Start with the SMF
    fig, ax = plt.subplots(1,1)
    smf(gals, sim_props, ax)

    # Overplot relevant observational data
    # PG08
    if ((redshift>0.8) & (redshift<=1.0)):
        # Read in and prepare the data
        with h5.File(os.path.join(observations_dir, "SMFs-PG08.hdf5"), "r") as fin:
            obs = fin["Table"][:]
        obs["logM"] += 2.*np.log10(0.7) - 2.*np.log10(hubble_h) 
        for col in ("logIRAC", "E_logIRAC", "e_logIRAC"):
            obs[col] += 3.*np.log10(hubble_h) - 3.*np.log10(0.7)
        m_err = 10.**obs["logIRAC"] - 10.**(obs["logIRAC"]-obs["e_logIRAC"])
        p_err = 10.**(obs["E_logIRAC"]+obs["logIRAC"]) - 10.**obs["logIRAC"]

        # Plot it
        sel = (obs["Range"]=="0.8<z<1.0")
        ax.errorbar(obs[sel]["logM"], 10.**obs[sel]["logIRAC"], yerr=(m_err[sel], p_err[sel]),
                    ls='none', marker='s', capthick=2, mec='none',
                    label="PG08", zorder=0)


    plt.legend(loc='lower left', numpoints=1, fontsize='small', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "smf."+fig_format))

    # Plot the halo mass function
    fig, ax = plt.subplots(1,1)
    halo_mf(gals, sim_props, ax)

    plt.legend(loc='lower left', numpoints=1, fontsize='small', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "halo_mf."+fig_format))


if __name__ == '__main__':
    args = docopt(__doc__)
    allresults(args['<master_file>'], output_dir=args['--output'],
               fig_format=args['--format'])
