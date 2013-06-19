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
from samtools.io import read_gals
from samtools.plots import *
from docopt import docopt
from astropy import log

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

    # Start with the SMF
    fig, ax = plt.subplots(1,1)
    print fig.get_facecolor()
    smf(gals, sim_props, ax)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "smf."+fig_format))


if __name__ == '__main__':
    args = docopt(__doc__)
    allresults(args['<master_file>'], output_dir=args['--output'],
               fig_format=args['--format'])
