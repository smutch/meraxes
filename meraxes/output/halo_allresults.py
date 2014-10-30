"""Plot all results for halos...

Usage: allresults.py <master_file> [--snapshot=<snap> --output=<dir_path> --format=<ext>]

Options:
    --output=<dir_path>   Target directory for output figures [default: ./plots]
    --format=<ext>        Image format (e.g. png, pdf) [default: png]
    --snapshot=<snap>     Snapshot to read in [default: -1]

"""

import os
import sys
__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.join(__script_dir__, "../utils"))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc_file(__script_dir__+"/matplotlibrc")
from docopt import docopt
from astropy import log
from astropy.table import Table

from dragons.meraxes.io import read_gals
from dragons import munge

__author__ = "Simon Mutch"
__date__   = "2013/08/08"

def Rvir_distrib(gals, redshift, sim_props, output_dir, fig_format):
    """Plot the Rvir distribution for this snapshot"""

    type0_sel = gals["Type"]==0
    type1_sel = gals["Type"]==1

    fig, ax = plt.subplots(1,1)
    ax.set_xlabel(r"$\log_{10}(R_{\rm vir}/{\rm Mpc})$")
    ax.set_ylabel(r"$\phi\ [{\rm dex^{-1}\,Mpc^{-3}}]$")
    ax.set_yscale('log', nonposy='clip')
    dashes = (8,2)

    # Start with all halos (that hold galaxies)
    rvir = np.log10(gals["Rvir"][(type0_sel|type1_sel) & (gals["Rvir"]>0)])
    mf = munge.mass_function(rvir, sim_props["Volume"], bins=20)
    ax.plot(mf[:,0], mf[:,1], lw=4, label="All haloes hosting galaxies")

    # Type 0
    rvir = np.log10(gals["Rvir"][type0_sel & (gals["Rvir"]>0)])
    mf = munge.mass_function(rvir, sim_props["Volume"], bins=20)
    l, = ax.plot(mf[:,0], mf[:,1], ls='--', alpha=0.7, label="Type 0")
    l.set_dashes(dashes)

    # Type 1
    rvir = np.log10(gals["Rvir"][type1_sel & (gals["Rvir"]>0)])
    mf = munge.mass_function(rvir, sim_props["Volume"], bins=20)
    l, = ax.plot(mf[:,0], mf[:,1], ls='--', alpha=0.7, label="Type 1")
    l.set_dashes(dashes)

    # Finish up plot and save
    ax.text(0.97,0.95, "z=%.2f"%redshift,
           horizontalalignment='right',
           verticalalignment='top',
           transform=ax.transAxes)
    ax.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"Rvir_distrib."+fig_format))

if __name__ == '__main__':

    # deal with the command line arguments
    args        = docopt(__doc__)
    master_file = args['<master_file>']
    snapshot    = int(args['--snapshot'])
    output_dir  = args['--output']
    fig_format  = args['--format']

    # Read in the galaxies
    gals, sim_props = read_gals(master_file, snapshot=snapshot, sim_props=True)
    redshift = sim_props['Redshift']

    # Make the output directory if it doesn't already exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Generate all of the plots which we can for this redshift
    Rvir_distrib(gals, redshift, sim_props, output_dir, fig_format)
