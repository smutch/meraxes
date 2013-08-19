"""Plot all results...

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

from samtools.io import read_gals
from samtools import munge

__author__ = "Simon Mutch"
__date__   = "2013/08/06"

def smf_z0(gals, sim_props, output_dir, fig_format):
    """Plot the redshift 0 stellar mass function against observational data."""

    hubble_h = sim_props["Hubble_h"]

    # xlim = (9.4,12)
    # ylim = (1e-5,2e-2)
    xlim = (8,12)
    ylim = (1e-6,1e-1)
    dashes = [8,2]

    # Generate the model mass function
    stars = np.log10(gals["StellarMass"][gals["StellarMass"]>0]*1.e10)
    mf = munge.mass_function(stars, sim_props["Volume"], range=(xlim[0]-0.1,xlim[1]+0.1), bins='knuth')

    # Setup the plot
    fig, ax = plt.subplots(1,1)
    ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(r"$\log_{10}(M_*/{\rm M_{\odot}})$")
    ax.set_ylabel(r"$\phi\ [{\rm dex^{-1}\,Mpc^{-3}}]$")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Read in, munge and plot the observational data one by one
    # Baldry+ 2008
    obs = Table.read(os.path.join(__script_dir__,"../utils/obs_datasets/smf/Baldry2008-total.txt"), format='ascii')
    obs_mass  = obs['col1']-np.log10(0.7) # IMF correction to Salpeter
    obs_mass  = obs_mass+(2.*np.log10(0.7))-(2.*np.log10(hubble_h))  # hubble conversion
    hubble_cor = 1./(0.7**3)*(hubble_h**3) # hubble conversion
    obs_phi   = obs['col3']*hubble_cor
    obs_merr  = (obs['col3']-obs['col5'])*hubble_cor
    obs_perr  = (obs['col6']-obs['col3'])*hubble_cor
    l, _, _ = ax.errorbar(obs_mass, obs_phi, yerr=(obs_merr, obs_perr), ls="--", marker='s', label="Baldry 2008")
    l.set_dashes([8,2])

    # Bell+ 2003
    obs = Table.read(os.path.join(__script_dir__,"../utils/obs_datasets/smf/Bell2003-total.txt"), format='ascii')
    obs_mass = obs['col1']-np.log10(0.7)-2.*np.log10(hubble_h) # -log10(0.7) -> IMF correction to Salpeter
    obs_phi  = obs['col2']*(hubble_h**3)
    obs_merr = (obs['col2']-obs['col3'])*(hubble_h**3)
    obs_perr = (obs['col4']-obs['col2'])*(hubble_h**3)
    l, _, _ = ax.errorbar(obs_mass, obs_phi, yerr=(obs_merr, obs_perr), ls="--", marker='s', label="Bell 2003")
    l.set_dashes([8,2])

    # Plot the model
    ax.plot(mf[:,0], mf[:,1], lw=4, label="Meraxes")

    # Last few plotting bits and pieces
    ax.legend(loc="lower left", numpoints=1)
    ax.text(0.96,0.95, 
            "z=0\nAssumed Salpeter IMF\nh={:.3f}".format(hubble_h),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    plt.tight_layout()

    # Save the plot
    figname = os.path.join(output_dir, "smf_z0.")+fig_format
    print "Saving figure: %s" % figname
    plt.savefig(figname)



def smf_colorsplit_z0(gals, sim_props, output_dir, fig_format):
    """Plot the redshift 0 color split stellar mass function against observational data."""

    hubble_h = sim_props["Hubble_h"]

    # xlim = (9.4,12)
    # ylim = (1e-5,2e-2)
    xlim = (8,12)
    ylim = (1e-6,1e-1)
    dashes = [8,2]

    # Generate the model mass function
    sel = (gals["StellarMass"]>0)
    stars = np.log10(gals["StellarMass"][sel]*1.e10)
    colors = gals["Mag"][sel,0]-gals["Mag"][sel,1]
    red_sel = (colors>=0.8)
    blue_sel = (colors<0.8)
    mf, edges = munge.mass_function(stars, sim_props["Volume"], range=(xlim[0]-0.1,xlim[1]+0.1), bins='knuth', return_edges=True)
    mf_red = munge.mass_function(stars[red_sel], sim_props["Volume"], range=(xlim[0]-0.1,xlim[1]+0.1), bins=edges)
    mf_blue = munge.mass_function(stars[blue_sel], sim_props["Volume"], range=(xlim[0]-0.1,xlim[1]+0.1), bins=edges)

    # import IPython; IPython.embed()

    # Setup the plot
    fig, ax = plt.subplots(1,1)
    ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(r"$\log_{10}(M_*/{\rm M_{\odot}})$")
    ax.set_ylabel(r"$\phi\ [{\rm dex^{-1}\,Mpc^{-3}}]$")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Read in, munge and plot the observational data one by one
    # Baldry+ 2008
    obs = Table.read(os.path.join(__script_dir__,"../utils/obs_datasets/smf/Baldry2008-total.txt"), format='ascii')
    obs_mass  = obs['col1']-np.log10(0.7) # IMF correction to Salpeter
    obs_mass  = obs_mass+(2.*np.log10(0.7))-(2.*np.log10(hubble_h))  # hubble conversion
    hubble_cor = 1./(0.7**3)*(hubble_h**3) # hubble conversion
    obs_phi   = obs['col3']*hubble_cor
    obs_merr  = (obs['col3']-obs['col5'])*hubble_cor
    obs_perr  = (obs['col6']-obs['col3'])*hubble_cor
    l, _, _ = ax.errorbar(obs_mass, obs_phi, yerr=(obs_merr, obs_perr), ls="--", marker='s', label="Baldry 2008")
    l.set_dashes([8,2])

    # Bell+ 2003
    obs = Table.read(os.path.join(__script_dir__,"../utils/obs_datasets/smf/Bell2003-total.txt"), format='ascii')
    obs_mass = obs['col1']-np.log10(0.7)-2.*np.log10(hubble_h) # -log10(0.7) -> IMF correction to Salpeter
    obs_phi  = obs['col2']*(hubble_h**3)
    obs_merr = (obs['col2']-obs['col3'])*(hubble_h**3)
    obs_perr = (obs['col4']-obs['col2'])*(hubble_h**3)
    l, _, _ = ax.errorbar(obs_mass, obs_phi, yerr=(obs_merr, obs_perr), ls="--", marker='s', label="Bell 2003")
    l.set_dashes([8,2])

    # Plot the model
    ax.plot(mf[:,0], mf[:,1], lw=4, label="Meraxes total")
    ax.plot(mf_red[:,0], mf_red[:,1], lw=4, label="red")
    ax.plot(mf_blue[:,0], mf_blue[:,1], lw=4, label="blue")

    # Last few plotting bits and pieces
    ax.legend(loc="lower left", numpoints=1)
    ax.text(0.96,0.95, 
            "z=0\nAssumed Salpeter IMF\nh={:.3f}".format(hubble_h),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    plt.tight_layout()

    # Save the plot
    figname = os.path.join(output_dir, "csplit_smf_z0.")+fig_format
    print "Saving figure: %s" % figname
    plt.savefig(figname)


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
    if redshift<=0.1:
        smf_z0(gals, sim_props, output_dir, fig_format)
        smf_colorsplit_z0(gals, sim_props, output_dir, fig_format)

