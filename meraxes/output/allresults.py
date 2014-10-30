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

from dragons.meraxes.io import read_gals
from dragons import munge

__author__ = "Simon Mutch"
__date__   = "2013/08/06"

def smf_z0(gals, sim_props, output_dir, fig_format):
    """Plot the redshift 0 stellar mass function against observational data."""

    hubble_h = sim_props["Hubble_h"]

    # xlim = (9.4,12)
    # ylim = (1e-5,2e-2)
    xlim = (-25,13.5)
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
    colors = gals["MagDust"][sel,0]-gals["MagDust"][sel,1]
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
    line_colors = plt.rcParams["axes.color_cycle"]

    # Read in, munge and plot the observational data one by one
    # Baldry+ 2008
    obs = Table.read(os.path.join(__script_dir__,"../utils/obs_datasets/smf/Baldry2008-total.txt"), format='ascii')
    obs_mass  = obs['col1']-np.log10(0.7) # IMF correction to Salpeter
    obs_mass  = obs_mass+(2.*np.log10(0.7))-(2.*np.log10(hubble_h))  # hubble conversion
    hubble_cor = 1./(0.7**3)*(hubble_h**3) # hubble conversion
    obs_phi   = obs['col3']*hubble_cor
    obs_merr  = (obs['col3']-obs['col5'])*hubble_cor
    obs_perr  = (obs['col6']-obs['col3'])*hubble_cor
    l, _, _ = ax.errorbar(obs_mass, obs_phi, yerr=(obs_merr, obs_perr), ls="--", color='0.75', marker='s', label="Baldry 2008")
    l.set_dashes([8,2])

    # Bell+ 2003
    obs = Table.read(os.path.join(__script_dir__,"../utils/obs_datasets/smf/Bell2003-total.txt"), format='ascii')
    obs_mass = obs['col1']-np.log10(0.7)-2.*np.log10(hubble_h) # -log10(0.7) -> IMF correction to Salpeter
    obs_phi  = obs['col2']*(hubble_h**3)
    obs_merr = (obs['col2']-obs['col3'])*(hubble_h**3)
    obs_perr = (obs['col4']-obs['col2'])*(hubble_h**3)
    l, _, _ = ax.errorbar(obs_mass, obs_phi, yerr=(obs_merr, obs_perr), ls="--", color='0.5', marker='^', label="Bell 2003")
    l.set_dashes([8,2])

    # Plot the model
    ax.plot(mf[:,0], mf[:,1], lw=4, color='k', label="Meraxes total")
    ax.plot(mf_red[:,0], mf_red[:,1], lw=4, color=line_colors[1], label="red")
    ax.plot(mf_blue[:,0], mf_blue[:,1], lw=4, color=line_colors[0], label="blue")

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



def bJ_LF_z0(gals, sim_props, output_dir, fig_format):
    """Plot the redshift 0 bJ band luminosity function against observational data."""

    hubble_h = sim_props["Hubble_h"]

    # xlim = (9.4,12)
    # ylim = (1e-5,2e-2)
    xlim = (-24,-13)
    ylim = (1e-8,1e-1)
    dashes = [8,2]

    # Generate the model mass function
    mag = gals["MagDust"][:, 0] - 0.28 * (gals["MagDust"][:, 0] - gals["MagDust"][:,1])
    mf = munge.mass_function(mag, sim_props["Volume"], range=(xlim[0]-0.1,xlim[1]+0.1), bins='knuth')

    # Setup the plot
    fig, ax = plt.subplots(1,1)
    ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(r"$M_{\rm b_J}$")
    ax.set_ylabel(r"$\phi\ [{\rm dex^{-1}\,Mpc^{-3}}]$")
    ax.set_xlim(np.array(xlim)[::-1])
    ax.set_ylim(ylim)
    line_colors = plt.rcParams["axes.color_cycle"]

    # Data of Norberg et al. (2002) bJ band 2dFGRS LF
    obs_phi = np.array([0.0, 4.9295418E-02, 4.8414569E-02, 6.1575811E-02, 4.0561352E-02,
               4.4041574E-02, 3.4943096E-02, 3.4614738E-02, 3.2133184E-02,
               2.9635418E-02, 2.3947729E-02, 2.4405226E-02, 2.3334850E-02,
               2.1415474E-02, 1.9061793E-02, 1.8999882E-02, 1.6375385E-02,
               1.5369940E-02, 1.4364687E-02, 1.3468309E-02, 1.2573769E-02,
               1.1145771E-02, 9.4299037E-03, 8.0003440E-03, 6.1532548E-03,
               4.5839567E-03, 3.0883602E-03, 1.9342416E-03, 1.0532817E-03,
               5.1049568E-04, 2.1174403E-04, 7.1619586E-05, 2.1268283E-05,
               5.9436675E-06, 1.8153961E-06, 5.1464536E-07, 3.0028801E-07,])
    obs_phi *= (hubble_h**3.)

    obs_err_phi = np.array([1.0000000E-10, 2.4690846E-02, 1.2388939E-02, 9.7498214E-03,
                   5.7149841E-03, 4.5108106E-03, 2.9827626E-03, 2.3375382E-03,
                   1.8226152E-03, 1.4552969E-03, 1.0482093E-03, 8.9600182E-04,
                   7.4631697E-04, 6.2729424E-04, 4.9186783E-04, 4.0296704E-04,
                   2.9608142E-04, 2.3880773E-04, 1.9196754E-04, 1.5552285E-04,
                   1.2866473E-04, 9.9682424E-05, 7.6998745E-05, 5.9886646E-05,
                   4.4532266E-05, 3.3105953E-05, 2.3582392E-05, 1.6317146E-05,
                   1.0632872E-05, 6.6587650E-06, 4.0021373E-06, 2.2571021E-06,
                   1.2226549E-06, 6.4833682E-07, 3.6088997E-07, 1.9372916E-07,
                   1.4984870E-07,])
    obs_err_phi *= (hubble_h**3.)

    obs_mag = np.array([-13.00000, -13.27500, -13.55000, -13.82500, -14.10000,
               -14.37500, -14.65000, -14.92500, -15.20000, -15.47500,
               -15.75000, -16.02500, -16.30000, -16.57500, -16.85000,
               -17.12500, -17.40000, -17.67500, -17.95000, -18.22500,
               -18.50000, -18.77500, -19.05000, -19.32500, -19.60000,
               -19.87500, -20.15000, -20.42500, -20.699999999999999,
               -20.97500, -21.25000, -21.52500, -21.80000, -22.07500,
               -22.35000, -22.62500, -22.90000,])
    obs_mag += 5.*np.log10(hubble_h)

    # Plot the observations
    plt.errorbar(obs_mag, obs_phi, yerr=obs_err_phi, color=line_colors[1],
                 marker='s', ms=7, ls='none',
                 label='Norberg et al. (2002)')

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
    figname = os.path.join(output_dir, "bJ_LF_z0.")+fig_format
    print "Saving figure: %s" % figname
    plt.savefig(figname)


def K_LF_z0(gals, sim_props, output_dir, fig_format):
    """Plot the redshift 0 K-band luminosity function against observational data."""

    hubble_h = sim_props["Hubble_h"]

    xlim = (-28,-18)
    ylim = (1e-8,1e-1)
    dashes = [8,2]

    # Generate the model mass function
    mag = gals["MagDust"][:, 4]
    print mag[(mag>(xlim[0]-0.1)) & (mag<(xlim[1]+0.1))].shape
    mf = munge.mass_function(mag, sim_props["Volume"], range=(xlim[0]-0.1,xlim[1]+0.1), bins='knuth')

    # Setup the plot
    fig, ax = plt.subplots(1,1)
    ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(r"$M_{\rm K}$")
    ax.set_ylabel(r"$\phi\ [{\rm dex^{-1}\,Mpc^{-3}}]$")
    ax.set_xlim(np.array(xlim)[::-1])
    ax.set_ylim(ylim)
    line_colors = plt.rcParams["axes.color_cycle"]

    # Cole et al. K band 2dFGRS LF
    obs_phi = np.array([3.1315561E-03, 8.2625253E-03, 0.0, 4.6483092E-03,
                        5.7576019E-03, 9.1649834E-03, 1.1232893E-02,
                        1.0536440E-02, 8.5763102E-03, 8.8181989E-03,
                        6.9448259E-03, 6.0896124E-03, 9.2596142E-03,
                        6.9631678E-03, 7.2867479E-03, 6.9923755E-03,
                        5.9844730E-03, 5.9305103E-03, 5.3865365E-03,
                        5.8525647E-03, 5.2373926E-03, 4.9635037E-03,
                        4.1801766E-03, 2.7171015E-03, 1.8800517E-03,
                        1.2136410E-03, 6.5419916E-04, 3.4594961E-04,
                        1.4771589E-04, 5.5521199E-05, 2.1283222E-05,
                        9.4211919E-06, 1.0871951E-06, 2.7923562E-07,])
    obs_phi = obs_phi * (hubble_h**3.)

    obs_err_phi = np.array([ 3.6377162E-03, 6.6833422E-03, 1.0000000E-10,
                            4.0996978E-03, 4.3155681E-03, 5.6722397E-03,
                            6.4211683E-03, 5.7120644E-03, 4.6346937E-03,
                            3.8633577E-03, 2.4383855E-03, 1.6279118E-03,
                            1.6941463E-03, 1.1781409E-03, 9.7785855E-04,
                            7.9027453E-04, 6.0649612E-04, 5.1598746E-04,
                            4.2267537E-04, 3.7395130E-04, 2.8177485E-04,
                            2.1805518E-04, 1.6829016E-04, 1.1366483E-04,
                            8.1871600E-05, 5.7472309E-05, 3.6554517E-05,
                            2.3141622E-05, 1.2801432E-05, 6.5092854E-06,
                            3.3540452E-06, 1.9559407E-06, 5.6035748E-07,
                            2.8150106E-07, ])
    obs_err_phi = obs_err_phi * (hubble_h**3.)

    obs_mag = np.array([ -18.00000, -18.25000, -18.50000, -18.75000,
                          -19.00000, -19.25000, -19.50000, -19.75000,
                          -20.00000, -20.25000, -20.50000, -20.75000, -21.0,
                          -21.25000, -21.50000, -21.75000, -22.00000,
                          -22.25000, -22.5, -22.75000, -23.00000, -23.25000,
                          -23.50000, -23.75000, -24.0, -24.25000, -24.50000,
                          -24.75000, -25.00000, -25.25000, -25.5, -25.75000,
                          -26.00000, -26.25000, ])
    obs_mag += 5.*np.log10(hubble_h)

    # Plot the observations
    plt.errorbar(obs_mag, obs_phi, yerr=obs_err_phi, color=line_colors[1],
                 marker='s', ms=7, ls='none',
                 label='Cole et al. (2001)')

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
    figname = os.path.join(output_dir, "K_LF_z0.")+fig_format
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
        # smf_z0(gals, sim_props, output_dir, fig_format)
        smf_colorsplit_z0(gals, sim_props, output_dir, fig_format)
        bJ_LF_z0(gals, sim_props, output_dir, fig_format)
        K_LF_z0(gals, sim_props, output_dir, fig_format)

