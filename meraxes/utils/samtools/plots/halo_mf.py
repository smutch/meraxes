#!/usr/bin/env python

"""Plot the halo mass function"""

import numpy as np
import matplotlib.pyplot as plt

__author__ = "Simon Mutch"
__date__   = "2013/06/21"

def halo_mf(gals, sim_props, ax, xlim=(9,14), ylim=None, bins=20, ls='-', color=None):

    # Set up the axis
    ax.set_ylabel(r"$h^3_{%.2f}\,\phi\ [Mpc^{-3}\,dex^{-1}]$" %
                  sim_props["Hubble_h"])
    ax.set_xlabel(r"$\log_{10}(h_{%.2f}\,M_{\rm vir}/M_{\odot})$" %
                  sim_props["Hubble_h"])
    ax.set_yscale('log', nonposy='clip')
    ax.set_xlim(xlim)
    ax.xaxis.set_minor_locator(plt.MaxNLocator(20))

    # Calculate the mass function
    sel = (gals['Type']==0)
    mvir = np.log10(gals[sel]['Mvir']*1.e10)
    mvir = mvir[(mvir>=xlim[0]) & (mvir<xlim[1])]
    n, edges = np.histogram(mvir, bins=bins, normed=False) 
    bin_width = edges[1]-edges[0]
    phi = n/sim_props["Volume"]/bin_width
    bin_centers = edges[:-1]+0.5*bin_width

    # Plot the result
    line, = ax.plot(bin_centers, phi, ls=ls, lw=3, label="Meraxes")
    if color!=None:
        line.set_color(color)

    return line

