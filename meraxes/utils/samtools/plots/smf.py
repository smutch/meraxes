#!/usr/bin/env python

"""Plot the stellar mass function"""

import numpy as np
import matplotlib.pyplot as plt

__author__ = "Simon Mutch"
__date__   = "2013/06/19"

def smf(gals, sim_props, ax, xlim=(8.5,12.5), ylim=None, bins=20, ls='-', color=None):

    # Set up the axis
    ax.set_ylabel(r"$h^3_{%.2f}\,\phi\ [Mpc^{-3}\,dex^{-1}]$" %
                  sim_props["Hubble_h"])
    ax.set_xlabel(r"$\log_{10}(h_{%.2f}\,M_{*}/M_{\odot})$" %
                  sim_props["Hubble_h"])
    ax.set_yscale('log', nonposy='clip')
    ax.set_xlim(xlim)

    # Calculate the mass function
    sel = (gals['StellarMass']>0)
    stellar = np.log10(gals[sel]['StellarMass']*1.e10)
    stellar = stellar[(stellar>=xlim[0]) & (stellar<xlim[1])]
    n, edges = np.histogram(stellar, bins=bins, normed=False) 
    bin_width = edges[1]-edges[0]
    phi = n/sim_props["Volume"]/bin_width
    bin_centers = edges[:-1]+0.5*bin_width

    # Plot
    line, = ax.plot(bin_centers, phi, ls=ls, lw=3, label="Meraxes")
    if color!=None:
        line.set_color(color)

    return line
