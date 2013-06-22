#!/usr/bin/env python

"""Plots specific to the formation history model"""

import numpy as np
import matplotlib.pyplot as plt

__author__ = "Simon Mutch"
__date__   = "2013/06/22"

def fhmodel_dmvirdt(gals, sim_props, ax, xlim=(-4,8), ylim=None, bins=40, ls='-', color=None):

    # Set up the axis
    ax.set_ylabel(r"$h^3_{%.2f}\,\phi\ [Mpc^{-3}\,dex^{-1}]$" %
                  sim_props["Hubble_h"])
    ax.set_xlabel(r"$\log_{10}(\frac{dM_{\rm vir}}{dt}/(M_{\odot}/yr))$")
    ax.set_yscale('log', nonposy='clip')
    ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
    if xlim!=None:
        ax.set_xlim(xlim)

    # Calculate the distribution
    sel = (gals['dMdt']>0)
    dmdt = np.log10(gals[sel]['dMdt'])
    if xlim!=None:
        dmdt = dmdt[(dmdt>=xlim[0]) & (dmdt<xlim[1])]
    n, edges = np.histogram(dmdt, bins=bins, normed=False) 
    bin_width = edges[1]-edges[0]
    phi = n/sim_props["Volume"]/bin_width
    bin_centers = edges[:-1]+0.5*bin_width

    # Plot
    line, = ax.plot(bin_centers, phi, ls=ls, lw=3, label="Meraxes")
    if color!=None:
        line.set_color(color)

    return line
