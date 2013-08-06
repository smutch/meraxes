#!/usr/bin/env python

"""Various mass history plots."""

import numpy as np
import matplotlib.pyplot as plt

__author__ = "Simon Mutch"
__date__   = "2013/06/25"

def mvir(gals, redshift, sim_props, ax, ylim=(8,15), xlim=None, bins=20, **kwargs):
    
    # Set up the axis
    ax.set_xlabel(r"$z+1$")
    ax.set_ylabel(r"$\log_{10}(h_{%.2f}\,M_{vir}/M_{\odot})$" %
                  sim_props["Hubble_h"])
    if xlim!=None:
        ax.set_xlim(xlim)
    if ylim!=None:
        ax.set_ylim(ylim)
    # ax.xaxis.set_minor_locator(plt.MaxNLocator(20))
    ax.set_xscale('log') 

    # Plot
    lines = []
    for g in gals:
        sel = [(g['Type']>-1) & (g['Type']<2)]
        z = redshift[sel]
        m = np.log10(g[sel]['Mvir']*1.e10)
        line, = ax.plot(z+1, m, label=g[sel][0]['ID'], **kwargs)
        lines.append(line)

    return lines

