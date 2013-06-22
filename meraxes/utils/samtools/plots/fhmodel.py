#!/usr/bin/env python

"""Plots specific to the formation history model"""

import numpy as np
import matplotlib.pyplot as plt

__author__ = "Simon Mutch"
__date__   = "2013/06/22"

def dmvirdt(gals, sim_props, ax, xlim=(-4,8), ylim=None, bins=40, ls='-', color=None):

    # Set up the axis
    ax.set_ylabel(r"$h^3_{%.2f}\,\phi\ [Mpc^{-3}\,dex^{-1}]$" %
                  sim_props["Hubble_h"])
    ax.set_xlabel(r"$\log_{10}(\frac{dM_{\rm vir}}{dt}/(M_{\odot}/yr))$")
    ax.set_yscale('log', nonposy='clip')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    if xlim!=None:
        ax.set_xlim(xlim)
    if ylim!=None:
        ax.set_ylim(ylim)

    # Calculate the distribution
    sel = (gals['dMdt']>0) & (gals['Type']==0)
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


def dmvirdt_vs_mvir(gals, sim_props, ax, xlim=(8,15), ylim=(1,1e7), bins=30, ls='-', color=None):

    # Set up the axis
    ax.set_ylabel(r"mean $\frac{dM_{vir}}{dt}\, [M_{\odot}/yr]$")
    ax.set_xlabel(r"$h^{-1}_{%.2f}\,M_{vir}\,[M_{\odot}]$" %
                  sim_props["Hubble_h"])
    ax.set_yscale('log', nonposy='clip')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    if xlim!=None:
        ax.set_xlim(xlim)
    if ylim!=None:
        ax.set_ylim(ylim)

    # Calculate the distribution
    sel = (gals['dMdt']>0) & (gals['Type']==0) & (gals['Mvir']>0)
    if xlim!=None:
        mvir = np.log10(gals['Mvir']*1.e10)
        sel = sel & (mvir>=xlim[0]) & (mvir<xlim[1])
    mvir = np.log10(gals[sel]['Mvir']*1.e10)
    dmdt = gals[sel]['dMdt']
    n, edges = np.histogram(mvir, bins=bins, normed=False) 
    edges[-1]=edges[-1]+1e-6
    ind = np.digitize(mvir, edges)-1
    mean_dmdt = np.array([dmdt[ind==ii].mean() for ii in range(bins)])
    std_dmdt = np.array([dmdt[ind==ii].std() for ii in range(bins)])
    bin_width = edges[1]-edges[0]
    bin_centers = edges[:-1]+0.5*bin_width

    # Plot
    line, = ax.plot(bin_centers, mean_dmdt, ls=ls, lw=3, label="Meraxes")
    poly  = ax.fill_between(bin_centers, mean_dmdt+std_dmdt,
                             mean_dmdt-std_dmdt, facecolor=line.get_color(),
                             edgecolor='none',
                             alpha=0.5)
    
    if color!=None:
        line.set_color(color)
        for p in poly:
            p.set_color(color)

    return line, poly
