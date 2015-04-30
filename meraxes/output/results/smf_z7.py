#!/usr/bin/env python

"""Plot the z=7 SMF on it's own.

Usage: smf_z7.py <fname>

"""

import numpy as np
import matplotlib.pyplot as plt
from dragons import meraxes, plotutils
from docopt import docopt
from smf import Smf, schechter
import pandas as pd

import os
import sys
__script_dir__ = os.path.dirname(os.path.realpath(__file__))
sys.path.append(__script_dir__)

__author__ = "Simon Mutch"
__date__ = "2015-04-14"

COLS = plt.rcParams["axes.color_cycle"]
OBS_DATASETS_DIR = "/home/smutch/models/21cm_sam/meraxes/utils/obs_datasets"
XLIM = [6, 12]
YLIM = [-6, 1]
HUBBLE = 0.7


def plot_obs(ax):
    # Start with Duncan+ 2014
    obs = pd.read_table(os.path.join(__script_dir__,
        "{:s}/smf/Duncan14_MF_z7.txt".format(OBS_DATASETS_DIR)),
        delim_whitespace=True,
        header=None,
        skiprows=8,
        names=["sm", "phi", "merr", "perr"])
    obs.merr[obs.merr >= obs.phi] = obs.phi - 1e-10

    # convert obs to same hubble value and IMF
    obs.sm += 2.0*np.log10(0.702/HUBBLE)
    obs.sm += 0.25  # IMF correction Chabrier -> Salpeter
    for col in ["phi", "merr", "perr"]:
        obs[col] /= (0.7**3/HUBBLE**3)

    # plot the observations
    ax.errorbar(obs.sm, np.log10(obs.phi),
                yerr=[np.log10(obs.phi / (obs.phi - obs.merr)),
                      np.log10(1.0 + (obs.perr / obs.phi))],
                label="Duncan et al. 2014", ls="None", color=COLS[1],
                lw=2, capsize=2.5, marker='o', mec='None')

    # now plot Grazian+ 2015 (h=0.7, Salpeter IMF)
    mass = np.linspace(9.76 + 2.0*np.log10(0.7/HUBBLE), XLIM[1], 30)
    phi = np.log10(schechter(10.69, -5.24, -1.88,
                             mass + 2.0*np.log10(HUBBLE/0.7)))\
            + 3.0*np.log10(HUBBLE/0.7)
    ax.plot(mass, phi, ls="--", label="Grazian et al. 2015", color=COLS[3],
            lw=3)

    # now plot Gonzalez+ 2011 (h=0.7, Chabrier IMF)
    obs = pd.read_table(os.path.join(__script_dir__,
        "{:s}/smf/Gonzalez11_z7_smf.txt".format(OBS_DATASETS_DIR)),
        delim_whitespace=True,
        header = None,
        skiprows = 3,
        names = ["sm", "log_phi", "m_err", "p_err"])

    # convert obs to same hubble value and IMF
    obs.sm += 2.0*np.log10(0.7/HUBBLE)
    obs.sm += 0.25  # IMF correction Chabrier -> Salpeter
    for col in ["log_phi", "m_err", "p_err"]:
        obs[col] -= 3.0*np.log10(0.7/HUBBLE)

    # plot the observations
    ax.errorbar(obs.sm, obs.log_phi, yerr=[obs.m_err, obs.p_err],
                label="Gonzalez et al. 2011", ls="None", color=COLS[2],
                lw=2, capsize=2.5, marker='s', mec='None')


def plot(model, ax):
    model.generate(limits=XLIM)
    # model.gen_uncert(20)
    model.plot(ax, label="Meraxes", lw=4)
    plot_obs(ax)

    model.set_axlabel(ax)
    ax.axis(XLIM + YLIM)
    ax.legend(loc=1)

    # add some text
    ax.text(0.05, 0.05, "z=7\nh=%.1f" % HUBBLE,
            horizontalalignment="left", verticalalignment="bottom",
            transform=ax.transAxes)


if __name__ == '__main__':
    args = docopt(__doc__)
    fname = args["<fname>"]

    # set up
    meraxes.set_little_h(HUBBLE)
    plt.style.use(["dragons", ])
    fig, ax = plt.subplots(1, 1)

    model = Smf(fname, 7.0)
    plot(model, ax)

    # save
    plt.tight_layout()
    plt.savefig(os.path.dirname(fname)+"/plots/smf-z7.pdf")
