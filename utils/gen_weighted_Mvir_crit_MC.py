#!/usr/bin/env python

"""Usage: gen_weighted_Mvir_crit_MC.py <fname>"""

import numpy as np
import matplotlib.pyplot as plt
from dragons import meraxes, plotutils
from docopt import docopt
from tqdm import tqdm
import h5py as h5

__author__ = "Emanuele Ventura"
__date__ = "2022-06-06"

plt.style.use(['dragons', 'paper'])


def gen_weighted_Mvir_crit_MC(fname):
    props = ("FOFMvir", "Type", "GhostFlag", "MvirCrit_MC")

    snaplist, zlist, _ = meraxes.read_snaplist(fname)
    mvir_crit_MC = np.zeros(snaplist.size)

    for ii, snap in tqdm(enumerate(snaplist), total=snaplist.size):
        try:
            gals = meraxes.read_gals(fname, snapshot=snap, props=props,
                                     quiet=True).view(np.recarray)
            gals = gals[(gals.Type == 0) & (gals.GhostFlag == 0)]
            mvir_crit_MC[ii] = gals.MvirCrit_MC.mean()
        except IndexError:
            pass

    return mvir_crit_MC, zlist, snaplist


def plot_mvir_crit_MC(mvir_crit_MC, zlist, snapslist):
    fig, ax = plt.subplots(1, 1)
    ax.plot(zlist, mvir_crit_MC)
    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$M_{\rm filt}$")
    ax.set_xlim(5, 15)
    plt.savefig('mean_Mvir_crit_MC_gal.png')


def save(mvir_crit_MC, zlist, snaplist):
    fout = h5.File("mean_Mvir_crit_MC_gal.hdf5", "w")
    fout.create_dataset("snapshot", data=snaplist)
    fout.create_dataset("redshift", data=zlist)
    fout.create_dataset("mean_Mvir_crit_MC", data=mvir_crit_MC)
    fout.close()


if __name__ == '__main__':
    args = docopt(__doc__)
    fname = args["<fname>"]
    mvir_crit_MC, zlist, snaplist = gen_weighted_Mvir_crit_MC(fname)
    plot_mvir_crit_MC(mvir_crit_MC, zlist, snaplist)
    save(mvir_crit_MC, zlist, snaplist)
