#!/usr/bin/env python

"""Usage: gen_weighted_Mvir_crit.py <fname>"""

import numpy as np
import matplotlib.pyplot as plt
from dragons import meraxes, plotutils
from docopt import docopt
from tqdm import tqdm
import h5py as h5

__author__ = "Simon Mutch"
__date__ = "2015-06-22"

plt.style.use(['dragons', 'paper'])


def gen_weighted_Mvir_crit(fname):
    props = ("FOFMvir", "Type", "GhostFlag", "MvirCrit")

    snaplist, zlist, _ = meraxes.read_snaplist(fname)
    mvir_crit = np.zeros(snaplist.size)

    for ii, snap in tqdm(enumerate(snaplist), total=snaplist.size):
        try:
            gals = meraxes.read_gals(fname, snapshot=snap, props=props,
                                     quiet=True).view(np.recarray)
            gals = gals[(gals.Type == 0) & (gals.GhostFlag == 0)]
            mvir_crit[ii] = gals.MvirCrit.mean()
        except IndexError:
            pass

    return mvir_crit, zlist, snaplist


def plot_mvir_crit(mvir_crit, zlist, snapslist):
    fig, ax = plt.subplots(1, 1)
    ax.plot(zlist, mvir_crit)
    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$M_{\rm filt}$")
    ax.set_xlim(5, 15)
    plt.savefig('mean_Mvir_crit_gal.png')


def save(mvir_crit, zlist, snaplist):
    fout = h5.File("mean_Mvir_crit_gal.hdf5", "w")
    fout.create_dataset("snapshot", data=snaplist)
    fout.create_dataset("redshift", data=zlist)
    fout.create_dataset("mean_Mvir_crit", data=mvir_crit)
    fout.close()


if __name__ == '__main__':
    args = docopt(__doc__)
    fname = args["<fname>"]
    mvir_crit, zlist, snaplist = gen_weighted_Mvir_crit(fname)
    plot_mvir_crit(mvir_crit, zlist, snaplist)
    save(mvir_crit, zlist, snaplist)
