#!/usr/bin/env python

"""Randomly sample a forests_info file to generate a list of forest_ids.

Usage: sample_forests_ids.py <sample_frac> <seed> <forests_info_file> <box_size>
"""

import numpy as np
import h5py as h5
from docopt import docopt
from random import sample, seed
import matplotlib.pylab as plt
import pandas as pd
from dragons import munge, nbody
import os

__author__ = "Simon Mutch"
__date__ = "2014-05-04"

args = docopt(__doc__)


def sample_forest_ids(sample_frac, forests_info_file):

    with h5.File(forests_info_file, "r") as fd:
        forest_ids = fd["info"]["forest_id"][:]

    n_to_sample = int(forest_ids.shape[0] * sample_frac)

    print "Sampling {:d} of {:d} forest_ids...".format(n_to_sample,
                                                       forest_ids.shape[0])

    sampled = forest_ids[sample(xrange(forest_ids.shape[0]), n_to_sample)]
    sampled.sort()

    output_file = "sampled_forest_ids-{:.2f}.txt".format(sample_frac)
    with open(output_file, "w") as fd:
        fd.write("{:d}\n".format(sampled.shape[0]))
        np.savetxt(fd, sampled, fmt='%d')

    return sampled


def plot_hmf(sampled_ids, sample_frac, forests_info_file, volume):

    print "Plotting sampled & unsampled mass functions"

    # read in the trees
    tree_fname = os.path.join(os.path.dirname(forests_info_file),
                              "horizontal_trees_099.hdf5")
    with h5.File(tree_fname, "r") as fd:
        trees = pd.DataFrame(fd["trees"][:])

    # read in the catalogues
    cat_dir = os.path.join(os.path.dirname(forests_info_file),
                           "../catalogs")
    mvir = {}
    for cat_type in ["groups", "subgroups"]:
        mvir[cat_type] = nbody.read_halo_catalog(
            cat_dir + "/subfind_099.catalog_{:s}_properties"
            .format(cat_type))["M_vir"]

    # select only type 0 halos
    # halos = halos[halos["central_index"] == halos.index]

    # plot the mass functions
    fig, ax = plt.subplots(1, 1)
    mf, edges = munge.mass_function(np.log10(mvir["groups"]), volume, 50,
                                    return_edges=True)
    l, = ax.plot(mf[:, 0], np.log10(mf[:, 1]), label="groups")
    mf = munge.mass_function(np.log10(mvir["subgroups"]), volume, bins=edges)
    ax.plot(mf[:, 0], np.log10(mf[:, 1]),
            color=l.get_color(), ls="--", label="subgroups")

    # plot the sampled mass functions
    sampled_ids = pd.DataFrame(sampled_ids, columns=["forest_id", ])
    trees["Mvir"] = mvir["subgroups"]
    trees = trees.merge(sampled_ids, on="forest_id")

    group_ind = np.unique(trees.group_index.values)
    mf = munge.mass_function(
        np.log10(mvir["groups"][group_ind]), sample_frac * volume, bins=edges)
    l, = ax.plot(mf[:, 0], np.log10(mf[:, 1]), label="sampled")
    mf = munge.mass_function(
        np.log10(trees.Mvir.values), sample_frac * volume, bins=edges)
    ax.plot(mf[:, 0], np.log10(mf[:, 1]), ls="--", color=l.get_color())

    ax.legend(loc=1)

    plt.savefig("sampled_mf.png")


if __name__ == '__main__':
    sample_frac = float(args["<sample_frac>"])
    forests_info_file = args["<forests_info_file>"]
    volume = float(args["<box_size>"])**3
    rseed = int(args["<seed>"])

    seed(rseed)
    sampled_ids = sample_forest_ids(sample_frac, forests_info_file)
    plot_hmf(sampled_ids, sample_frac, forests_info_file, volume)
