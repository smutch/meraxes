#!/usr/bin/env python

"""Subsample merger trees to maintain the halo mass function at a particular
redshift.

Usage: sample_trees.py <frac> <tree_dir> <snapshot> <prefix>
"""

import os
from ssimpl import munge
import numpy as np
import pandas as pd
import h5py as h5
import matplotlib.pyplot as plt
from docopt import docopt

__author__ = "Simon Mutch"
__date__   = "2014/05/02"

args = docopt(__doc__)


def sample_trees(frac, tree_dir, snapshot, prefix):

    # read in halos from tree file
    tree_file = os.path.join(tree_dir, prefix+"_{:03d}.hdf5".format(snapshot))
    with h5.File(tree_file, "r") as fd:
        halos = fd["trees"][:]
        halos = halos[["central_index", "forest_id", "fof_mvir"]]

    # select centrals only
    halos = halos[halos["central_index"] == np.arange(halos.shape[0])]
    halos = pd.DataFrame(halos)

    # bin halos by fof_mvir
    nbins = 100
    edges = np.linspace(halos.fof_mvir.min(), halos.fof_mvir.max()+0.0001, nbins+1)
    bind = np.digitize(halos.fof_mvir, edges)-1
    halos["bind"] = bind

    # group by bin index
    grouped = halos.groupby("bind")

    # loop through each bin and randomly drop requested fraction of
    # galaxies
    for b, g in grouped:
        n_to_drop = int(np.rint(float(g.shape[0]) * 1.0-frac))
        halos = halos.drop(np.random.choice(g.index, n_to_drop))

    # now write out the unique list of remaining forest ids
    with h5.File("sampled_forest_ids.hdf5", "w") as fd:
        fd.create_dataset("forest_id", data=halos["forest_id"].unique())



if __name__ == '__main__':
    tree_dir = args["<tree_dir>"]
    snapshot = int(args["<snapshot>"])
    prefix = args["<prefix>"]
    frac = float(args["<frac>"])
    sample_trees(frac, tree_dir, snapshot, prefix)
