#!/usr/bin/env python

"""Track the location of the most dense grid cell from a simulation at each
snapshot.

    Usage: track_density_cell.py <sim> <first_snap> <last_snap> <dim>

"""

import numpy as np
from docopt import docopt
import h5py as h5
from astropy.utils.console import ProgressBar

__author__ = "Simon Mutch"
__date__   = "2013/05/24"

args = docopt(__doc__)

# cmd line args
SIM = args['<sim>']
FIRST_SNAP = int(args['<first_snap>'])
LAST_SNAP = int(args['<last_snap>'])
DIM = int(args['<dim>'])
DIR = "/home/gpoole/GiggleZ/{:s}/grids/".format(SIM)

cell_locs = np.ones((LAST_SNAP-FIRST_SNAP+1, 3), 'i4')*-1

# loop through each snapshot
with ProgressBar(LAST_SNAP-FIRST_SNAP+1) as bar:
    for ii, snap in enumerate(xrange(FIRST_SNAP, LAST_SNAP+1)):
        with open("{:s}snapshot_{:03d}_dark_grid.dat".format(DIR,snap), "rb") as fin:
            # read the header info
            n_cell = np.fromfile(fin, 'i4', 3)
            grid_size = np.fromfile(fin, 'f8', 3)
            n_grids = np.fromfile(fin, 'i4', 1)[0]
            ma_scheme = np.fromfile(fin, 'i4', 1)[0]

            # read in the identifier
            ident = np.fromfile(fin, '|S32', 1)

            # read in the grid
            grid = np.fromfile(fin, '<f4', DIM**3)
            grid.shape = [DIM,]*3

            # find the most dense cell
            w = np.array(np.where(grid==grid.max()), 'i4')

            # save the result
            cell_locs[ii,:] = w[:,0]
        bar.update()

with h5.File("locs.hdf5", "w") as fout:
    fout.create_dataset("cell_locs", data=cell_locs)

