#cython: boundscheck=False
#cython: wraparound=False

from __future__ import division
import numpy as np
from astropy.utils.console import ProgressBar
cimport numpy as np

def regrid(np.ndarray[np.float32_t, ndim=3] old_grid not None,
           float box_size,
           int new_dim):

    cdef int old_dim = old_grid.shape[0]
    cdef float resample_factor = float(new_dim/old_dim)
    cdef float old_cell_volume = (box_size/float(old_dim))**3
    cdef float new_cell_volume = (box_size/float(new_dim))**3

    cdef np.ndarray[np.float32_t, ndim=3] new_grid = np.zeros([new_dim,new_dim,new_dim], np.float32)

    cdef unsigned int i,j,k, rsi, rsj, rsk
    with ProgressBar(old_dim) as bar:
        for i in xrange(old_dim):
            for j in xrange(old_dim):
                for k in xrange(old_dim):
                    rsi = int(i*resample_factor)
                    rsj = int(j*resample_factor)
                    rsk = int(k*resample_factor)
                    new_grid[rsi, rsj, rsk] += old_grid[i,j,k]
            bar.update()

    new_grid *= old_cell_volume / new_cell_volume

    return new_grid
