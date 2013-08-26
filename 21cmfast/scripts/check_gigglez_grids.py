#!/usr/bin/env python

"""Check the GiggleZ density grids and compare with theoretical
expectations.

Usage: check_gigglez_density_grids.py <fname> <particle_mass>
"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt

__author__ = "Simon Mutch"
__date__   = "2013/07/04"

cosmology = {
    "omega_l"    : 0.727,
    "omega_m"    : 0.273,
    "omega_k"    : 0.,
    "omega_b"    : 0.0456,
    "h"          : 0.705,
    "sigma_8"    : 0.812,
    "n_spectral" : 0.960,
}

def read_grid(fname):
    """Read in the density grid."""

    # print
    # print "Reading grid..."
    with open(fname, "rb") as fin:
        # read the header info
        n_cell = np.fromfile(fin, 'i4', 3)
        box_size = np.fromfile(fin, 'f8', 3)
        n_grids = np.fromfile(fin, 'i4', 1)[0]
        ma_scheme = np.fromfile(fin, 'i4', 1)[0]
        # read in the identifier
        ident = np.fromfile(fin, 'S32', 1)
        # read in the grid
        # N.B. We are assuming that the first grid here is the density grid
        grid = np.fromfile(fin, '<f4', n_cell.cumprod()[-1])
        grid.shape = n_cell
    # print "...done"

    return grid, box_size, n_cell


def main(fname, particle_mass):

    print
    print "File = {:s}".format(fname)
    print "Particle mass = {:.6e}/h Msol".format(particle_mass)

    grid, box_size, n_cell = read_grid(fname)

    print
    print "RESULTS"
    print "-------"

    # mean density
    mean = grid.mean()
    print "{:40s} : {:.2e}h^2 Msol Mpc^-3".format("mean comoving density",
                                                  mean*1.e10)

    # theoretical mean comoving density
    Hubble = 100.  # km/s/Mpc
    Grav = 6.67384e-11 # m^3 kg^-1 s^-2
    rho_crit = 3.*Hubble*Hubble/(8.*np.pi*Grav)*0.015514351  # h^2 Msol Mpc^-3
    rho_mean = cosmology["omega_m"]*rho_crit  # h^2 Msol Mpc^-3
    print "{:40s} : {:.2e}h^2 Msol Mpc^-3".format("THEORETICAL mean comoving density"
                                                  ,rho_mean)

    print "{:40s} : {:.2f}".format("fraction",mean*1.e10/rho_mean)
    print

    # total mass
    cell_volume = (box_size[0]/float(n_cell[0]))**3
    mass = grid.sum()*cell_volume
    print "{:40s} : {:.2e}/h Msol".format("total mass",mass*1.e10)

    # number of particles
    npart = int(mass/particle_mass*1.e10)
    print "{:40s} : {:d}^3".format("number of particles",
                                   int(npart**(1./3.)))

    # max
    max_npart = int(grid[grid==grid.max()]*cell_volume/particle_mass*1.e10)
    print "{:40s} : {:.2e}h^2 Msol Mpc^-3".format("maximum density",grid.max()*1.e10)
    print "{:40s} : {:d}".format("maximum particle number",max_npart)
    
    print


if __name__ == '__main__':
    args = docopt(__doc__)
    fname = args['<fname>']
    particle_mass = float(args['<particle_mass>'])
    main(fname, particle_mass)
