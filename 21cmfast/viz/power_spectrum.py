#!/usr/bin/env python

""" Plot the power spectrum of a cube.

    Usage: power_spectrum.py <input_file> <output_file> 

    Args:
        input_file            Input filename
        output_file           Ouput filename

"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from show_slice import parse_filename

__author__ = "Simon Mutch"
__date__   = "2013/05/22"

def calc_power_spectrum(box, size, klim=None, nbins=20):
    """Calculate the power spectrum of a cube."""

    # Calculate the 3d power spectrum
    ps_3d = np.abs(np.fft.fftn(box))**2.
    cell_size = size/box.shape[0]
    volume = size**3.
    k_mgrid = np.meshgrid(*([np.fft.fftfreq(box.shape[0], cell_size),]*3))
    k_mag = np.sqrt(k_mgrid[0]**2 + k_mgrid[1]**2 + k_mgrid[2]**2)

    # Convert to 'radial' binning
    if klim==None:
        k_edges = np.linspace(0, k_mag.max(), nbins+1)
    else:
        k_edges = np.linspace(klim[0], klim[1], nbins+1)
    hdk = (k_edges[1]-k_edges[0])/2.
    k = np.zeros(nbins, 'f8')
    ps = np.zeros(nbins, 'f8')
    for annulus in xrange(nbins):
        k[annulus] = k_edges[annulus]+hdk
        sel = ((k_mag>=k_edges[annulus]) & (k_mag<k_edges[annulus+1]))
        ps[annulus] = (k[annulus]**3 /volume) * ps_3d[sel].mean()

    return ps, k


def setup_fig(redshift, neutral_fraction):

    fig = plt.figure(0)
    fig.clf()
    ax = plt.subplot(111)

    ax.set_ylabel(r"$\Delta^2({\bf k})$")
    ax.set_xlabel(r"$\rm k\cdot h\ (Mpc^{-1})$")
    plt.title('z=%.2f; nf=%.3f' % (redshift, neutral_fraction))
    return fig, ax


if __name__ == '__main__':
    
    # Deal with command line arguments etc.
    args = docopt(__doc__)
    INPUT_FILENAME = args['<input_file>']
    OUTPUT_FILENAME = args['<output_file>']
    redshift, neutral_fraction, dim, box_size = parse_filename(INPUT_FILENAME)

    # Read in the box
    box = np.fromfile(INPUT_FILENAME, '<f4', count=dim**3)
    box.shape = [dim,]*3

    # Convert to "overdensity"
    mean = box.mean()
    box = box/mean -1.

    # Calculate the power spectrum
    ps, k = calc_power_spectrum(box, box_size)

    fig, ax = setup_fig(redshift, neutral_fraction)
    ax.loglog(k, ps)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILENAME)
