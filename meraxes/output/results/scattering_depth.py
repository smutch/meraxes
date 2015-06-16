#!/usr/bin/env python

"""Plot the electron scattering depth for a full 21cmFAST coupled run.

Usage: scattering_depth.py <fname>

"""

import numpy as np
import matplotlib.pyplot as plt
from dragons import meraxes, plotutils
import h5py as h5
from docopt import docopt
import os

__author__ = "Simon Mutch"
__date__ = "2014-12-08"

PLOT_RANGE = [6, 20, 0.02, 0.08]
HUBBLE = 0.7


def plot(model, ax):

    ax.plot(model[0], model[1], label="Meraxes", lw=4)

    # Planck 2015 constraint
    best, uncert = 0.066, 0.012
    ax.axhline(best, color='0.3', zorder=0, ls='-', lw=1)
    ax.fill_between([PLOT_RANGE[0], PLOT_RANGE[1]],
                    best-uncert, best+uncert,
                    color='0.5', alpha=0.3, zorder=0)

    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"$\tau_e$")

    ax.text(PLOT_RANGE[0]+0.1, 0.067, "Planck 2015", horizontalalignment="left",
            verticalalignment="bottom")

    ax.axis(PLOT_RANGE)
    ax.legend(loc="lower right", ncol=2)


if __name__ == '__main__':
    
    args = docopt(__doc__)

    fname = args['<fname>']

    # set up
    meraxes.set_little_h(HUBBLE)
    plt.style.use(["dragons", "white_background"])
    fig, ax = plt.subplots(1, 1)
    ax.grid('off')

    model = meraxes.electron_optical_depth(fname)
    plot(model, ax)

    # save
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(fname), "plots/scattering_depth.pdf"))
