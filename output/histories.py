#!/usr/bin/env python

"""Plot histories...

Usage: histories.py <master_file> <gal_id> <last_snapshot> [--output=<dir_path> --format=<ext>]

Options:
    --output=<dir_path>   Target directory for output figures [default: ./plots]
    --format=<ext>        Image format (e.g. png, pdf) [default: png]

"""

from os import path
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt

from dragons import meraxes, munge, plotutils

__author__ = "Simon Mutch"
__date__   = "2014-03-17"


def plot_histories(fname, gal_id, last_snapshot, output_dir, fig_format):

    # read in the galaxy history
    gal_hist = meraxes.galaxy_history(fname, gal_id, last_snapshot, pandas=True)

    # initialise the dragons plotting style
    plotutils.init_style(context='inline')

    # create the figure & axes
    n_axes = 7
    fig, axes = plt.subplots(n_axes, 1, figsize=(16,3*n_axes))

    snapshots = np.arange(last_snapshot+1)

    # baryonic total mass components
    ax = axes[0]
    ax.plot(snapshots, np.log10(1e10 * gal_hist.StellarMass), label="stellar")
    ax.plot(snapshots, np.log10(1e10 * gal_hist.ColdGas), label="cold")
    ax.plot(snapshots, np.log10(1e10 * gal_hist.HotGas), label="hot")
    ax.plot(snapshots, np.log10(1e10 * gal_hist.EjectedGas), label="ejected")
    ax.set_xlim(0, last_snapshot)
    ax.set_xlabel("snapshot")
    ax.set_ylabel(r"$\log_{10}(M/{\rm M_{\odot}})$")
    ax.legend(loc="upper left", frameon=True)

    # metallicity
    ax = axes[1]
    ax.plot(snapshots, gal_hist.MetalsStellarMass/gal_hist.StellarMass, label="stellar")
    ax.plot(snapshots, gal_hist.MetalsColdGas/gal_hist.ColdGas, label="cold")
    ax.plot(snapshots, gal_hist.MetalsHotGas/gal_hist.HotGas, label="hot")
    ax.plot(snapshots, gal_hist.MetalsEjectedGas/gal_hist.EjectedGas, label="ejected")
    ax.axhline(0.02, color="SlateGrey", ls="--", label="solar")
    ax.set_xlim(0, last_snapshot)
    ax.set_xlabel("snapshot")
    ax.set_ylabel(r"$Z_{\rm total}$")
    ax.legend(loc="upper left", frameon=True)

    # sfr
    ax = axes[2]
    ax.plot(snapshots, gal_hist.Sfr)
    ax.set_xlim(0, last_snapshot)
    ax.set_xlabel("snapshot")
    ax.set_ylabel(r"${\rm SFR\ (M_{\odot}/yr^{-1})}$")

    # cooling mass
    ax = axes[3]
    ax.plot(snapshots, np.log10(1.e10 * gal_hist.Mcool))
    ax.set_xlim(0, last_snapshot)
    ax.set_xlabel("snapshot")
    ax.set_ylabel(r"$\log_{10}(M_{\rm cool}/{\rm M_{\odot}})$")

    # Mvir
    ax = axes[4]
    ax.plot(snapshots, np.log10(1.e10 * gal_hist.Mvir))
    ax.set_xlim(0, last_snapshot)
    ax.set_xlabel("snapshot")
    ax.set_ylabel(r"$\log_{10}(M_{\rm vir}/{\rm M_{\odot}})$")

    # Type / ghost flag
    ax = axes[5]
    ax.plot(snapshots, gal_hist.Type, alpha=0.6, label="type")
    ax.plot(snapshots, gal_hist.GhostFlag, alpha=0.6, label="ghost flag")
    ax.set_xlim(0, last_snapshot)
    ax.set_ylim(-1, 3)
    ax.set_xlabel("snapshot")
    ax.set_ylabel("value")
    ax.legend(loc="upper right", frameon=True)

    # gas metallicity [O/H]
    ax = axes[6]
    s = ax.scatter(np.log10(1.e10*gal_hist.StellarMass),
                   np.log10(gal_hist.MetalsColdGas/gal_hist.ColdGas/0.02)+9.0, alpha=0.8,
                   s=60, edgecolors='none', cmap=plt.cm.Blues,
                   c=np.arange(100))
    cax = plt.colorbar(s, ax=ax)
    cax.set_label("snapshot")
    ax.set_xlabel(r"$\log_{10}(M_*/{\rm M_{\odot}})$")
    ax.set_ylabel(r"$\log_{10}[{\rm O/H}] + 12$")

    # save
    plt.tight_layout()
    plt.savefig(path.join(output_dir,
                          "galaxy_history_ID{:d}.{:s}".format(gal_id,
                                                             fig_format)))


if __name__ == '__main__':
    args = docopt(__doc__)

    plot_histories(args['<master_file>'], int(args['<gal_id>']),
                   int(args['<last_snapshot>']), args['--output'],
                   args['--format'])
