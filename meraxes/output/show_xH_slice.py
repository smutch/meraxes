#!/usr/bin/env python

""" Usage:
        show_xH_slice.py <input_file> <snapshot> <output_file> <slice> [--color_bar] [--cmap=<name>] [--galaxies]

    Options:
        --color_bar           Show color bar.
        --cmap=<name>         Colormap name [default: gist_earth_r]
        --galaxies            Overplot the galaxies

    Args:
        input_file            Input filename
        snapshot              Snapshot to plot
        output_file           Ouput filename
        slice                 Box slice (e.g. [:,10,:] or [:,:,5])
"""

import sys, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docopt import docopt
import h5py as h5

from dragons import meraxes, nbody, plotutils

def setup_fig(redshift, neutral_fraction, slice_str, slice_axis):

    fig, ax = plt.subplots(1,1, facecolor='w')
    ax.grid('off')
    # ax.grid(color='0.8', alpha=0.5)

    axis_names = np.array(('x','y','z'))
    sel = np.ones(3, bool)
    sel[slice_axis] = False
    labels = np.array(['x', 'y', 'z'])[sel]
    ax.set_xlabel(labels[0]+' (Mpc/h)')
    ax.set_ylabel(labels[1]+' (Mpc/h)')
    plt.title('z=%.2f; nf=%.3f; slice=[%s]' % (redshift, neutral_fraction,
                                             slice_str),
              size="medium")
    return fig, ax

def plot_slice(slice_img, ax, dim, slice_axis, box_size, galaxies=False, dm=False, color_bar=False, cmap='gist_earth_r'):

    final_size = (slice_img.shape[0]/dim)*box_size
    extent = (0,final_size,0,final_size)

    if dm is not False:
        ax.imshow(np.log10(dm.T),
                  interpolation='bilinear',
                  cmap=plt.cm.Spectral,
                  extent = extent,
                  origin='lower')

    alpha_cmap = plt.cm.get_cmap("Blues_r")
    alpha_cmap._init()
    alphas = np.abs(np.linspace(0.7, 0.0, alpha_cmap.N))
    alpha_cmap._lut[:-3,-1] = alphas
    cax = ax.imshow(1.0-slice_img.T, interpolation='bilinear',
                    cmap=alpha_cmap,
                    extent = extent,
                    origin='lower',
                    vmin=0,
                    vmax=1)

    ax.contour(slice_img.T, extent=extent, extend='both',
               linewidths=0.5, alpha=0.5,
               cmap=plt.cm.Greys_r)

    if color_bar:
        cbar = fig.colorbar(cax,)
        cbar.set_label(r'${\rm x_H}$')

    if galaxies is not False:
        i_axis = np.arange(3, dtype=int)
        plot_axis = np.argwhere(i_axis!=slice_axis).squeeze()
        ax.scatter(galaxies["Pos"][:,plot_axis[0]],
                   galaxies["Pos"][:,plot_axis[1]],
                   s=np.log10(galaxies['StellarMass']*1.e10)**9*3e-7,
                   c=np.log10(galaxies['StellarMass']*1.e10),
                   cmap=plt.cm.Blues,
                   marker='o',
                   # edgecolors='none',
                   linewidths=0.1,
                   alpha=0.5,
                   zorder=3)

    ax.set_xlim(0,final_size)
    ax.set_ylim(0,final_size)

def parse_slice(slice_str, max_dim):

    pattern = r'(?P<sx>.+),(?P<sy>.+),(?P<sz>.+)'
    axis_pattern = r'(?P<from>[0-9]+)?:?(?P<to>[0-9]+)?'
    slice_dict = re.match(pattern, slice_str).groupdict()

    default = [0,max_dim,None]
    for ii, k in enumerate(slice_dict.iterkeys()):
        if slice_dict[k]==':':
            slice_dict[k]=default[:]
        else:
            axis_dict = re.match(axis_pattern, slice_dict[k]).groupdict()
            slice_dict[k]=default[:]
            if axis_dict['from'] is not None:
                slice_dict[k][0]=int(axis_dict['from'])
            else:
                slice_dict[k][0]=0
            if axis_dict['to'] is not None:
                slice_dict[k][1]=int(axis_dict['to'])
            else:
                slice_dict[k][1]=slice_dict[k][0]+1

    # for i in range(3):
    #     if slice_sel[i].start is None:
    #         slice_sel[i]=slice(0,slice_sel[i].stop,None)
    #     if slice_sel[i].stop is None:
    #         slice_sel[i]=slice(slice_sel[i].start,slice_dim,None)

    s = [slice(*slice_dict['sx']),
         slice(*slice_dict['sy']),
         slice(*slice_dict['sz'])]

    return s


if __name__ == '__main__':

    # deal with the input args
    args = docopt(__doc__)
    input_file = args['<input_file>']
    snapshot = int(args['<snapshot>'])
    output_file = args['<output_file>']
    redshift = meraxes.io.grab_redshift(input_file, snapshot)
    simprops = meraxes.io.read_input_params(input_file)

    # read the grid and parse the requested slice
    grid = meraxes.io.read_grid(input_file, snapshot, "xH")
    dm = meraxes.io.read_grid(input_file, snapshot, "deltax")
    dm /= dm.mean()
    slice_dim = grid.shape[0]
    slice_sel = parse_slice(args['<slice>'], slice_dim)
    box_len = simprops["BoxSize"]
    with h5.File(input_file, "r") as fd:
        global_xH = fd["Snap{:03d}/Grids/xH".format(snapshot)].attrs["global_xH"][0]

    # if requested read in the galaxies and select those within our slice
    if(args['--galaxies']):
        gals = meraxes.io.read_gals(input_file, snapshot)
        edges = np.linspace(0, box_len, slice_dim+1)
        sel = np.prod([((gals["Pos"][:,i] <= edges[slice_sel[i].stop]) &
                     (gals["Pos"][:,i] > edges[slice_sel[i].start]))
                    for i in range(3) ], axis=0)
        gals = gals[np.nonzero(sel)]
    else:
        gals = False

    # Select out the grid slice
    grid_slice = grid[slice_sel]
    slice_axis = np.argmin(grid_slice.shape)
    grid_slice = grid_slice.squeeze()
    dm_slice = dm[slice_sel].squeeze()
    if grid_slice.ndim==3:
        grid_slice = grid_slice.mean(axis=slice_axis)
        dm_slice = dm_slice.mean(axis=slice_axis)

    plotutils.init_style(theme="white_bg")
    fig, ax = setup_fig(redshift, global_xH, args['<slice>'], slice_axis)
    plot_slice(grid_slice, ax, slice_dim, slice_axis, box_len,
               galaxies=gals,
               dm=dm_slice,
               color_bar=args['--color_bar'],
               cmap=args['--cmap'])

    plt.tight_layout()
    plt.savefig(args['<output_file>'], dpi=200)
