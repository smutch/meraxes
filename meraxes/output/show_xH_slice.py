#!/usr/bin/env python

""" Usage: 
        show_xH_slice.py <input_file> <snapshot> <output_file> <slice> [--color_bar] [--cmap=<name>]

    Options:
        --color_bar           Show color bar.
        --cmap=<name>         Colormap name [default: BuPu]

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

from ssimpl import meraxes, nbody

def setup_fig(redshift, neutral_fraction, slice_str):

    fig, ax = plt.subplots(1,1, facecolor='w')
    ax.grid('off')
    # ax.grid(color='0.8', alpha=0.5)

    axis_names = np.array(('x','y','z'))
    s = np.array(parse_slice(slice_str))
    labels = axis_names[s==slice(None, None, None)]
    ax.set_xlabel(labels[0]+' (Mpc)')
    ax.set_ylabel(labels[1]+' (Mpc)')
    plt.title('z=%.2f; nf=%.3f; slice=[%s]' % (redshift, neutral_fraction,
                                             slice_str))
    return fig, ax


def plot_slice(slice_img, ax, dim, box_size, galaxies=None, color_bar=False, cmap='jet'):
   
    final_size = (slice_img.shape[0]/dim)*box_size
    extent = (0,final_size,0,final_size)

    cax = ax.imshow(slice_img.T, interpolation='nearest',
                    cmap=getattr(plt.cm, cmap),
                    extent = extent,
                    origin='lower')

    ax.contour(slice_img.T, extent=extent, extend='both',
               linewidths=1, alpha=0.5, cmap=getattr(plt.cm, cmap+'_r'))

    if color_bar:
        cbar = fig.colorbar(cax)
        cbar.set_label(r'${\rm x_H}$')


def parse_slice(slice_str):

    pattern = r'(?P<sx>.+?),(?P<sy>.+?),(?P<sz>.+?)'
    slice_dict = re.match(pattern, slice_str).groupdict()

    default = [None,None,None]
    for ii, k in enumerate(slice_dict.iterkeys()):
        if slice_dict[k]==':':
            slice_dict[k]=default[:]
        else:
            tmp = int(slice_dict[k])
            slice_dict[k]=default[:]
            slice_dict[k][0]=tmp
            slice_dict[k][1]=tmp+1

    s = (slice(*slice_dict['sx']),
         slice(*slice_dict['sy']),
         slice(*slice_dict['sz']))

    return s

if __name__ == '__main__':
   
    # deal with the input args
    args = docopt(__doc__)
    input_file = args['<input_file>']
    snapshot = int(args['<snapshot>'])
    output_file = args['<output_file>']
    slice_sel = parse_slice(args['<slice>'])
    redshift = meraxes.io.grab_redshift(input_file, snapshot)

    # read the grid
    grid, grid_props = meraxes.io.read_xH_grid(input_file, snapshot)
    grid_slice = grid[slice_sel].squeeze()
    slice_dim = grid_slice.shape[0]

    # if requested also read in the galaxies
    # if(args['--galaxies']):
    #     gals = meraxes.io.read_gals(input_file, snapshot)
    # else:
    gals = None
    

    fig, ax = setup_fig(redshift, grid_props['global_xH'], args['<slice>'])
    plot_slice(grid_slice, ax, slice_dim, grid_props['box_len'],
               galaxies=gals,
               color_bar=args['--color_bar'],
               cmap=args['--cmap'])

    plt.tight_layout()
    plt.savefig(args['<output_file>'])
