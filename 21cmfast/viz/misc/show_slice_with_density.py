#!/usr/bin/env python

""" Usage: show_slice.py <input_file> <output_file> <slice> <slice_depth> <sim> <snapshot> [--color_bar] [--dt_type=<type>]

    --color_bar               Show color bar.
    --dt_type=<type>          char or float type for output box [default: float]

    Args:
        input_file            Input filename
        output_file           Ouput filename
        slice                 Box slice centre (e.g. [:,10,:] or [:,:,5])
        slice_depth           +/- slice_depth cells

"""

import sys, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docopt import docopt

import pyximport
pyximport.install()
from regrid import regrid

HUBBLE_h = 0.705

def setup_fig(redshift, neutral_fraction, slice_str, slice_depth):

    fig = plt.figure(0)
    fig.clf()
    ax = plt.subplot(111)

    axis_names = np.array(('x','y','z'))
    s = np.array(parse_slice(slice_str, slice_depth))
    labels = axis_names[s==slice(None, None, None)]
    ax.set_xlabel(labels[0]+' $(h^{-1}\,Mpc)$')
    ax.set_ylabel(labels[1]+' $(h^{-1}\,Mpc)$')
    plt.title('z=%.2f; nf=%.3f; slice=[%s]+-%d' % (redshift, neutral_fraction,
                                                   slice_str, slice_depth))
    return fig, ax


def plot_slice(slice_img, ax, dim, box_size, color_bar=False):
   
    final_size = (float(grid_slice.shape[0])/float(dim))*box_size

    cax = ax.imshow(1.-slice_img.T, interpolation='nearest',
                    cmap=getattr(plt.cm, 'Reds'),
                    extent=(0,final_size,0,final_size),
                    origin='lower',
                    vmin=0.,
                    vmax=1.,
                    alpha=0.3)

    if color_bar:
        cbar = fig.colorbar(cax)
        cbar.solids.set_edgecolor("face")
        cbar.set_label('ionised fraction')

def plot_dgrid(grid_slice, ax, dim, box_size):
   
    final_size = (float(grid_slice.shape[0])/float(dim))*box_size

    cax = ax.imshow(np.log10(grid_slice.T*1e10), interpolation='bilinear',
                    cmap=getattr(plt.cm, 'Greys'),
                    extent=(0,final_size,0,final_size),
                    origin='lower')


def parse_filename(fname):

    pattern = r'.*_z(?P<redshift>\d+\.\d+)_nf(?P<neutral_fraction>\d+\.\d+).*_(?P<dim>\d+)_(?P<box_size>\d+)Mpc'
    params = re.match(pattern, fname).groupdict()

    return (float(params['redshift']),
            float(params['neutral_fraction']), int(params['dim']),
            float(params['box_size']))
    

def parse_slice(slice_str, slice_depth):

    pattern = r'(?P<sx>.+?),(?P<sy>.+?),(?P<sz>.+?)$'
    slice_dict = re.match(pattern, slice_str).groupdict()

    default = [None,None,None]
    for ii, k in enumerate(slice_dict.iterkeys()):
        if slice_dict[k]==':':
            slice_dict[k]=default[:]
        else:
            tmp = int(slice_dict[k])
            slice_dict[k]=default[:]
            slice_dict[k][0]=tmp-slice_depth
            slice_dict[k][1]=tmp+slice_depth+1

    s = (slice(*slice_dict['sx']),
         slice(*slice_dict['sy']),
         slice(*slice_dict['sz']))

    return s


if __name__ == '__main__':
   
    args = docopt(__doc__)

    redshift, neutral_fraction, dim, box_size = parse_filename(args['<input_file>'])
    slice_str = args['<slice>']
    slice_depth = int(args['<slice_depth>'])

    sim_dir = "/home/smutch/Tiamat/{:s}/grids/".format(args['<sim>'])
    snapshot = int(args['<snapshot>'])

    # read the binary
    if args['--dt_type']=='char':
        dt_type='<S1'
    else:
        dt_type='<f4'
    box = np.fromfile(args['<input_file>'], dtype=np.dtype((dt_type,(dim, dim, dim))), count=1)[0]
  
    # read the density grid
    # print "Reading grid..."
    with open("{:s}grid_nompi_{:d}_1024_dark_grid.dat".format(sim_dir,snapshot), "rb") as fin:
        # read the header info
        n_cell = np.fromfile(fin, 'i4', 3)
        box_size_grid = np.fromfile(fin, 'f8', 3)
        n_grids = np.fromfile(fin, 'i4', 1)[0]
        ma_scheme = np.fromfile(fin, 'i4', 1)[0]
        # read in the identifier
        ident = np.fromfile(fin, 'S32', 1)
        # read in the grid
        grid = np.fromfile(fin, '<f4', n_cell.cumprod()[-1])
        grid.shape = n_cell
    # print "...done"

    # resample the grid if necessary
    print "Resampling density grid to match ionization grid..."
    grid = regrid(grid, box_size_grid[0], dim)
    print "...done"

    slice_img = box[parse_slice(slice_str, slice_depth)].squeeze()
    grid_slice = grid[parse_slice(slice_str, slice_depth)].squeeze()

    # If we have to, sum through the slice
    if(slice_depth>0):
        i_axis = np.where(slice_img.shape == np.array(slice_img.shape, int).min())[0]
        slice_img = slice_img.sum(axis=i_axis)/float(slice_depth*2+1)
        grid_slice = grid_slice.sum(axis=i_axis)

    min_val = np.min(grid_slice)
    max_val = np.max(grid_slice)

    # print "min_val, max_val = %.3f, %.3f" % (min_val, max_val)

    slice_dim = slice_img.shape

    fig, ax = setup_fig(redshift, neutral_fraction, slice_str, slice_depth)
    plot_dgrid(grid_slice, ax, dim, box_size_grid[0])
    plot_slice(slice_img, ax, dim, box_size_grid[0],
               color_bar=args['--color_bar'])

    plt.tight_layout()
    plt.savefig(args['<output_file>'], dpi=plt.rcParams['savefig.dpi']*2)
