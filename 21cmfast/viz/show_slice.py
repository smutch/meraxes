#!/usr/bin/env python

""" Usage: show_slice.py <input_file> <output_file> <slice> [--expansion_factor=<val>] [--color_bar] [--dt_type=<type>] [--cmap=<name>]

    --expansion_factor=<val>  Ratio of final image to input array dims [default: 1].
    --color_bar               Show color bar.
    --dt_type=<type>          char or float type for output box [default: float]
    --cmap=<name>             Colormap name [default: BuPu]

    Args:
        input_file            Input filename
        output_file           Ouput filename
        slice                 Box slice (e.g. [:,10,:] or [:,:,5])

"""

import sys, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docopt import docopt

def setup_fig(redshift, neutral_fraction, slice_str):

    fig = plt.figure(0)
    fig.clf()
    ax = plt.subplot(111)

    axis_names = np.array(('x','y','z'))
    s = np.array(parse_slice(slice_str))
    labels = axis_names[s==slice(None, None, None)]
    ax.set_xlabel(labels[0]+' (Mpc)')
    ax.set_ylabel(labels[1]+' (Mpc)')
    plt.title('z=%.2f; nf=%.3f; slice=[%s]' % (redshift, neutral_fraction,
                                             slice_str))
    return fig, ax


def plot_slice(slice_img, ax, dim, box_size, color_bar=False, cmap='jet'):
   
    final_size = (slice_img.shape[0]/dim)*box_size

    cax = ax.imshow(slice_img.T, interpolation='nearest',
                    cmap=getattr(plt.cm, cmap),
                    extent=(0,final_size,0,final_size),
                    origin='lower')

    if color_bar:
        cbar = fig.colorbar(cax)
        cbar.set_label('mK')

def parse_filename(fname):

    pattern = r'.*_z(?P<redshift>\d+\.\d+)_nf(?P<neutral_fraction>\d+\.\d+).*_(?P<dim>\d+)_(?P<box_size>\d+)Mpc'
    params = re.match(pattern, fname).groupdict()

    return (float(params['redshift']),
            float(params['neutral_fraction']), int(params['dim']),
            float(params['box_size']))
    

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
   
    args = docopt(__doc__)

    redshift, neutral_fraction, dim, box_size = parse_filename(args['<input_file>'])
    slice_str = args['<slice>']

    expansion_factor = float(args['--expansion_factor'])

    # read the binary
    if args['--dt_type']=='char':
        dt_type='<S1'
    else:
        dt_type='<f4'
    box = np.fromfile(args['<input_file>'], dtype=np.dtype((dt_type,(dim, dim, dim))), count=1)[0]
   
    slice_img = box[parse_slice(slice_str)].squeeze()
    min_val = np.min(slice_img)
    max_val = np.max(slice_img)

    print "min_val, max_val = %.3f, %.3f" % (min_val, max_val)

    slice_dim = slice_img.shape

    fig, ax = setup_fig(redshift, neutral_fraction, slice_str)
    plot_slice(slice_img, ax, dim, box_size,
               color_bar=args['--color_bar'],
               cmap=args['--cmap'])

    plt.tight_layout()
    plt.savefig(args['<output_file>'])
