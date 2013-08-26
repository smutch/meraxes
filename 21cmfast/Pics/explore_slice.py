#!/usr/bin/env python

""" Usage: explore.py <input_file> [--expansion_factor=<val>] [--color_bar] [--dt_type=<type>] [--cmap=<name>]

    --expansion_factor=<val>  Ratio of final image to input array dims [default: 1].
    --color_bar               Show color bar.
    --dt_type=<type>          char or float type for output box [default: float]
    --cmap=<name>             Colormap name [default: BuPu]

    Args:
        input_file            Input filename

"""

import numpy as np
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from docopt import docopt

from show_slice import setup_fig, plot_slice, parse_filename

__author__ = "Simon Mutch"
__date__   = "23 Mar 2013"

args = docopt(__doc__)

redshift, neutral_fraction, dim, box_size = parse_filename(args['<input_file>'])

expansion_factor = float(args['--expansion_factor'])

# read the binary
if args['--dt_type']=='char':
    dt_type='<S1'
else:
    dt_type='<f4'
box = np.fromfile(args['<input_file>'], dtype=np.dtype((dt_type,(dim, dim, dim))), count=1)[0]

slice_img = box[:,:,int(dim/2)].squeeze()
box_min = np.min(box)
box_max = np.max(box)

# Plot
fig, ax_slice = setup_fig(redshift, neutral_fraction, ":,:,1")
plt.title('z=%.2f; nf=%.3f' % (redshift, neutral_fraction))
cax = plot_slice(slice_img, ax_slice, dim, box_size,
                 cmap=args['--cmap'])

if args['--color_bar']:
    cbar = fig.colorbar(cax)

# Add the slider
fig.subplots_adjust(bottom=0.2)
ax_slider = plt.axes((0.1,0.1,0.8,0.02)) 

def slider_action(val):
    slice_img = box[:,:,val].squeeze()
    plot_slice(slice_img, ax_slice, dim, box_size,
               cmap=args['--cmap'])
slider = Slider(ax_slider, "z", 0, dim, valinit=int(dim/2), color='#AAAAAA')
slider.on_changed(slider_action)

plt.show(block=True)
