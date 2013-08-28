#!/usr/bin/env python

"""Usage: show_box.py <input_file> <dim>

"""

import numpy as np
from mayavi import mlab
from docopt import docopt

__author__ = "Simon Mutch"
__date__   = "13 Mar 2013"

args = docopt(__doc__)

dim = int(args['<dim>'])
box = np.fromfile(args['<input_file>'], dtype=np.dtype(('<f4',(dim, dim, dim))), count=1)[0]
print "Box size = %.2fMB"%(box.nbytes/1024/1024)
# box = np.log10(box)
# box[np.isinf(box)]=0.

# box_slice = box[10:40,10:40,10:40].copy()
# del(box)
mlab.contour3d(box, contours=10, transparent=True, opacity=0.22)
mlab.show()
