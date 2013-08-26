import numpy as np
from yt.mods import *
from yt.frontends.stream.api import load_uniform_grid

# Load data
dim = 256
n_elem = dim**3
grid = np.fromfile("./xH_z005.84_nf0.178367_eff31.5_HIIfilter1_Mmin1.7e+08_RHIImax30_256_177Mpc", dtype='<f4', count=n_elem) 
grid.shape = [dim,dim,dim]

# Loading data into yt structure
data = dict(Density = grid)
pf = load_uniform_grid(data, grid.shape, 1, nprocs=2)

# Setting up parameters for volume rendering. See the following 
# link for more details on the parameters:
# http://yt-project.org/doc/cookbook/simple_plots.html#cookbook-simple-volume-rendering
tf = ColorTransferFunction((0.,1.))
# tf.add_layers(1, colormap='Purples')
tf.sample_colormap(1., 0.5, colormap='PuRd')
c = [0.5, 0.5, 0.5]
W = 1.5
Nvec = 512
L = [0.1, 0.1, 1.0]
cam = pf.h.camera(c, L, W, Nvec, tf)
# im = cam.snapshot()
# im = cam.draw_domain(im)
# im.write_png('test.png')

for i, snapshot in enumerate(cam.rotation(2.*np.pi, 360)):
    im = cam.snapshot()
    im = cam.draw_domain(im)
    im.write_png('rotation_%04d.png' % i)
