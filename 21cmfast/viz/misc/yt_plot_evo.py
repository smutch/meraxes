import numpy as np
from yt.mods import *
from yt.frontends.stream.api import load_uniform_grid
from glob import glob

# Load data
files = np.sort(np.array(glob('../GiggleZ/Boxes/xH*')))
dim = 256
n_elem = dim**3
grid = np.fromfile(files[0], dtype='<f4', count=n_elem) 
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
Nvec = 100 #512
L = [0.1, 0.1, 1.0]
cam = pf.h.camera(c, L, W, Nvec, tf)
# im = cam.snapshot()
# im = cam.draw_domain(im)
# im.write_png('frame_%04d.png' % 0)

n_files = len(files)
frames_per_file = 2
n_frames = n_files*frames_per_file
i_file = 0
for i, angle in enumerate(np.linspace(0, 2.*np.pi, n_frames)):
    if i>0:
        if i%frames_per_file==0:
            i_file+=1
            grid = np.fromfile(files[i_file], dtype='<f4', count=n_elem) 
            grid.shape = [dim,dim,dim]
            data = dict(Density = grid)
            pf = load_uniform_grid(data, grid.shape, 1, nprocs=2)
            cam = pf.h.camera(c, L, W, Nvec, tf)
            cam.rotate(angle)
        else:
            cam.rotate(angle-last_angle)
    im = cam.snapshot()
    im = cam.draw_domain(im)
    im.write_png('frame_%04d.png' % i)
    last_angle = angle

