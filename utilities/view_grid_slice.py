"""
TODO
====

- add command line wrappers
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm
from astropy.utils.console import ProgressBar

def read_grid_header(fin):
    """Read in the grid file header from the open file handle `fin`."""
    
    n_cell = np.fromfile(fin, 'i4', 3)
    box_size = np.fromfile(fin, 'f8', 3)
    n_grids = np.fromfile(fin, 'i4', 1)[0]
    ma_scheme = np.fromfile(fin, 'i4', 1)[0]
    
    print "-"*20
    print "n_cell    = ",n_cell
    print "box_size  = ",box_size
    print "n_grids   = ",n_grids
    print "ma_scheme = ", ma_scheme
    print "-"*20
    
    return n_cell, box_size, n_grids, ma_scheme


def read_grids(fname):
    """Read in a density grid from file `fname`."""
    
    grid = {}
    
    with open(fname, "rb") as fin:
        n_cell, box_size, n_grids, ma_scheme = read_grid_header(fin)
        data_desc = '('+str(n_cell[0])+','+str(n_cell[1])+','+str(n_cell[2])+')f4'
        for i_grid in xrange(n_grids):
            ident = np.fromfile(fin, 'S32', 1)[0]
            print("Reading grid "+ident)
            grid[ident] = np.fromfile(fin, data_desc, 1)[0]
    
    return grid, n_cell, box_size


def animate(index):
    """This function is called for each new frame and updates the data being
    plotted."""

    # create a new slice
    s = list(np.s_[:,:,:])
    s[axis] = index

    # update the data in the density image and quivers
    cax.set_data(grid['rho_r_dark'][s])
    quiver.set_UVC(grid['v_'+labels[(axis-1)%3]+'_r_dark'][s],
                   grid['v_'+labels[(axis+1)%3]+'_r_dark'][s])

    # update the figure title
    slice_pos = index*(cell_half_width[axis]*2.)+cell_half_width[axis]
    plt.title(labels[axis]+" = {:03.2f}".format(slice_pos)+r"h$^{-1}$ Mpc")

    # redraw
    fig.canvas.draw()
    
    # update the progress bar
    progressbar.update()



MOVIE = True
GRID_FNAME = "/Users/smutch/Work/data/grids/GiggleZ_LR/snapshot_441_dark_grid.dat"
FPS = 5
labels = ['x','y','z']
axis = 2
index = 0

# read in the grids
grid, n_cell, box_size = read_grids(GRID_FNAME)

# We will use a log10 normalisation for the density field so multiply the
# densities by 1e10 to get h^-1 Msol units.
grid['rho_r_dark'] = grid['rho_r_dark']*1.e10
min_density, max_density = grid['rho_r_dark'][grid['rho_r_dark']>0].min(), grid['rho_r_dark'].max()

# set up the figure
fig = plt.figure(0)
ax = plt.subplot(111)

# calculcate useful quantities
cell_half_width = box_size/n_cell/2.
extent = (cell_half_width[0],box_size[0]-cell_half_width[0],cell_half_width[1],box_size[1]-cell_half_width[1])
sliced_grid = {}

# create the requested slice
s = list(np.s_[:,:,:])
s[axis] = index

# plot the density grid
cax = plt.imshow(grid['rho_r_dark'][s], origin='lower', extent=extent,
                 interpolation='bicubic',
                 vmin=min_density, vmax=max_density,
                 norm=LogNorm())

if not MOVIE:
    contours = plt.contour(grid['rho_r_dark'][s], origin='lower', linewidths=1, extent=extent, colors='w', alpha=0.3)

# ass a colorbar
cb = plt.colorbar(cax)

# add quivers for the velocity field
X,Y = np.meshgrid(np.linspace(extent[0]+cell_half_width[0],extent[1]-cell_half_width[0],n_cell[0]),
                np.linspace(extent[2]+cell_half_width[1],extent[3]-cell_half_width[1],n_cell[2]))
quiver = ax.quiver(X,Y,grid['v_'+labels[(axis-1)%3]+'_r_dark'][s],grid['v_'+labels[(axis+1)%3]+'_r_dark'][s],
           units='xy', scale=200,
           color='w', alpha=0.7,
           pivot='tail',
           headaxislength=5,
           headwidth=2.5)

# add a key for the quivers that indicates the normailsation of the lengths
quiverkey = ax.quiverkey(quiver, 0.92, 0.93, 500, r'500 ks$^{-1}$', fontproperties={'weight': 'bold'}, labelcolor='w', alpha=1, coordinates='axes')

# set the axis labels and figure title
units = r"h$^{-1}$ [Mpc]"
plt.xlabel(labels[(axis+1)%3]+units)
plt.ylabel(labels[(axis-1)%3]+units)
cb.set_label(r"h$^{-1}$ [M$_{\odot}$ Mpc$^{-3}$]")
slice_pos = index*(cell_half_width[axis]*2.)+cell_half_width[axis]
plt.title(labels[axis]+" = {:03.2f}".format(slice_pos)+r"h$^{-1}$ Mpc")


# save the image or create a movie by recursively calling the `animate` function
if not MOVIE:
    plt.savefig("slice_"+labels[axis]+"{:04d}.png".format(index))
else:
    print "Generating movie..."
    progressbar = ProgressBar(n_cell[axis]+1)
    anim = animation.FuncAnimation(fig, animate, frames=n_cell[axis], blit=True)
    # Save as mp4. This requires ffmpeg to be installed.
    anim.save('test.mp4', fps=FPS, bitrate=1800*2)
    progressbar.__exit__(None,None,None)

