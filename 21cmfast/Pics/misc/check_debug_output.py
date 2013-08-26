from __future__ import division
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from glob import glob

BOXES_DIR = "/home/smutch/data/21cmfast_runs/GiggleZ/MR/Boxes"
LITTLE_H = 0.705
BOX_SIDE = 125./LITTLE_H  # Mpc/h

# define the slice that we will plot
slice_start = 10
slice_width = 5
s = np.s_[:,:,slice_start:slice_start+slice_width]

plt.rcParams['figure.figsize']=(16,12)

# get the present snapshot numbers and sort them
with h5.File(BOXES_DIR+"/halos.hdf5", "r") as fin:
    snaps = sorted([int(k[-4:]) for k in fin.keys()])

# get the xH files and similarly sort them to match up
xHfiles = sorted(glob(BOXES_DIR+"/xH*"))

if len(snaps) != len(xHfiles):
    raise ValueError("Length of snaps and xHfiles do not match!")


# read in the density grid, halos and xH grid for each snapshot then plot
for i_snap, snap in enumerate(snaps):

    # read in the data and reshape where necessary
    # Note that the read data should already be in h=0.705 units
    snap_group = 'snap{:04d}'.format(snap)
    with h5.File(BOXES_DIR+"/halos.hdf5", "r") as fin:
        halo_mass = fin[snap_group+"/mass"][:]
        halo_pos = fin[snap_group+"/pos"][:]
    with h5.File(BOXES_DIR+"/grids.hdf5", "r") as fin:
        density = fin[snap_group+"/density"][:]
        n_cells = density.size
        dim = int(np.round(n_cells**(1./3.)))
        density.shape = [dim,]*3
    xH = np.fromfile(xHfiles[i_snap], 'f4', n_cells)
    xH.shape = [dim,]*3

    # munge the data 
    zpos_range = (BOX_SIDE/float(dim)*slice_start,
                  BOX_SIDE/float(dim)*(slice_start+slice_width))
    halo_pos *= BOX_SIDE/float(dim)
    halo_sel = ((halo_pos[:,2]>=zpos_range[0]) & (halo_pos[:,2]<zpos_range[1]))
    halo_mass = np.log10(halo_mass[halo_sel])
    print halo_mass
    print "selected {:d} halos".format(halo_mass.size)
    density = density[s].sum(axis=2)
    halo_pos = halo_pos[halo_sel]
    xH = xH[s].sum(axis=2)

    # TODO: plot the results...
    fig, ax = plt.subplots(1,1)
    ax.imshow(density.T, origin='lower', extent=(0,BOX_SIDE,0,BOX_SIDE))
    ax.imshow(xH.T, origin='lower', alpha=0.5, cmap=plt.cm.Greys_r, extent=(0,BOX_SIDE,0,BOX_SIDE))
    ax.scatter(halo_pos[:,0], halo_pos[:,1], alpha=0.5, c=halo_mass,
               s=(10.**halo_mass)/1e10*5)
    ax.set_xlabel("x (Mpc)")
    ax.set_ylabel("y (Mpc)")
    ax.set_ylim((0, BOX_SIDE))
    ax.set_xlim((0, BOX_SIDE))
    plt.tight_layout()
    plt.savefig("./frames/snap{:04d}.png".format(snap))

