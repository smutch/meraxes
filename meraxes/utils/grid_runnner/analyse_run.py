import sys
import os
import glob
import numpy as np
from ssimpl import meraxes, munge

def analysis(meraxes_fname):

    h = 0.702
    snaps = [99, 77, 62]

    for snap in snaps:
        fname_out = os.path.join(os.path.dirname(meraxes_fname),
                                 "smf_{:d}.npy".format(snap))

        gals, sim_props = meraxes.io.read_gals(meraxes_fname, h=h, snapshot=snap, sim_props=True)
        volume = sim_props["Volume"]

        edges = np.linspace(6, 11.5, 51)
        gals = np.compress(gals["StellarMass"] > 0, gals)
        smf = munge.mass_function(np.log10(gals["StellarMass"]*1.e10), volume,
                                  bins=edges)

        fout = open(fname_out, "ab")
        np.save(fout, smf[:,1])
        fout.close()


def delete_input(meraxes_fname):

    dirname, fname = os.path.split(meraxes_fname)
    fname_base = fname[:-5]
    glob_pattern = dirname+'/'+fname_base+'*.hdf5'
    files = glob.glob(glob_pattern)
    for f in files:
        os.remove(f)


if __name__ == '__main__':
    fname = sys.argv[1]
    analysis(fname)
    delete_input(fname)
