import sys
import os
import glob
import numpy as np
from ssimpl import meraxes, munge

def analysis(meraxes_fname):

    h = 0.702
    snaps = [99, 77, 62]
    sfr_limit = np.array([-0.44, -0.265, -0.18])
    sfr_limit *= (0.7**2)/(h**2)

    for snap, lim in zip(snaps, sfr_limit):
        # SMF
        fname_out = os.path.join(os.path.dirname(meraxes_fname),
                                 "smf_{:d}.npy".format(snap))

        gals, sim_props = meraxes.io.read_gals(meraxes_fname, h=h,
                                               snapshot=snap, sim_props=True)
        volume = sim_props["Volume"]

        edges = np.linspace(6, 11.5, 51)
        gals = np.compress(gals["StellarMass"] > 0, gals)
        smf = munge.mass_function(np.log10(gals["StellarMass"]*1.e10), volume,
                                  bins=edges)

        fout = open(fname_out, "ab")
        np.save(fout, smf[:, 1])
        fout.close()

        # SFRF
        fname_out = os.path.join(os.path.dirname(meraxes_fname),
                                 "sfrf_{:d}.npy".format(snap))

        edges = np.linspace(-2, 4, 51)
        gals = np.compress((gals["Sfr"] > 0) & (gals["GhostFlag"] == 0), gals)
        sfrf = munge.mass_function(np.log10(gals["Sfr"]), volume, bins=edges)

        fout = open(fname_out, "ab")
        np.save(fout, sfrf[:, 1])
        fout.close()

        # SMF with limit
        fname_out = os.path.join(os.path.dirname(meraxes_fname),
                                 "smf_wlim_{:d}.npy".format(snap))

        gals, sim_props = meraxes.io.read_gals(meraxes_fname, h=h,
                                               snapshot=snap, sim_props=True)
        volume = sim_props["Volume"]

        edges = np.linspace(6, 11.5, 51)
        gals = np.compress(np.log10(gals["Sfr"]) > lim, gals)
        smf = munge.mass_function(np.log10(gals["StellarMass"]*1.e10), volume,
                                  bins=edges)

        fout = open(fname_out, "ab")
        np.save(fout, smf[:, 1])
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
