#!/usr/bin/env python

"""Convert the vertical tree binary files to hdf5 format."""

import numpy as np
import h5py as h5

__author__ = "Simon Mutch"
__date__   = "12 Mar 2013"

FNAME = "/home/gpoole/GiggleZ/GiggleZ_LR/trees/GiggleZ_LR_step_008_scan_016/vertical/subgroups/GiggleZ_LR_step_008_scan_016.subgroup_trees.0"

halo_desc = ([
  # merger tree pointers and match type
  ('descendant', 'i4'),
  ('progenitor_first', 'i4'),
  ('progenitor_next', 'i4'),
  ('group_halo_first', 'i4'),
  ('group_halo_next', 'i4'),
  ('tree_flags', 'i4'),

  # properties of halo
  ('n_particles', 'i4'),
  ('M_vir', 'f4'),  # This is the FoF mass for the most massive substructure
  ('R_vir', 'f4'),
  ('position_MBP', ('f4', 3)),
  ('velocity_COM', ('f4', 3)),
  ('sigma_v', 'f4'),
  ('V_max', 'f4'),
  ('spin', ('f4', 3)),
  ('most_bound_id', 'q'),

  # original position in halo-finder output
  ('snap_num', 'i4'),
  ('halo_index', 'i4'),
  ('id', 'i4'),
  ('group_id', 'i4'),
  ('desc_id', 'i4'),
])
halo_dtype = np.dtype({'names':[halo_desc[ii][0] for ii in
                                xrange(len(halo_desc))], 
                       'formats':[halo_desc[ii][1] for ii in
                                xrange(len(halo_desc))]},
                      align=True)

# read the input tree file
with open(FNAME, 'rb') as fin:
    # grab the header info
    n_trees = np.fromfile(fin, dtype='i4', count=1)[0]
    tot_n_halos = np.fromfile(fin, dtype='i4', count=1)[0]
    tree_n_halos = np.fromfile(fin, dtype='i4', count=n_trees)
    # read in all of the halos
    halos = np.fromfile(fin, dtype=halo_dtype, count=tot_n_halos)

# Work out the number of snapshots
snap_min = halos['snap_num'].min()
snap_max = halos['snap_num'].max()

# Generate the depth-first indices
dfi = np.zeros(tot_n_halos, 'i4')
ind=0
for tree in xrange(n_trees):
    n_halos = tree_n_halos[tree]
    dfi[ind:ind+n_halos] = range(n_halos)

def test():
    tot = 0
    for ii in xrange(117):
        n = len(np.where(halos['snap_num']==ii)[0])
        print "SNAPSHOT %d: %d halos" % (ii, n)
        tot+=n
    print "TOTAL: %d" % tot
    print "tot_n_halos = %d" % tot_n_halos
# test()

# create the hdf5 files
for ii in xrange(snap_min, snap_max+1):
    selection = halos['snap_num']==ii
    snap_halos = halos[selection]
    with h5.File('output/vertical/snap%03d.hdf5'%ii, 'w') as fout:
        for prop in halo_dtype.names:
            fout.create_dataset(prop, data=snap_halos[prop])
        # add the type field
        type = np.zeros(snap_halos.shape[0], 'i4')
        type[snap_halos['group_halo_first']!=dfi[selection]]=1
        fout.create_dataset('type', data=type)
    
