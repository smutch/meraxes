# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""Generate a series of diagnostic plots for the horizontal merger trees.

Note - This code is rather brute force and certainly will not scale well with
high temporal or spatial resolution runs...
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from glob import glob
import h5py as h5
from astropy.utils.console import ProgressBar
from multiprocessing import Pool
from functools import partial

__author__ = "Simon Mutch"
__date__   = "26 Feb 2013"

# store the default colors being used
colors = plt.rcParams['axes.color_cycle']

# bitwise flags for the merger trees
tree_flags = {
    # 'simple'                        : 1,
    'strayed'                       : 2,
    'sputtered'                     : 4,
    'dropped'                       : 8,
    'merger'                        : 16,
    'bridged'                       : 32,
    'emerged'                       : 64,
    # 'bridge_progenitor'             : 128,
    # 'bridge_progenitor_unprocessed' : 256,
    # 'bridge_finalize'               : 512,
    # 'bridge_default'                : 1024,
    # 'found'                         : 2048,
    # 'main_progenitor'               : 4096,
    # 'unprocessed'                   : 8192,
    # 'invalid'                       : 16384,
}

def get_flags(value):
    """Given a tree_flags value, work out what actual flags are set.
    
    Args:
        value  -  the tree_flags value

    Returns:
        flag_list  -  a list of strings with the set flag names
    """

    flag_list = []

    for name, flag in tree_flags.iteritems():
        if (value & flag) == flag:
            flag_list.append(name)

    return flag_list


def check_main_progenitor_flag(value):
    """Check if the main_progenitor flag is set.

    Args:
        value  -  the tree_flags value

    Return:
        bool representing whether or not the main_progenitor flag is set
    """

    if (value & 4096)==4096:
        return True
    else:
        return False


def plot_events(ax, halo, snapshot):
    """Plot events marked by tree_flags value.
    
    Args:
        ax       - the plot axis
        halo     - halo history array
        snapshot - array of snapshots corresponding to each entry of `halo`
    """

    trans = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    flags = []
    for s in xrange(snapshot.size):
        halo_flags = get_flags(halo[s]['tree_flags'])
        if len(halo_flags)>0:
            flags.extend([(f, snapshot[s]) for f in halo_flags])

    for f in flags:
        ax.axvline(f[1], ls='--', color=colors[2])
        t = ax.text(f[1], 1.02, f[0], horizontalalignment='center',
                    verticalalignment='bottom', size='xx-small', rotation=90,
                    transform=trans)


def read_snap(fname, props=None):
    """Read in a snapshot from a simple hdf5 file.

    Args:
        fname  -  input filename
        props  -  requested properties (=None for all)

    Returns:
        halo  -  the halos from the input file
    """
    
    fin = h5.File(fname, 'r')
    
    n_halos = fin['id'].shape[0]
    
    # Generate the datatype
    halo_desc = []
    if props==None:
        props = fin.keys()
    for p in props:
        if len(fin[p].shape)==2:
            descr = (str(p), (fin[p].dtype.descr[0][1], fin[p].shape[1]))
        elif len(fin[p].shape)==1:
            descr = (str(p), fin[p].dtype.descr[0][1])
        halo_desc.append(descr)
    halo_dtype = np.dtype(halo_desc)
    
    # Allocate the memory
    halo = np.empty(n_halos, dtype=halo_dtype)
    
    # Read in the halo properties
    for p in props:
        halo[p][:] = fin[p].value
        
    # close the file!
    fin.close()
        
    return halo

def gen_plotset(init_id, step=[], last_M_vir=[], dead_halos=[]):
    """Generate a plot for a given halo ID."""

    snapshot = []
    halo = []
    id = init_id

    # trace the history of the halo with id=init_id
    for s in xrange(n_steps-1,-1,-1):
        boolarr = (step[s]['id']==id)
        if any(boolarr==True):
            halo.append(step[s][boolarr][0])
            snapshot.append(int(snap_files[s][14:-5]))
            if any(dead_halos[s][boolarr]):
                break
            if not check_main_progenitor_flag(step[s][boolarr][0]['tree_flags']):
                sys.stderr.write("ID=%d, step=%d :: Not the main progenitor by score...\n" % (id, s))
            id = step[s]['desc_id'][boolarr]
    snapshot = np.array(snapshot)
    halo = np.array(halo)

    # select/create the figure and clear it
    fig = plt.figure(0, figsize=(12,12))
    fig.clf()

    # work out the masks to select only type=0 and type=1 halos
    type0 = np.ma.masked_array(halo, mask=(halo['type']!=0))
    type0_snap = np.ma.masked_array(snapshot, mask=(halo['type']!=0))
    type1 = np.ma.masked_array(halo, mask=(halo['type']!=1))
    type1_snap = np.ma.masked_array(snapshot, mask=(halo['type']!=1))
    
    # V_max
    ax = plt.subplot(321)
    ax.set_xlim((30,116))
    ax.set_ylim((100,5000))
    ax.semilogy(snapshot, halo['V_max'], color='0.5')
    ax.semilogy(type0_snap, type0['V_max'], label='type 0')
    ax.semilogy(type1_snap, type1['V_max'], label='type 1')
    plot_events(ax, halo, snapshot)
    ax.set_ylabel('V_max')
    leg = ax.legend(loc='upper left', labelspacing=0.2)
    plt.setp(leg.get_texts(), size='x-small')
    
    # M_vir
    ax = plt.subplot(322)
    ax.semilogy(snapshot, halo['M_vir'], color='0.5')
    ax.semilogy(type0_snap, type0['M_vir'], label='type 0')
    ax.semilogy(type1_snap, type1['M_vir'], label='type 1')
    plot_events(ax, halo, snapshot)
    ax.set_ylabel('M_vir')
    leg = ax.legend(loc='upper left', labelspacing=0.2)
    plt.setp(leg.get_texts(), size='x-small')
    ax.set_xlim((30,116))
    ax.set_ylim((1e10, 1e15))
    
    # R_halo and R_vir
    ax = plt.subplot(323)
    ax.semilogy(snapshot, halo['R_halo'], label='R_halo')
    ax.semilogy(snapshot, halo['R_vir'], label='R_vir')
    plot_events(ax, halo, snapshot)
    ax.set_ylabel('R')
    leg = ax.legend(loc='upper left', labelspacing=0.2)
    plt.setp(leg.get_texts(), size='x-small')
    ax.set_xlim((30,116))
    ax.set_ylim((0.01,100))
    
    # spin (specific angular momentum) magnitude
    ax = plt.subplot(324)
    spin = np.array([np.linalg.norm(v) for v in halo['spin']])
    spin0 = np.ma.masked_array(spin, mask=(halo['type']!=0))
    spin1 = np.ma.masked_array(spin, mask=(halo['type']!=1))
    ax.semilogy(snapshot, spin, color='0.5')
    ax.semilogy(type0_snap, spin0, label='type 0')
    ax.semilogy(type1_snap, spin1, label='type 1')
    plot_events(ax, halo, snapshot)
    ax.set_ylabel('spin magnitude')
    leg = ax.legend(loc='upper left', labelspacing=0.2)
    plt.setp(leg.get_texts(), size='x-small')
    ax.set_xlim((30,116))
    ax.set_ylim((1,1e3))
    
    # fix the positions for periodic boundaries
    initial = halo['position_MBP'][0,:]
    final = halo['position_MBP'][-1,:]
    diff = halo['position_MBP']-final
    wbool = diff<-(box_size/2.)
    halo['position_MBP'][wbool]+=box_size
    wbool = diff>(box_size/2.)
    halo['position_MBP'][wbool]-=box_size
    
    # positions
    ax = plt.subplot(325)
    ax.plot(halo['position_MBP'][:,0], halo['position_MBP'][:,1], color='0.5')
    l, = ax.plot(type0['position_MBP'][:,0], type0['position_MBP'][:,1])
    ax.plot(type1['position_MBP'][:,0], type1['position_MBP'][:,1])
    ax.scatter(halo['position_MBP'][0,0], halo['position_MBP'][0,1], marker='o', s=60, c=l.get_color())
    ax.set_xlabel('x (MBP)')
    ax.set_ylabel('y (MBP)')
    ax.grid(True)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=90)
    
    # positions
    ax = plt.subplot(326)
    ax.plot(halo['position_MBP'][:,0], halo['position_MBP'][:,2], color='0.5')
    ax.plot(type0['position_MBP'][:,0], type0['position_MBP'][:,2])
    ax.plot(type1['position_MBP'][:,0], type1['position_MBP'][:,2])
    ax.scatter(halo['position_MBP'][0,0], halo['position_MBP'][0,2], marker='o', s=60, c=l.get_color())
    ax.set_xlabel('x (MBP)')
    ax.set_ylabel('z (MBP)')
    ax.grid(True)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=90)

    # label the plot
    plt.text(0.5, 0.95, "last M_vir = %.2e; id = %04d"%(last_M_vir[init_id], init_id), 
         horizontalalignment='center',
         verticalalignment='center',
         transform=fig.transFigure)

    # save
    fname = 'plots/Mvir_%.3e.png'%last_M_vir[init_id]
    plt.savefig(fname)


if __name__ == '__main__':
    
    # read in the snapshots
    box_size = 125.
    snap_files = sorted(glob('../output/snap*.hdf5'))[::-1]
    n_steps = len(snap_files)

    print "Reading snapshots..."
    step = [read_snap(snap) for snap in snap_files]

    # find out how many unique ids there are
    id_max =0
    for s in xrange(0, n_steps):
        local_max = np.max(step[s]['id'])
        if local_max>id_max:
            id_max =local_max

    
    # identify dead halos
    print "Identifying 'dead' halos..."
    dead_halos = []
    with ProgressBar(n_steps) as bar:
        for s in xrange(0,n_steps):
            dead_halos.append(np.zeros(step[s].shape[0], np.bool))
            for i in xrange(step[s].shape[0]):
                if not dead_halos[s][i]:
                    family_members = np.where(step[s]['desc_id']==step[s][i]['desc_id'])[0]
                    if family_members.size>1:
                        sorted_ind = np.argsort(step[s][family_members], order='M_vir')[::-1]
                        dead_halos[s][family_members[sorted_ind[0]]] = False
                        dead_halos[s][family_members[sorted_ind[1:]]] = True
            bar.update()


    # loop through each halo, find it's final descendant and record it's M_vir value
    print "Finding final halo masses..."
    last_M_vir = np.zeros(id_max+1)
    update_every = int((id_max+1)/100)
    with ProgressBar(100) as bar:
        for id_init in xrange(id_max+1):
            id = id_init
            for s in xrange(n_steps-1,-1,-1):
                boolarr = (step[s]['id']==id)
                if any(boolarr):
                    id = step[s]['desc_id'][boolarr]
                    last_mass = step[s]['M_vir'][boolarr] 
                if any(dead_halos[s][boolarr]):
                    break
            last_M_vir[id_init] = last_mass[0]
            if ((id_init+1)%update_every)==0:
                bar.update()

    # do a sort on the ids by last_M_vir mass
    print "Sorting by final M_vir value..."
    sorted_ids = np.argsort(last_M_vir)[::-1]

    # use python's built in multiprocessing module to generate the plots
    print "Generating plots..."
    mapfunc = partial(gen_plotset, step=step, last_M_vir=last_M_vir, dead_halos=dead_halos)    
    worker_pool = Pool(6)
    worker_pool.map(mapfunc, sorted_ids)
