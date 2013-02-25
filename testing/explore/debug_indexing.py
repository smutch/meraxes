# -*- coding: utf-8 -*-
# imports
import numpy as np
from matplotlib import use as mpl_use
mpl_use('Agg')
from matplotlib import pyplot as plt
from glob import glob
import h5py as h5
# from astropy.utils.console import ProgressBar
from multiprocessing import Pool
from functools import partial

def read_snap(fname, props=None):
    
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

def gen_plotset(init_id, step=[]):

    print init_id

    snapshot = []
    halo = []
    id = init_id

    for s in xrange(n_steps-1,-1,-1):
        boolarr = (step[s]['id']==id)
        if any(boolarr==True):
            id = step[s]['desc_id'][boolarr]
            halo.append(step[s][boolarr][0])
            snapshot.append(int(snap_files[s][14:-5]))

    snapshot = np.array(snapshot)
    halo = np.array(halo)

    # with h5.File('tree.hdf5', 'w') as fout:
    #     fout['id'] = init_id
    #     fout.create_dataset('tree', data=halo)
    #     fout.create_dataset('snapshot', data=snapshot)

    fig = plt.figure(0, figsize=(12,12))

    type0 = np.ma.masked_array(halo, mask=(halo['type']==0))
    type0_snap = np.ma.masked_array(snapshot, mask=(halo['type']==0))
    type1 = np.ma.masked_array(halo, mask=(halo['type']==1))
    type1_snap = np.ma.masked_array(snapshot, mask=(halo['type']==1))
    
    ax = plt.subplot(321)
    ax.semilogy(snapshot, halo['V_max'], color='0.5')
    ax.semilogy(type0_snap, type0['V_max'])
    ax.semilogy(type1_snap, type1['V_max'])
    ax.set_ylabel('V_max')
    ax.set_xlim((30,116))
    ax.set_ylim((100,5000))
    
    ax = plt.subplot(322)
    ax.semilogy(snapshot, halo['M_vir'], color='0.5')
    ax.semilogy(type0_snap, type0['M_vir'])
    ax.semilogy(type1_snap, type1['M_vir'])
    ax.set_ylabel('M_vir')
    ax.set_xlim((30,116))
    ax.set_ylim((1e10, 1e15))
    
    ax = plt.subplot(323)
    ax.semilogy(snapshot, halo['R_halo'], label='R_halo')
    ax.semilogy(snapshot, halo['R_vir'], label='R_vir')
    ax.set_ylabel('R')
    leg = ax.legend(loc='upper left')
    plt.setp(leg.get_texts(), size='small')
    ax.set_xlim((30,116))
    ax.set_ylim((0.01,100))
    
    ax = plt.subplot(324)
    spin = np.array([np.linalg.norm(v) for v in halo['spin']])
    spin0 = np.ma.masked_array(spin, mask=(halo['type']==0))
    spin1 = np.ma.masked_array(spin, mask=(halo['type']==1))
    ax.semilogy(snapshot, spin, color='0.5')
    ax.semilogy(type0_snap, spin0)
    ax.semilogy(type1_snap, spin1)
    ax.set_ylabel('spin magnitude')
    ax.set_xlim((30,116))
    ax.set_ylim((1,1e3))
    
    #Fix the positions
    initial = halo['position_MBP'][0,:]
    final = halo['position_MBP'][-1,:]
    diff = halo['position_MBP']-final
    wbool = diff<-(box_size/2.)
    halo['position_MBP'][wbool]+=box_size
    wbool = diff>(box_size/2.)
    halo['position_MBP'][wbool]-=box_size
    
    ax = plt.subplot(325)
    ax.plot(halo['position_MBP'][:,0], halo['position_MBP'][:,1], color='0.5')
    ax.plot(type0['position_MBP'][:,0], type0['position_MBP'][:,1])
    ax.plot(type1['position_MBP'][:,0], type1['position_MBP'][:,1])
    ax.scatter(halo['position_MBP'][0,0], halo['position_MBP'][0,1], marker='o')
    ax.set_xlabel('x (MBP)')
    ax.set_ylabel('y (MBP)')
    ax.grid(True)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=90)
    # ax.set_xlim((-125,125))
    # ax.set_ylim((-125,125))
    
    ax = plt.subplot(326)
    ax.plot(halo['position_MBP'][:,0], halo['position_MBP'][:,2], color='0.5')
    ax.plot(type0['position_MBP'][:,0], type0['position_MBP'][:,2])
    ax.plot(type1['position_MBP'][:,0], type1['position_MBP'][:,2])
    ax.scatter(halo['position_MBP'][0,0], halo['position_MBP'][0,2], marker='o')
    ax.set_xlabel('x (MBP)')
    ax.set_ylabel('z (MBP)')
    ax.grid(True)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=90)
    # ax.set_xlim((-125,125))
    # ax.set_ylim((-125,125))

    plt.text(0.5, 0.95, "last M_vir = %.2e; id = %04d"%(z0_M_vir[init_id], init_id), 
         horizontalalignment='center',
         verticalalignment='center',
         transform=fig.transFigure)

    fname = 'plots/Mvir_%.3e.png'%z0_M_vir[init_id]
    plt.savefig(fname)
    del(fig)


if __name__ == '__main__':
    
    # read in the snapshots
    box_size = 125.
    snap_files = sorted(glob('../output/snap*.hdf5'))[::-1]
    n_steps = len(snap_files)

    print "Reading snapshots..."
    step = [read_snap(snap) for snap in snap_files]

    #loop through each halo, find it's z=0 descendant and record it's M_vir value
    print "Sorting by final M_vir value..."
    id_max =0
    for s in xrange(0, n_steps):
        local_max = np.max(step[s]['id'])
        if local_max>id_max:
            id_max =local_max

    z0_M_vir = np.zeros(id_max+1)
    for id_init in xrange(id_max+1):
        id = id_init
        for s in xrange(n_steps-1,-1,-1):
            boolarr = (step[s]['id']==id)
            if any(boolarr==True):
                id = step[s]['desc_id'][boolarr]
                last_mvir = step[s]['M_vir'][boolarr] 
        z0_M_vir[id_init] = last_mvir[0]
        
    #do a sort on the ids by z0_M_vir mass
    sorted_ids = np.argsort(z0_M_vir)[::-1]

    mapfunc = partial(gen_plotset, step=step)    
    worker_pool = Pool(6)
    worker_pool.map(mapfunc, sorted_ids)
