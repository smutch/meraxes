from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

dim = 256
edge = 125.0
volume = edge**3

# read in data
halos = {}
halos['109'] = np.loadtxt("./halos_mbp_109.txt")
halos['110'] = np.loadtxt("./halos_mbp_110.txt")

# munge
for k in halos.iterkeys():
    # sel = ((halos[k][:,3]>0) & (halos[k][:,3]<8))
    # halos[k] = halos[k][sel]
    # volume = volume/256*edge/float(dim)*7
    halos[k][:,0] = np.log10(halos[k][:,0])

# plot the mass functions
fig = plt.figure(0)
ax = plt.subplot(111)
for name, snap in halos.iteritems():
    phi, bin_edges = np.histogram(snap[:,0], bins=20)
    bin_width = bin_edges[1]-bin_edges[0]
    phi = np.log10(phi/(bin_width*(volume)))
    ax.step(bin_edges[:-1], phi, where='post', alpha=0.75, label=name) 

plt.legend(loc='lower left')
plt.savefig('./frames/mfunc.png')

