#!/usr/bin/env python

"""Plot the merger tree of galaxy with a given ID.

Usage: 
    merger_tree.py <meraxes_file> <galaxy_ID> <last_snapnum> [--output=<dir_path> --format=<ext>]

Options:
    --output=<dir_path>   Target directory for output figures [default: ./plots]
    --format=<ext>        Image format (e.g. png, pdf) [default: png]

"""

import sys
import os
from collections import OrderedDict
from random import sample

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docopt import docopt

__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.join(__script_dir__, "../utils"))
from samtools import io as samio

__author__ = "Simon Mutch"
__date__   = "2013/08/14"


class Graph:
    
    def __init__(self, n_snapshots):
        self.nodes = OrderedDict()
        self.edges = OrderedDict()
        self.snap_count = np.zeros(n_snapshots, int)
        self.type_first = "FirstProgenitor"
        self.type_next = "NextProgenitor"
        self.__counter__ = 0
        
    def add_node(self, index, snapshot, galaxy=None):
        name = str(self.__counter__)
        if galaxy!=None:
            self.nodes[name] = {"index":index, "snapshot":snapshot, "type":self.type_first, "galaxy":galaxy.copy()}
        else:
            self.nodes[name] = {"index":index, "snapshot":snapshot, "type":self.type_first, "galaxy":None}
        self.__counter__+=1
        self.snap_count[snapshot]+=1
        return name
        
    def add_edge(self, start, end, edge_type):
        self.edges["{:s}->{:s}".format(start, end)] = {"start" : start, "end" : end, "type" : edge_type}
        self.nodes[end]["type"]=edge_type
        
    def iter_nodes(self):
        return self.nodes.itervalues()
    
    def iter_edges(self):
        return self.edges.itervalues()
    
    def calc_positions(self):
        pos = np.zeros((self.__counter__, 2), float)
        x_i = 0
        for i, node in enumerate(self.iter_nodes()):
            pos[i,1] = node["snapshot"]
            if node["type"]=="NextProgenitor":
                x_i+=1
            pos[i,0] = x_i
        for i, node in enumerate(self.iter_nodes()):
            pos[i,0] =  pos[i,0]/x_i*100.0
            node["plot_pos"] = pos[i]
            
        # now loop back through by following edges
        return pos


def walk(snap, ind, fp_ind, np_ind):
    if ind>-1:
        node = G.add_node(ind, snap)
    
    first_ind = fp_ind[snap][ind]
    if first_ind>-1:
        first_node = walk(snap-1, first_ind, fp_ind, np_ind)
        G.add_edge(node, first_node, G.type_first)
        next_ind = np_ind[snap-1][first_ind]
        while next_ind>-1:
            next_node = walk(snap-1, next_ind, fp_ind, np_ind)
            G.add_edge(node, next_node, G.type_next)
            next_ind = np_ind[snap-1][next_ind]
    
    return node


if __name__ == '__main__':
    args = docopt(__doc__)
    fname_gals = args["<meraxes_file>"]
    galaxy_ID = int(args["<galaxy_ID>"])
    last_snapnum = int(args["<last_snapnum>"])
    fig_format = args["--format"]
    output_dir = args["--output"]
    
    # Check that we can construct the tree
    snaplist, zlist, _ = samio.read_snaplist(fname_gals)
    if last_snapnum not in snaplist:
        raise ValueError("last_snapnum not present in meraxes_file")
    if (np.diff(snaplist)!=1).any():
        raise ValueError("There are non-consecutive snapshots in meraxes_file")
    snaplist = np.arange(snaplist[0], last_snapnum+1, dtype=int)

    # Read in the walk indices 
    fp_ind = []
    np_ind = []
    for snap in snaplist:
        try:
            fp_ind.append(samio.read_firstprogenitor_indices(fname_gals, snap))
        except:
            fp_ind.append([])
        try:
            np_ind.append(samio.read_nextprogenitor_indices(fname_gals, snap))
        except:
            np_ind.append([])
 
    # Find the index of our requested galaxy
    gal = samio.read_gals(fname_gals, snapshot=last_snapnum, verbose=False)
    ind = np.argwhere(gal["ID"]==galaxy_ID)[0][0]

    # Walk the fp and np indices and construct the graph
    G = Graph(snaplist.size)
    first_node = walk(last_snapnum, ind, fp_ind, np_ind)

    # Attach the galaxies to the graph
    for snap in snaplist:
        gal = samio.read_gals(fname_gals, snapshot=snap, verbose=False)
        for node in G.iter_nodes():
            if node["snapshot"] == snap:
                node["galaxy"] = gal[node["index"]]

    # Plot
    colors = plt.rcParams['axes.color_cycle']
    fig, ax = plt.subplots(1,1, dpi=plt.rcParams['figure.dpi']*2,
                           figsize=np.array(plt.rcParams['figure.figsize'])*2,
                           facecolor=None)
    ax.set_ylabel("Snapshot")
    plt.rcParams['font.size']=18

    # start by drawing the edges
    pos = G.calc_positions()
    for e in G.iter_edges():
        start_pos = G.nodes[e["start"]]["plot_pos"]
        end_pos   = G.nodes[e["end"]]["plot_pos"]
        ax.plot([start_pos[0], end_pos[0]], [start_pos[1], end_pos[1]], 'k-', lw=2, alpha=0.3, zorder=0)

    # now draw the nodes
    mstar = np.array([n["galaxy"]["StellarMass"] for n in G.iter_nodes()])
    mstar = np.log10(mstar*1.e10)*10
    sel = (~np.isinf(mstar))
    mstar = mstar[sel]**5
    mstar = (mstar/mstar.max())*300

    mvir = np.array([n["galaxy"]["Mvir"] for n in G.iter_nodes()])
    mvir = np.log10(mvir[sel]*1.e10)
    mvir_max = mvir.max()
    mvir_min = mvir.min()

    pos = pos[sel]

    ghost_flag = np.array([n["galaxy"]["GhostFlag"] for n in G.iter_nodes()])
    ghost_sel = (ghost_flag[sel]==1) 
    gal_type = np.array([n["galaxy"]["Type"] for n in G.iter_nodes()])

    type1_sel = (gal_type[sel]==1) & ~ghost_sel
    type2_sel = (gal_type[sel]==2) & ~ghost_sel
    type0_sel = (~type1_sel & ~type2_sel)

    cmap = plt.cm.winter
    sc = ax.scatter(pos[type0_sel,0], pos[type0_sel,1], vmin=mvir_min,
                    vmax=mvir_max, marker='o', c=mvir[type0_sel],
                    s=mstar[type0_sel], cmap=cmap, label="Type 0")
    ax.scatter(pos[type1_sel,0], pos[type1_sel,1], vmin=mvir_min,
               vmax=mvir_max, marker='^', c=mvir[type1_sel],
               s=mstar[type1_sel], cmap=cmap, label="Type 1")
    ax.scatter(pos[type2_sel,0], pos[type2_sel,1], vmin=mvir_min,
               vmax=mvir_max, marker='v', c=mvir[type2_sel],
               s=mstar[type2_sel], cmap=cmap, label="Type 2")
    ax.scatter(pos[ghost_sel,0], pos[ghost_sel,1], vmin=mvir_min,
               vmax=mvir_max, marker='s', facecolor='0.5', s=mstar[ghost_sel],
               label="Ghosts")

    # add labels for each branch ID
    for n in G.iter_nodes():
        if (n["type"]==G.type_next) or (n==G.nodes[first_node]):
            plot_pos = n["plot_pos"]
            ax.text(plot_pos[0], plot_pos[1]+1,
                    "{:d}".format(n["galaxy"]["ID"]),
                    horizontalalignment="center",
                    verticalalignment="bottom",
                    size="x-small",
                    rotation='vertical',
                    color='0.5')

    # add a color bar
    # import IPython; IPython.embed()
    cb = fig.colorbar(sc)
    cb.set_clim(mvir.min(), mvir.max())
    cb.set_label(r"$\log_{10}(M_{\rm vir}/{\rm M_{\odot}})$")
    # labels = cb.ax.get_yticklabels()
    # for i, t in enumerate(labels):
    #         labels[i] = "{:.1f}".format(float(t.get_text())/5.)
    #     cb.ax.set_yticklabels(labels)

    # add a legend
    leg = ax.legend(loc='upper right', ncol=4)
    for t in leg.get_texts():
        t.set_size("medium")

    # tidy up
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_xticklines(), visible=False)
    plt.setp(ax.get_yticklines(), visible=False)
    ax.grid(True, axis='y', color='0.9')
    ax.grid(False, axis='x')
    ax.set_ylim((np.argwhere(G.snap_count==0)[-1], G.snap_count.size+1))
    ax.set_xlim((-5, 101))
    ax.set_frame_on(False)

    # TODO: Add the redshift axis
    # # add the redshift y-axis
    # axz = ax.twinx()
    # axz.set_ylim(ax.get_ylim())
    # axz_labels = [l.get_text() for l in axz.get_yticklabels()]
    # snaplist_max = snaplist.max()
    # snaplist_min = snaplist.min()
    # # for i, t in enumerate(axz_labels):
    # #     if len(t)>0:
    # #         if (int(t)<=snaplist_max) & (int(t)>=snaplist_min):
    # #             axz_labels = "{:.1f}".format(zlist[int(t)])
    # #         else:
    # #             axz_labels = ""
    # # axz.set_yticklabels(axz_labels)

    # save
    plt.tight_layout()
    plt.savefig("{:s}/merger_tree_ID{:d}.{:s}".format(output_dir, galaxy_ID, fig_format))

