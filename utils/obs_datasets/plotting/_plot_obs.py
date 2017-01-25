#!/usr/bin/env python
from _plot_obsMagorrian import *
from _plot_obsGSMF import *
from _plot_obsQLF import *
from _plot_obsBHMF import *
from _plot_obsGLF import *

__author__ = ('Yuxiang Qin')                     
__email__ = ("Yuxiang.L.Qin@Gmail.com")          
__all__  = ('plot_obs', 'plot_obs_legend',)

def plot_obs(ax, name, z,markersize=10,legend=True,hubble_h=0.678,silent=False):
    '''
    Notes:
    GSMF      -> z:[ 7.0, 5.0, 4.0, 2.0, 1.75, 1.3, 0.95, 0.55, 0.0]
    GLF       -> z:[10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0]
    QLF       -> z:[ 6.0, 5.0, 4.0, 3.0, 2.0, 1.5, 1.3, 0.5]
    BHMF      -> z:[ 0.5, 0.0]
    Magorrian -> z:[0.0] !!hubble_h doesn't work for Magorrian!!
    '''
    _plot_obs = {'GSMF':      plot_obsGSMF,\
                 'GLF':       plot_obsGLF,\
                 'QLF':       plot_obsQLF,\
                 'BHMF':      plot_obsBHMF,\
                 'Magorrian': plot_obsMagorrian}
    try:
        (_plot_obs[name])(ax,z,markersize=markersize,legend=legend,hubble_h=hubble_h,silent=silent)
    except:
        print 'available data:', _plot_obs.keys()
        return -1

def plot_obs_legend(ax, name, markersize=10):
    _plot_obs_legend = {'GSMF':      plot_obsGSMF_legend,\
                        'GLF':       plot_obsGLF_legend,\
                        'QLF':       plot_obsQLF_legend,\
                        'BHMF':      plot_obsBHMF_legend,\
                        'Magorrian': plot_obsMagorrian_legend}
    try:
        (_plot_obs_legend[name])(ax, markersize=markersize)
    except:
        print 'available data:', _plot_obs_legend.keys()
        return -1
