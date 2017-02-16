#!/usr/bin/env python
import numpy as np

__author__ = ('Yuxiang Qin')                     
__email__ = ("Yuxiang.L.Qin@Gmail.com")          
__all__  = ('plot_obsMagorrian','plot_obsMagorrian_legend')

####################################################################################
# Scott et al. 2013, 2015
def plot_scott13(x,xerrl,xerrr,y,yerrl,yerrr,ax,marker='>',label=None,color='gray',markersize=10):
    ax.errorbar([np.log10(x)+10.],[np.log10(y)+8.],\
                xerr=([xerrl/x/2.3026],[xerrr/x/2.3026]),\
                yerr=([yerrl/y/2.3026],[yerrr/y/2.3026]),\
                marker=marker,label=label,color=color,markersize=markersize)
    
# all can be found in Graham & Scott 2015
def plot_scott15(x,xerrl,xerrr,y,yerrl,yerrr,ax,marker='o',label=None,color='gray',markersize=10):
    ax.errorbar([np.log10(x)+9.],[np.log10(y)+5.],\
                xerr=([xerrl/x/2.3026],[xerrr/x/2.3026]),\
                yerr=([yerrl/y/2.3026],[yerrr/y/2.3026]),\
                marker=marker,label=label,color=color,markersize=markersize)
def plot_scott15_empty(x,xerrl,xerrr,y,yerrl,yerrr,ax,marker='o',label=None,color='black',markersize=10):
    ax.errorbar([np.log10(x)+9.],[np.log10(y)+5.],\
                xerr=([xerrl/x/2.3026],[xerrr/x/2.3026]),\
                yerr=([yerrl/y/2.3026],[yerrr/y/2.3026]),\
                fmt=marker,label=label,color=color,mfc="None",markersize=markersize)

    
def plot_jiang11(x,y,ax,marker='o',label=None,color='black',markersize=10):
    ax.errorbar([x],[y],[0],fmt=marker,label=label,color=color,mfc="None",markersize=markersize)
    
# Sanghvi et al. 2014
def plot_sanghvi14(x,y,ax,marker='*',label=None,color='gray',markersize=10):
    ax.errorbar([x],[y],[0],marker=marker,label=label,color=color,markersize=markersize)

####################################################################################
# Wyithe 2006
#def wyithe06(x):
#    alpha, beta1, beta2 = 8.05, 1.15, 0.12
#    return alpha+beta1*(x-11)+beta2*(x-11)**2.


def plot_obsMagorrian_z0pt0(ax,markersize=10,legend=True,hubble_h=None):
    if legend:
        plot_scott13(69,32,59,39,5,4,ax=ax,label='Scott et al. 2015; z:0',markersize=markersize)
    else:
        plot_scott13(69,32,59,39,5,4,ax=ax,markersize=markersize)
    plot_scott13(37,17,32,11,2,2,ax=ax,markersize=markersize)
    plot_scott13(1.4,0.8,2.0,0.45,0.01,0.17,ax=ax,markersize=markersize)
    plot_scott13(55,33,80,25,7,7,ax=ax,markersize=markersize)
    plot_scott13(27,12,23,24,10,10,ax=ax,markersize=markersize)
    plot_scott13(2.4,1.4,3.5,0.044,0.022,0.044,ax=ax,markersize=markersize)
    plot_scott13(0.46,0.28,0.68,1.4,0.3,0.9,ax=ax,markersize=markersize)
    plot_scott13(1.0,0.6,1.5,0.73,0,0,ax=ax,markersize=markersize)
    plot_scott13(19,9,16,9,0.8,0.9,ax=ax,markersize=markersize)
    plot_scott13(23,10,19,58,3.5,3.5,ax=ax,markersize=markersize)
    plot_scott13(0.61,0.36,0.89,0.1,0.05,0.1,ax=ax,markersize=markersize)
    plot_scott13(4.6,2.7,6.6,8.3,1.3,2.7,ax=ax,markersize=markersize)
    plot_scott13(11,5,9,0.39,0.09,0.26,ax=ax,markersize=markersize)
    plot_scott13(1.9,1.1,2.7,0.42,0.04,0.04,ax=ax,markersize=markersize)
    plot_scott13(4.5,2.7,6.6,0.084,0.003,0.003,ax=ax,markersize=markersize)
    plot_scott13(1.4,0.8,2.1,0.66,0.03,0.03,ax=ax,markersize=markersize)
    plot_scott13(0.66,0.4,0.97,0.73,0.35,0.69,ax=ax,markersize=markersize)
    plot_scott13(4.7,2.8,6.9,15,2,2,ax=ax,markersize=markersize)
    plot_scott13(26,12,22,4.7,0.6,0.6,ax=ax,markersize=markersize)
    plot_scott13(2.0,1.2,2.9,0.083,0.004,0.004,ax=ax,markersize=markersize)
    plot_scott13(0.39,0.23,0.57,0.14,0.13,0.02,ax=ax,markersize=markersize)
    plot_scott13(0.35,0.21,0.52,0.15,0.1,0.09,ax=ax,markersize=markersize)
    plot_scott13(0.30,0.18,0.45,0.40,0.05,0.04,ax=ax,markersize=markersize)
    plot_scott13(3.5,2.1,5.1,0.12,0.005,0.005,ax=ax,markersize=markersize)
    plot_scott13(1.9,1.1,2.4,1.7,0.2,0.2,ax=ax,markersize=markersize)
    plot_scott13(0.88,0.52,1.28,0.024,0.012,0.024,ax=ax,markersize=markersize)
    plot_scott13(1.9,1.1,2.7,8.8,2.7,10.0,ax=ax,markersize=markersize)
    plot_scott13(0.93,0.56,1.37,0.14,0.06,0.10,ax=ax,markersize=markersize)
    plot_scott13(1.24,0.7,1.8,2.0,0.5,0.5,ax=ax,markersize=markersize)
    plot_scott13(0.86,0.51,1.26,0.073,0.015,0.015,ax=ax,markersize=markersize)
    plot_scott13(2.0,0.9,1.7,0.77,0.06,0.04,ax=ax,markersize=markersize)
    plot_scott13(5.4,2.5,4.7,4.0,1.0,1.0,ax=ax,markersize=markersize)
    plot_scott13(1.2,0.7,1.7,0.17,0.02,0.01,ax=ax,markersize=markersize)
    plot_scott13(4.9,2.9,7.1,0.34,0.02,0.02,ax=ax,markersize=markersize)
    plot_scott13(2.0,1.2,2.9,0.34,0.02,0.02,ax=ax,markersize=markersize)
    plot_scott13(0.66,0.4,0.97,0.058,0.008,0.008,ax=ax,markersize=markersize)
    plot_scott13(5.1,3.0,7.4,3.1,0.6,1.4,ax=ax,markersize=markersize)
    plot_scott13(2.6,1.5,3.8,1.3,0.5,0.5,ax=ax,markersize=markersize)
    plot_scott13(3.2,1.5,2.7,2.0,0.6,1.1,ax=ax,markersize=markersize)
    plot_scott13(100,46,86,97,26,30,ax=ax,markersize=markersize)
    plot_scott13(1.4,0.9,2.1,8.1,1.9,2.0,ax=ax,markersize=markersize)
    plot_scott13(0.88,0.53,1.3,1.8,0.3,0.6,ax=ax,markersize=markersize)
    plot_scott13(1.3,0.8,1.9,0.65,0.07,0.07,ax=ax,markersize=markersize)
    plot_scott13(0.56,0.34,0.82,0.39,0.01,0.01,ax=ax,markersize=markersize)
    plot_scott13(29,13,25,5,1,1,ax=ax,markersize=markersize)
    plot_scott13(6.1,2.8,5.2,3.3,2.5,0.9,ax=ax,markersize=markersize)
    plot_scott13(0.65,0.39,0.96,4.5,1.5,2.3,ax=ax,markersize=markersize)
    plot_scott13(1,0.6,1.3,0.075,0.002,0.002,ax=ax,markersize=markersize)
    plot_scott13(2.0,1.2,3.0,0.68,0.13,0.13,ax=ax,markersize=markersize)
    plot_scott13(6.9,3.2,5.9,1.2,0.9,0.4,ax=ax,markersize=markersize)
    plot_scott13(1.4,0.6,1.2,0.13,0.08,0.08,ax=ax,markersize=markersize)
    plot_scott13(7.7,3.6,6.6,4.7,0.5,0.5,ax=ax,markersize=markersize)
    plot_scott13(0.9,0.5,1.3,0.59,0.09,0.03,ax=ax,markersize=markersize)
    plot_scott13(3.9,2.3,5.7,6.4,0.4,0.4,ax=ax,markersize=markersize)
    plot_scott13(1.8,1.1,2.7,0.79,0.33,0.38,ax=ax,markersize=markersize)
    plot_scott13(8.4,3.9,7.2,3.9,0.4,0.4,ax=ax,markersize=markersize)
    plot_scott13(27,12,23,47,10,10,ax=ax,markersize=markersize)
    plot_scott13(6.0,2.8,5.2,1.8,0.1,0.2,ax=ax,markersize=markersize)
    plot_scott13(0.43,0.26,0.64,0.06,0.014,0.014,ax=ax,markersize=markersize)
    plot_scott13(1.0,0.6,1.5,0.016,0.004,0.004,ax=ax,markersize=markersize)
    plot_scott13(122,57,105,210,160,160,ax=ax,markersize=markersize)
    plot_scott13(0.3,0.18,0.45,0.014,0.007,0.014,ax=ax,markersize=markersize)
    plot_scott13(29,13,25,7.4,3.0,4.7,ax=ax,markersize=markersize)
    plot_scott13(11,5,10,1.6,0.4,0.3,ax=ax,markersize=markersize)
    plot_scott13(20,9,17,6.8,0.7,0.7,ax=ax,markersize=markersize)
    plot_scott13(0.8,0.4,1.0,2.6,1.5,0.4,ax=ax,markersize=markersize)
    plot_scott13(24,11,20,11,1,1,ax=ax,markersize=markersize)
    plot_scott13(78,36,67,37,11,18,ax=ax,markersize=markersize)
    plot_scott13(96,44,83,5.9,2,2,ax=ax,markersize=markersize)
    plot_scott13(3.6,2.1,5.2,0.31,0.004,0.004,ax=ax,markersize=markersize)
    plot_scott13(2.6,1.5,3.8,0.1,0.001,0.001,ax=ax,markersize=markersize)
    plot_scott13(55,26,48,3.7,1.5,2.6,ax=ax,markersize=markersize)
    plot_scott13(1.4,0.8,2.0,0.55,0.19,0.26,ax=ax,markersize=markersize)
    plot_scott13(64,30,55,13,4,5,ax=ax,markersize=markersize)
    plot_scott13(1.2,0.7,1.8,0.11,0.005,0.005,ax=ax,markersize=markersize)
    plot_scott13(0.9,0.5,1.2,1.5,0.8,0.75,ax=ax,markersize=markersize)
    
    # Thornton et al. 2008
    if legend:
        plot_scott15(1.2,0,0,2.0,0.3,0.3,ax=ax,marker='^',label='Thornton et al. 2008; z:0',markersize=markersize)
    else:
        plot_scott15(1.2,0,0,2.0,0.3,0.3,ax=ax,marker='^',markersize=markersize)
    # Jiang et al. 2011
    if legend:
        plot_jiang11(9.25,5.7,ax=ax,marker='o',label='Jiang et al. 2011; z:0',markersize=markersize)
    else:
        plot_jiang11(9.25,5.7,ax=ax,marker='o',markersize=markersize)
    plot_jiang11(9.73,6.2,ax=ax,marker='o',markersize=markersize)
    plot_jiang11(10.09,5.8,ax=ax,marker='o',markersize=markersize)
    plot_jiang11(9.77,6.3,ax=ax,marker='o',markersize=markersize)
    plot_jiang11(9.92,5.9,ax=ax,marker='o',markersize=markersize)
    # Jiang et al. 2013
    if legend:
        plot_scott15(5.4,0,0,16,0,0,ax=ax,marker='o',label='Jiang et al. 2013; z:0',markersize=markersize)
    else:
        plot_scott15(5.4,0,0,16,0,0,ax=ax,marker='o',markersize=markersize)
    # Yuan et al. 2014
    if legend:
        plot_scott15(8.4,3.5,3.6,12.2,0,0,ax=ax,marker='D',label='Yuan et al. 2014; z:0',markersize=markersize)
    else:
        plot_scott15(8.4,3.5,3.6,12.2,0,0,ax=ax,marker='D',markersize=markersize)
    plot_scott15(14.7,6.2,6.1,5.85,0.75,0.75,ax=ax,marker='D',markersize=markersize)
    # Reines et al. 2013
    if legend:
        plot_scott15_empty(2.57,0,0,5.0,0,0,ax=ax,marker='s',label='Reines et al. 2013; z:0 BPT AGNs',markersize=markersize)
    else:
        plot_scott15_empty(2.57,0,0,5.0,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15_empty(2.29,0,0,2.5,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15_empty(1.32,0,0,0.8,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15_empty(2.95,0,0,12.6,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15_empty(1.26,0,0,1.0,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15_empty(2.88,0,0,1.6,0,0,ax=ax,marker='s',markersize=markersize)
    if legend:
        plot_scott15(2.57,0,0,2.5,0,0,ax=ax,marker='s',label='Reines et al. 2013; z:0 BPT composites',markersize=markersize)
    else:
        plot_scott15(2.57,0,0,2.5,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15(2.14,0,0,5.0,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15(1.32,0,0,1.3,0,0,ax=ax,marker='s',markersize=markersize)
    plot_scott15(1.74,0,0,1.6,0,0,ax=ax,marker='s',markersize=markersize)
    # Mathur et al. 2012
    if legend:
        plot_scott15(24.2,10.4,10.4,71,0,0,ax=ax,marker='v',label='Mathur et al. 2012; z:0',markersize=markersize)
    else:
        plot_scott15(24.2,10.4,10.4,71,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(16.8,7.2,7.1,210,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(18.4,7.9,7.8,54,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(26.6,11.4,11.3,45,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(66.7,28.6,28.6,269,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(18.4,7.9,7.8,217,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(608,26,26.1,167,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(15.3,6.6,6.5,124,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(9.64,4.1,4.2,39,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(42.1,18.1,18,39,0,0,ax=ax,marker='v',markersize=markersize)
    plot_scott15(42.1,18.1,18.0,100,0,0,ax=ax,marker='v',markersize=markersize)
    # Busch et al. 2014
    if legend:
        plot_scott15(15.2,0,0,6,0,0,ax=ax,marker='<',label='Busch et al. 2014; z:0',markersize=markersize)
    else:
        plot_scott15(15.2,0,0,6,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(53.5,0,0,39.8,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(102,0,0,311,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(11.7,0,0,39.8,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(32.7,0,0,166,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(63.5,0,0,155,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(10.3,0,0,15.5,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(27.2,0,0,135,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(26.2,0,0,195,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(32.1,0,0,490,0,0,ax=ax,marker='<',markersize=markersize)
    plot_scott15(164,0,0,1000,0,0,ax=ax,marker='<',markersize=markersize)
    
    if legend:
        plot_sanghvi14(10.22,7.86,ax=ax,label='Sanghvi et al. 2014; z:0.5-1',markersize=markersize)
    else:
        plot_sanghvi14(10.22,7.86,ax=ax,markersize=markersize)
    plot_sanghvi14(10.53,8.15,ax=ax,markersize=markersize)
    plot_sanghvi14(11.22,8.02,ax=ax,markersize=markersize)
    plot_sanghvi14(11.27,8.07,ax=ax,markersize=markersize)
    plot_sanghvi14(10.58,8.02,ax=ax,markersize=markersize)
    plot_sanghvi14(10.11,7.92,ax=ax,markersize=markersize)
    plot_sanghvi14(10.84,7.83,ax=ax,markersize=markersize)
    plot_sanghvi14(9.86,7.57,ax=ax,markersize=markersize)
    plot_sanghvi14(10.11,7.69,ax=ax,markersize=markersize)
    plot_sanghvi14(10.97,7.99,ax=ax,markersize=markersize)
    plot_sanghvi14(10.39,7.29,ax=ax,markersize=markersize)
    plot_sanghvi14(10.57,7.50,ax=ax,markersize=markersize)
    plot_sanghvi14(10.13,7.75,ax=ax,markersize=markersize)
    plot_sanghvi14(9.77,7.68,ax=ax,markersize=markersize)
    plot_sanghvi14(10.67,7.69,ax=ax,markersize=markersize)
    plot_sanghvi14(10.17,7.46,ax=ax,markersize=markersize)
    plot_sanghvi14(11.44,7.98,ax=ax,markersize=markersize)
    plot_sanghvi14(9.88,7.38,ax=ax,markersize=markersize)
    plot_sanghvi14(10.47,7.83,ax=ax,markersize=markersize)
    
    if legend:
        plot_scott13(0.91,0.07,0.07,0.0426,0.0014,0.0014,ax=ax,marker='>',color='red',label='Mikly Way',markersize=markersize)
    else:
        plot_scott13(0.91,0.07,0.07,0.0426,0.0014,0.0014,ax=ax,marker='>',color='red',markersize=markersize)
    #x = np.linspace(8,12,markersize=markersize)
    #ax.plot(x,wyithe06(x),'k--')#,label='Wyithe 2006; z:0',markersize=markersize)

def plot_obsMagorrian_legend(ax,markersize=10):
    plot_scott13(69,32,59,39,5,4,ax=ax,label='Scott et al. 2015; z:0',markersize=markersize)
    plot_scott15(1.2,0,0,2.0,0.3,0.3,ax=ax,marker='^',label='Thornton et al. 2008; z:0',markersize=markersize)
    plot_jiang11(9.25,5.7,ax=ax,marker='o',label='Jiang et al. 2011; z:0',markersize=markersize)
    plot_scott15(5.4,0,0,16,0,0,ax=ax,marker='o',label='Jiang et al. 2013; z:0',markersize=markersize)
    plot_scott15(8.4,3.5,3.6,12.2,0,0,ax=ax,marker='D',label='Yuan et al. 2014; z:0',markersize=markersize)
    plot_scott15_empty(2.57,0,0,5.0,0,0,ax=ax,marker='s',label='Reines et al. 2013; z:0 BPT AGNs',markersize=markersize)
    plot_scott15(2.57,0,0,2.5,0,0,ax=ax,marker='s',label='Reines et al. 2013; z:0 BPT composites',markersize=markersize)
    plot_scott15(24.2,10.4,10.4,71,0,0,ax=ax,marker='v',label='Mathur et al. 2012; z:0',markersize=markersize)
    plot_scott15(15.2,0,0,6,0,0,ax=ax,marker='<',label='Busch et al. 2014; z:0',markersize=markersize)
    plot_sanghvi14(10.22,7.86,ax=ax,label='Sanghvi et al. 2014; z:0.5-1',markersize=markersize)
    plot_scott13(0.91,0.07,0.07,0.0426,0.0014,0.0014,ax=ax,marker='>',color='red',label='Mikly Way',markersize=markersize)

def plot_obsMagorrian(ax,z,markersize=8,legend=True,hubble_h=None,silent=False):
    _plot_obsMagorrian = {0.0: plot_obsMagorrian_z0pt0}
    try:
        (_plot_obsMagorrian[z])(ax,markersize=markersize,legend=legend,hubble_h=hubble_h)
    except:
        if not silent:
            print 'available Magorrian relation at redshift:',_plot_obsMagorrian.keys()
        return -1
