#!/usr/bin/env python
import numpy as np

__author__ = ('Yuxiang Qin')                     
__email__ = ("Yuxiang.L.Qin@Gmail.com")          
__all__  = ('plot_obsBHMF', 'plot_obsBHMF_legend',)

def _h_convertor_y(hubble_h_obs,hubble_h=0.678):
    return (hubble_h/hubble_h_obs)**3.

def plot_obsBHMF_z0pt0(ax,hubble_h=0.678,markersize=8,legend=True):

    # Shankar et al. 2009
    file = '/home/yqin/dragons/data/bhmf_Shankar2009.txt'
    data0 = np.recfromtxt(file)
    data = data0[data0[:,0]==0.02][:,1:]
    if legend:
        ax.plot(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),color='lightgray',lw=5,alpha=0.8,label='Shankar et al. 2009')
    #ax.errorbar(data[:,0],10**data[:,1],yerr=[(10**data[:,1]-10**data[:,2])*_h_convertor_y(0.7,hubble_h),(10**data[:,3]-10**data[:,1]])*_h_convertor_y(0.7,hubble_h),color='green')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,2]*_h_convertor_y(0.7,hubble_h),lw=0,color='lightgray',alpha=0.8)#,hatch='///',edgecolor='k')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,3]*_h_convertor_y(0.7,hubble_h),lw=0,color='lightgray',alpha=0.8)#,hatch='///',edgecolor='k')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,4]*_h_convertor_y(0.7,hubble_h),lw=0,color='lightgray',alpha=0.8)#,hatch='///',edgecolor='k')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,5]*_h_convertor_y(0.7,hubble_h),lw=0,color='lightgray',alpha=0.8)#,hatch='///',edgecolor='k')

    # Graham et al. 2007
    file = '/home/yqin/dragons/data/bhmf_Graham2007.txt'
    data = np.recfromtxt(file)
    if legend:
        ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.7,hubble_h),\
                     yerr=[data[:,3]*1e-4*_h_convertor_y(0.7,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.7,hubble_h)],\
                     markersize=markersize,color='gray',fmt='o',label='Graham et al. 2007')
    else:
        ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.7,hubble_h),\
                     yerr=[data[:,3]*1e-4*_h_convertor_y(0.7,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.7,hubble_h)],\
                     markersize=markersize,color='gray',fmt='o')

    # Vika et al. 2009
    file = '/home/yqin/dragons/data/bhmf_Vika2009.txt'
    data = np.recfromtxt(file)
    if legend:
        ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.7,hubble_h),\
                     yerr=[data[:,3]*1e-4*_h_convertor_y(0.7,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.7,hubble_h)],\
                     markersize=markersize,color='gray',fmt='o',mfc='None',label='Vika et al. 2009')
    else:
        ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.7,hubble_h),\
                     yerr=[data[:,3]*1e-4*_h_convertor_y(0.7,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.7,hubble_h)],\
                     markersize=markersize,color='gray',fmt='o',mfc='None')

    # Davis et al. 2014
    file = '/home/yqin/dragons/data/bhmf_Davis2014.txt'
    data = np.recfromtxt(file)
    if legend:
        ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.6777,hubble_h),\
                     yerr=[data[:,3]*1e-4*_h_convertor_y(0.6777,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.6777,hubble_h)],\
                     markersize=markersize,color='gray',fmt='v',label='Davis et al. 2014')
    else:
        ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.6777,hubble_h),\
                     yerr=[data[:,3]*1e-4*_h_convertor_y(0.6777,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.6777,hubble_h)],\
                     markersize=markersize,color='gray',fmt='v')

    # Pakdil et al. 2016
    file = '/home/yqin/dragons/data/bhmf_Pakdil2016.txt'
    data = np.recfromtxt(file)
    if legend:
        ax.errorbar(data[:,0],10**data[:,4]*_h_convertor_y(0.6777,hubble_h),\
                     yerr=[(10**(data[:,4]+data[:,5])-10**data[:,4])*_h_convertor_y(0.6777,hubble_h),\
                           (10**data[:,4]-10**(data[:,4]-data[:,6]))*_h_convertor_y(0.6777,hubble_h)],\
                     color='gray',fmt='s',markersize=markersize,label='Mutlu-Pakdil et al. 2016')
    else:
        ax.errorbar(data[:,0],10**data[:,4]*_h_convertor_y(0.6777,hubble_h),\
                     yerr=[(10**(data[:,4]+data[:,5])-10**data[:,4])*_h_convertor_y(0.6777,hubble_h),\
                           (10**data[:,4]-10**(data[:,4]-data[:,6]))*_h_convertor_y(0.6777,hubble_h)],\
                     color='gray',fmt='s',markersize=markersize)
    
    
    
def plot_obsBHMF_z0pt5(ax,hubble_h=0.678,markersize=10,legend=True):

    # Shankar et al. 2009
    file = '/home/yqin/dragons/data/bhmf_Shankar2009.txt'
    data0 = np.recfromtxt(file)
    data = data0[data0[:,0]==0.50][:,1:]
    if legend:
        ax.plot(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),color='gray',alpha=0.5,lw=5)
    #ax.errorbar(data[:,0],10**data[:,1],yerr=[(10**data[:,1]-10**data[:,2])*_h_convertor_y(0.7,hubble_h),(10**data[:,3]-10**data[:,1]])*_h_convertor_y(0.7,hubble_h),color='green')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,2]*_h_convertor_y(0.7,hubble_h),lw=0,color='gray',alpha=0.5)#,hatch='\\\\\\',edgecolor='k')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,3]*_h_convertor_y(0.7,hubble_h),lw=0,color='gray',alpha=0.5)#,hatch='\\\\\\',edgecolor='k')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,4]*_h_convertor_y(0.7,hubble_h),lw=0,color='gray',alpha=0.5)#,hatch='\\\\\\',edgecolor='k')
    ax.fill_between(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),10**data[:,5]*_h_convertor_y(0.7,hubble_h),lw=0,color='gray',alpha=0.5)#,hatch='\\\\\\',edgecolor='k')

def plot_obsBHMF_legend(ax,markersize=10):
    hubble_h=0.678

    # Shankar et al. 2009
    file = '/home/yqin/dragons/data/bhmf_Shankar2009.txt'
    data0 = np.recfromtxt(file)
    data = data0[data0[:,0]==0.02][:,1:]
    ax.plot(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),color='lightgray',lw=5,alpha=0.8,label='Shankar et al. 2009; z=0')

    file = '/home/yqin/dragons/data/bhmf_Shankar2009.txt'
    data0 = np.recfromtxt(file)
    data = data0[data0[:,0]==0.50][:,1:]
    ax.plot(data[:,0],10**data[:,1]*_h_convertor_y(0.7,hubble_h),lw=5,color='gray',alpha=0.8,label='Shankar et al. 2009;z= 0.5')

    # Graham et al. 2007
    file = '/home/yqin/dragons/data/bhmf_Graham2007.txt'
    data = np.recfromtxt(file)
    ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.7,hubble_h),\
                yerr=[data[:,3]*1e-4*_h_convertor_y(0.7,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.7,hubble_h)],\
                markersize=markersize,color='gray',fmt='o',label='Graham et al. 2007')

    # Vika et al. 2009
    file = '/home/yqin/dragons/data/bhmf_Vika2009.txt'
    data = np.recfromtxt(file)
    ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.7,hubble_h),\
                yerr=[data[:,3]*1e-4*_h_convertor_y(0.7,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.7,hubble_h)],\
                markersize=markersize,color='gray',fmt='o',mfc='None',label='Vika et al. 2009')

    # Davis et al. 2014
    file = '/home/yqin/dragons/data/bhmf_Davis2014.txt'
    data = np.recfromtxt(file)
    ax.errorbar(data[:,0],data[:,1]*1e-4*_h_convertor_y(0.6777,hubble_h),\
                yerr=[data[:,3]*1e-4*_h_convertor_y(0.6777,hubble_h),data[:,2]*1e-4*_h_convertor_y(0.6777,hubble_h)],\
                markersize=markersize,color='gray',fmt='v',label='Davis et al. 2014')

    # Pakdil et al. 2016
    file = '/home/yqin/dragons/data/bhmf_Pakdil2016.txt'
    data = np.recfromtxt(file)
    ax.errorbar(data[:,0],10**data[:,4]*_h_convertor_y(0.6777,hubble_h),\
                yerr=[(10**(data[:,4]+data[:,5])-10**data[:,4])*_h_convertor_y(0.6777,hubble_h),\
                      (10**data[:,4]-10**(data[:,4]-data[:,6]))*_h_convertor_y(0.6777,hubble_h)],\
                color='gray',fmt='s',markersize=markersize,label='Mutlu-Pakdil et al. 2016')

def plot_obsBHMF(ax,z,hubble_h=0.678,markersize=10,legend=True,silent=False):
    _plot_obsBHMF = {0.5: plot_obsBHMF_z0pt5,\
                     0.0: plot_obsBHMF_z0pt0}
    try:
        (_plot_obsBHMF[z])(ax,hubble_h=hubble_h,markersize=markersize,legend=legend)
    except:
        if not silent:
            print 'available BHMF at redshift:', _plot_obsBHMF.keys()
        return -1
