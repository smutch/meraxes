#!/usr/bin/env python
import numpy as np
from astropy import constants
from _meraxes_util import _MBvega2M1450, _M14502MB

__author__ = ('Yuxiang Qin')                       
__email__ = ("Yuxiang.L.Qin@Gmail.com")
__all__  = ('plot_obsQLF', 'plot_obsQLF_legend',)

eta = 0.06
alpha_q = 1.57 
alpha_q_op = 0.44
Cha2Sal = 1.65
Sal2Cha = 1./Cha2Sal
Kro2Sal = 1.5
ergs_to_Lsun = 1./constants.L_sun.to('erg s**-1').value
_1450_to_912 = (1200./1450)**alpha_q_op*(912/1200.)**alpha_q
_Mi_to_Mg = -2.5*alpha_q_op*np.log10(4670./7471)
#Croom et al. 2009
_MB_to_Mg = -2.5*(1.+alpha_q_op)*np.log10(3)-0.14+0.209
_Mg_to_MB = -_MB_to_Mg
_MABB_to_MB = 0.09
_MB_to_MABB = -_MABB_to_MB


def _h_convertor_x(hubble_h_obs,hubble_h=0.678):   
    return -5. * np.log10(hubble_h/hubble_h_obs)
                                                    
def _h_convertor_y(hubble_h_obs,hubble_h=0.678):    
    return (hubble_h/hubble_h_obs)**3. 

# Wolf et al. 2003
def plot_wolf03(x,y,ax, markersize, marker='*',label=None,hubble_h=0.678):
    ax.errorbar([x+_MB_to_MABB+_h_convertor_x(0.65,hubble_h)],[10**y*_h_convertor_y(0.65,hubble_h)],fmt=marker,label=label,color='black',mfc='None',markersize=markersize)

# Hunt et al. 2004
def plot_hunt04(x,y,yerru,yerrl,ax, markersize, marker='o',label=None,hubble_h=0.678):
    ax.errorbar([x+_h_convertor_x(0.5,hubble_h)],[y*_h_convertor_y(0.5,hubble_h)],yerr=[[(yerru-y)*_h_convertor_y(0.5,hubble_h)],[(y-yerrl)*_h_convertor_y(0.5,hubble_h)]],fmt=marker,label=label,color='black',mfc='None',markersize=markersize)

# Dijkstra et al. 2006
# Ly-alpha constraints
def plot_dijkstra06(x,y,ax, markersize, marker='o',color='gray',label=None,hubble_h=0.678):
    MB = +5.48-2.5*np.log10(x*1e9)
    M1450 =_MBvega2M1450(MB)
    dndM1450 = y*1e-6/2.5*x*np.log(10)
    ax.errorbar([M1450+_h_convertor_x(0.7,hubble_h)], [dndM1450*_h_convertor_y(0.7,hubble_h)], yerr=[dndM1450*0.9], uplims=True, markersize=markersize, marker=marker,label=label,color=color)

# Bongiorno et al. 2007 (VVDS)
# AGN
def plot_bongiorno07(x,y,yerru,yerrl,ax, markersize, marker='o',label=None,hubble_h=0.678):
    ax.errorbar([x+_h_convertor_x(0.7,hubble_h)],[y*_h_convertor_y(0.7,hubble_h)],yerr=[[(yerru-y)*_h_convertor_y(0.7,hubble_h)],[(y-yerrl)*_h_convertor_y(0.7,hubble_h)]], markersize=markersize, marker=marker,label=label,color='gray')

# Croom et al. 2009
def plot_croom09(x,y,yerrl,yerru,ax, markersize, marker='.',label=None,hubble_h=0.678):
    ax.errorbar([x+_Mg_to_MB+_h_convertor_x(0.7,hubble_h)],[10**y*_h_convertor_y(0.7,hubble_h)],yerr=[[(10**(y+yerru)-10**y)*_h_convertor_y(0.7,hubble_h)],[(10**y-10**(y-yerrl))*_h_convertor_y(0.7,hubble_h)]], markersize=markersize, marker=marker,label=label,color='gray')

# Willott et al. 2010 (CHQS)
# quasars
def plot_willott10(x,y,yerru,yerrl,ax, markersize, marker='v',label=None,hubble_h=0.678):
    ax.errorbar([x]+_h_convertor_x(0.7,hubble_h),[y*_h_convertor_y(0.7,hubble_h)],yerr=[[(yerru-y)*_h_convertor_y(0.7,hubble_h)],[(y-yerrl)*_h_convertor_y(0.7,hubble_h)]],fmt=marker,label=label,color='black',mfc='None',markersize=markersize)

# Gilkman et al. 2011 
def plot_gilkman11(x,y,yerrl,yerrr,ax, markersize, marker='s',label=None,hubble_h=0.678):
    ax.errorbar([x+_h_convertor_x(0.7,hubble_h)],[y*1e-8*_h_convertor_y(0.7,hubble_h)],yerr=[[yerrl*1e-8*_h_convertor_y(0.7,hubble_h)],[yerrr*1e-8*_h_convertor_y(0.7,hubble_h)]],fmt=marker,label=label,color='black',mfc='None',markersize=markersize)

# Masters et al. 2012 (COSMOS)
# quasars
def plot_masters12(xl,xr,y,yerr,ax, markersize, marker='+',label=None,hubble_h=0.678):
    x = (xl+xr)/2.
    ax.errorbar([x+_h_convertor_x(0.7,hubble_h)],[y*1e-7*_h_convertor_y(0.7,hubble_h)],xerr=[x-xl],yerr=[yerr*1e-7*_h_convertor_y(0.7,hubble_h)], markersize=markersize, marker=marker,label=label,color='black')

# McGreer et al. 2013 (SDSS, UKIDSS, MMT) 
# Shen & Kelly 2012 (SDSS DR7)
# Quasars
def plot_mcgreer13(x,y,yerr,ax, markersize, marker='o',label=None,hubble_h=0.678):
    ax.errorbar([x+_h_convertor_x(0.7,hubble_h)],[10**y*_h_convertor_y(0.7,hubble_h)],yerr=[yerr*1e-9*_h_convertor_y(0.7,hubble_h)],fmt=marker,label=label,color='black',mfc='None',markersize=markersize)

# Palanque-Delabrouille et al. 2013
def plot_magneville13(x,y,yerr,ax, markersize, marker='^',label=None,hubble_h=0.678):
    ax.errorbar([x+_Mg_to_MB+_h_convertor_x(0.71,hubble_h)],[10**y*_h_convertor_y(0.71,hubble_h)],yerr=[(10**(y+yerr)-10**y)*_h_convertor_y(0.71,hubble_h)], markersize=markersize, marker=marker,label=label,color='gray')

# Giallongo et al. 2012 
def plot_giallongo12(x,y,yerru,yerrl,ax, markersize, marker='D',label=None,hubble_h=0.678):
    ax.errorbar([x+_h_convertor_x(0.7,hubble_h)],[10**y*_h_convertor_y(0.7,hubble_h)],yerr=[[10**((yerrl+y)/2.0)*(y-yerrl)*np.log(10)*_h_convertor_y(0.7,hubble_h)],[10**((yerru+y)/2.0)*(yerru-y)*np.log(10)*_h_convertor_y(0.7,hubble_h)]],markersize=markersize, marker=marker,label=label,color='gray')

# Giallongo et al. 2015 (Chandra, HST, Spitzer and various groundbased telescopes)
# AGN 
def plot_giallongo15(x,y,yerru,yerrl,ax, markersize, marker='s',label=None,hubble_h=0.678):
    ax.errorbar([x+_h_convertor_x(0.7,hubble_h)],[10**y*_h_convertor_y(0.7,hubble_h)],yerr=[[10**((yerrl+y)/2.0)*(y-yerrl)*np.log(10)*_h_convertor_y(0.7,hubble_h)],[10**((yerru+y)/2.0)*(yerru-y)*np.log(10)*_h_convertor_y(0.7,hubble_h)]],markersize=markersize, marker=marker,label=label,color='gray')

## Fitting
def fitting_giallongo15(M,beta,gamma,Mbreak,log_phi):
    return 10**log_phi/(10**(0.4*(Mbreak-M)*(beta-1))+10**(0.4*(Mbreak-M)*(gamma-1)))
#x = np.linspace(-15,-28.5)
#axs[2].plot(x,fitting_giallongo15(x,1.52,3.13,-23.2,-5.2),'k--')
#axs[2].plot(x,fitting_giallongo15(x,1.81,3.14,-23.6,-5.7),'k:')
#axs[3].plot(x,fitting_giallongo15(x,1.66,3.35,-23.4,-5.8),'k--')

def plot_obsQLF_z6pt0(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_dijkstra06(14,3.9,ax, markersize=markersize, marker='>',label='Dijkstra & Wyithe 2006; z:~5.7',hubble_h=hubble_h)
        plot_dijkstra06(20,1.2,ax, markersize=markersize, marker='^',label='Dijkstra & Wyithe 2006; z:~6.6',hubble_h=hubble_h)
        plot_willott10(-22.2, 5e-8, 1e-7, 1e-11, ax, markersize=markersize, marker='v',label='Willott et al. 2010; z:~6',hubble_h=hubble_h)
        plot_giallongo15(-19,-4.7,-4.2,-5.7,ax, markersize=markersize, marker='s',label='Giallongo et al. 2015',hubble_h=hubble_h)
    else:
        plot_dijkstra06(14,3.9,ax, markersize=markersize, marker='>',hubble_h=hubble_h)
        plot_dijkstra06(20,1.2,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
        plot_willott10(-22.2, 5e-8, 1e-7, 1e-11, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
        plot_giallongo15(-19,-4.7,-4.2,-5.7,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_willott10(-24.7, 7e-9, 9e-9, 5e-9,  ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_willott10(-25.2, 6e-9, 1e-8, 4e-9,  ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_willott10(-26.0, 2e-9, 3e-9, 1.5e-9,ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_willott10(-26.2, 1e-9, 1.8e-9,5e-10,ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_willott10(-26.8,1.1e-9,1.8e-9,6e-10,ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_willott10(-27.2,3.8e-10,5.5e-10,2e-10,ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_willott10(-27.7,2.1e-10,3.5e-10,1e-10,ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_giallongo15(-20,-4.7,-4.4,-5.0,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_giallongo15(-21,-5.7,-5.2,-6.5,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    
def plot_obsQLF_z5pt0(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_dijkstra06(20,1.6,ax, markersize=markersize, marker='<',label='Dijkstra & Wyithe 2006; z:~4.5',hubble_h=hubble_h)
        plot_dijkstra06(14,3.9,ax, markersize=markersize, marker='>',label='Dijkstra & Wyithe 2006; z:~5.7',hubble_h=hubble_h)
        plot_mcgreer13(-28.05,-9.45,0.21,ax, markersize=markersize, marker='D',label='Shen & Kelly 2012; z:~4.75',hubble_h=hubble_h)
        plot_mcgreer13(-27.00,-8.40,2.81,ax, markersize=markersize, marker='^',label='McGreer et al. 2013; z:4.7-5.1',hubble_h=hubble_h)
        plot_giallongo15(-19,-4.1,-3.80,-4.50,ax, markersize=markersize, marker='s',label='Giallongo et al. 2015',hubble_h=hubble_h)
    else:
        plot_dijkstra06(20,1.6,ax, markersize=markersize, marker='<',hubble_h=hubble_h)
        plot_dijkstra06(14,3.9,ax, markersize=markersize, marker='>',hubble_h=hubble_h)
        plot_mcgreer13(-28.05,-9.45,0.21,ax, markersize=markersize, marker='D',hubble_h=hubble_h)
        plot_mcgreer13(-27.00,-8.40,2.81,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
        plot_giallongo15(-19,-4.1,-3.80,-4.50,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_mcgreer13(-27.55,-9.24,0.26,ax, markersize=markersize, marker='D',hubble_h=hubble_h)
    plot_mcgreer13(-27.55,-9.24,0.26,ax, markersize=markersize, marker='D',hubble_h=hubble_h)
    plot_mcgreer13(-27.05,-8.51,0.58,ax, markersize=markersize, marker='D',hubble_h=hubble_h)
    plot_mcgreer13(-26.55,-8.20,0.92,ax, markersize=markersize, marker='D',hubble_h=hubble_h)
    plot_mcgreer13(-26.05,-7.90,1.89,ax, markersize=markersize, marker='D',hubble_h=hubble_h)
    plot_mcgreer13(-26.45,-7.84,6.97,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
    plot_mcgreer13(-25.90,-7.90,5.92,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
    plot_mcgreer13(-25.35,-7.53,10.23,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
    plot_mcgreer13(-24.80,-7.36,11.51,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
    plot_mcgreer13(-24.25,-7.14,19.90,ax, markersize=markersize, marker='^',hubble_h=hubble_h)
    plot_giallongo15(-20,-4.6,-4.20,-5.00,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_giallongo15(-21,-4.8,-4.50,-5.30,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_giallongo15(-22.5,-5.3,-4.75,-6.25,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    
def plot_obsQLF_z4pt0(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_wolf03(-24, -6.02, ax, markersize=markersize, marker='*',label='Wolf et al. 2003; z:3.6-4.2',hubble_h=hubble_h)
        plot_dijkstra06(20,1.6,ax, markersize=markersize, marker='<',label='Dijkstra & Wyithe 2006; z:~4.5',hubble_h=hubble_h)
        plot_gilkman11(-28.45,0.008,0.005,0.011,ax, markersize=markersize, marker='s',label='Gilkman et al. 2011; z:~4',hubble_h=hubble_h)
        plot_gilkman11(-25.37,24,9,13,ax, markersize=markersize, marker='o',label='Gilkman et al. 2011; z:~4',hubble_h=hubble_h)
        plot_masters12(-25.5,-23.5, 1.9,0.8, ax, markersize=markersize, marker='x',label='Masters et al. 2012; z:~4',hubble_h=hubble_h)
        plot_giallongo15(-19,-4.3,-3.80,-5.25,ax, markersize=markersize, marker='s',label='Giallongo et al. 2015',hubble_h=hubble_h)
    else:
        plot_wolf03(-24, -6.02, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
        plot_dijkstra06(20,1.6,ax, markersize=markersize, marker='<',hubble_h=hubble_h)
        plot_gilkman11(-28.45,0.008,0.005,0.011,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
        plot_gilkman11(-25.37,24,9,13,ax, markersize=markersize, marker='o',hubble_h=hubble_h)
        plot_masters12(-25.5,-23.5, 1.9,0.8, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
        plot_giallongo15(-19,-4.3,-3.80,-5.25,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_wolf03(-25, -6.02, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-26, -6.81, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-27, -6.81, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    # 1. SDSS
    plot_gilkman11(-27.33,0.20,0.03,0.04,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_gilkman11(-26.46,0.93,0.07,0.07,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_gilkman11(-25.70,4.3,0.5,0.5,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_gilkman11(-24.72,0.4,0.2,0.3,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    # 2. NOAO Deep Wide-Field Survey (NDWFS) and the Deep Lens Survey (DLS,hubble_h=hubble_h)
    plot_gilkman11(-24.71,8.8,4.8,8.5,ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_gilkman11(-23.47,143,53,77,ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_gilkman11(-22.61,307,133,208,ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_gilkman11(-21.61,434,280,572,ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_masters12(-23.5,-22.5, 3.5,1.2, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
    plot_masters12(-22.5,-21.5, 6.4,1.4, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
    plot_masters12(-21.5,-20.5, 20.6,2.6, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
    plot_giallongo15(-20,-4.6,-4.25,-4.90,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    plot_giallongo15(-21,-4.9,-4.50,-5.25,ax, markersize=markersize, marker='s',hubble_h=hubble_h)
    

def plot_obsQLF_z3pt0(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_wolf03(-24, -5.42, ax, markersize=markersize, marker='*',label='Wolf et al. 2003; z:2.4-3.0',hubble_h=hubble_h)
        plot_bongiorno07(-20.25,4.5e-6, 8e-6, 1.2e-6, ax, markersize=markersize, marker='o',label='Bongiorno et al. 2007; z:2.1-3.6',hubble_h=hubble_h)
        plot_hunt04(-20.5, 1.5e-6, 3.0e-6,4.5e-6, ax, markersize=markersize, marker='h',label='Hunt et al. 2004; z:3',hubble_h=hubble_h)
        plot_masters12(-25.5,-23.5, 5.3,  2.6, ax, markersize=markersize, marker='x',label='Masters et al. 2012; z:~3.2',hubble_h=hubble_h)
    else:
        plot_wolf03(-24, -5.42, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
        plot_bongiorno07(-20.25,4.5e-6, 8e-6, 1.2e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
        plot_hunt04(-20.5, 1.5e-6, 3.0e-6,4.5e-6, ax, markersize=markersize, marker='h',hubble_h=hubble_h)
        plot_masters12(-25.5,-23.5, 5.3,  2.6, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
    plot_wolf03(-25, -5.69, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-26, -6.01, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-27, -6.38, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-28, -6.38, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_bongiorno07(-21.25,3.5e-6,4.6e-6,2.1e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-22.25,2.0e-6,2.5e-6,1.5e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-23.25,1.8e-6,2.3e-6,1.5e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-24.25,5.0e-7,8.0e-7,2.0e-7, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-25.25,3.0e-7,5.0e-7,9.0e-8, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-26.25,1.8e-7,3.5e-7,1.5e-8, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_hunt04(-21.5, 2.8e-6, 4.5e-6,1.7e-6, ax, markersize=markersize, marker='h',hubble_h=hubble_h)
    plot_hunt04(-23.5, 1.8e-6, 2.6e-6,8.0e-7, ax, markersize=markersize, marker='h',hubble_h=hubble_h)
    plot_hunt04(-24.5, 5.8e-7, 1.3e-6,1.8e-7, ax, markersize=markersize, marker='h',hubble_h=hubble_h)
    plot_hunt04(-25.5, 9.8e-7, 2.0e-6,5.5e-7, ax, markersize=markersize, marker='h',hubble_h=hubble_h)
    plot_masters12(-23.5,-22.5, 11.9, 3.9, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
    plot_masters12(-22.5,-21.5, 39.0, 7.1, ax, markersize=markersize, marker='x',hubble_h=hubble_h)
    plot_masters12(-21.5,-20.5, 59.8, 8.4, ax, markersize=markersize, marker='x',hubble_h=hubble_h)

def plot_obsQLF_z2pt0(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_wolf03(-21.2, -5.53, ax, markersize=markersize, marker='*',label='Wolf et al. 2003; z:1.8-2.4',hubble_h=hubble_h)
        plot_croom09(-29.25, -9.08, 0.34, 0.29, ax, markersize=markersize, marker='.',label='Croom et al. 2009; z:1.82-2.20',hubble_h=hubble_h)
        plot_magneville13(-28.40, -7.40, 0.43, ax, markersize=markersize, marker='v',label='Palanque-Delabrouille et al. 2013; z:1.82-2.20',hubble_h=hubble_h)
    else:
        plot_wolf03(-21.2, -5.53, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
        plot_croom09(-29.25, -9.08, 0.34, 0.29, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
        plot_magneville13(-28.40, -7.40, 0.43, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_wolf03(-22.2, -5.36, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-23.2, -5.60, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-24.2, -5.85, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-25.2, -6.15, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-26.2, -6.85, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_croom09(-28.75, -8.38, 0.13, 0.12, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-28.25, -7.89, 0.07, 0.06, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-27.75, -7.30, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-27.25, -6.88, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-26.75, -6.47, 0.04, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-26.25, -6.18, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.75, -5.93, 0.03, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.25, -5.76, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.75, -5.67, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.25, -5.57, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.75, -5.58, 0.05, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_magneville13(-28.00, -6.89, 0.25, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-27.60, -7.37, 0.43, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-27.20, -6.68, 0.19, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-26.80, -6.33, 0.13, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-26.40, -5.98, 0.09, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-26.00, -5.86, 0.08, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-25.60, -5.88, 0.08, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-25.20, -5.76, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.80, -5.72, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.40, -5.67, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.00, -5.47, 0.06, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-23.60, -5.72, 0.13, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-23.20, -5.36, 0.32, ax, markersize=markersize, marker='v',hubble_h=hubble_h)

def plot_obsQLF_z1pt5(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_wolf03(-21.2, -5.78, ax, markersize=markersize, marker='*',label='Wolf et al. 2003; z:1.2-1.8',hubble_h=hubble_h)
        plot_croom09(-29.75, -9.40, 0.77, 0.52, ax, markersize=markersize, marker='.',label='Croom et al. 2009; z:1.44-1.82',hubble_h=hubble_h)
        plot_magneville13(-27.60, -6.88, 0.25, ax, markersize=markersize, marker='v',label='Palanque-Delabrouille et al. 2013; z:1.44-1.82',hubble_h=hubble_h)
    else:
        plot_wolf03(-21.2, -5.78, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
        plot_croom09(-29.75, -9.40, 0.77, 0.52, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
        plot_magneville13(-27.60, -6.88, 0.25, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_wolf03(-22.2, -5.58, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-23.2, -5.69, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-24.2, -6.02, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-25.2, -6.10, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_wolf03(-26.2, -6.80, ax, markersize=markersize, marker='*',hubble_h=hubble_h)
    plot_croom09(-29.25, -9.54, 0.77, 0.52, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-28.75, -8.69, 0.20, 0.19, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-28.25, -8.09, 0.09, 0.08, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-27.75, -7.58, 0.05, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-27.25, -7.11, 0.03, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-26.75, -6.68, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-26.25, -6.32, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.75, -6.06, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.25, -5.86, 0.03, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.75, -5.66, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.25, -5.54, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.75, -5.47, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.25, -5.49, 0.05, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_magneville13(-27.20, -6.88, 0.25, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-26.80, -7.34, 0.43, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-26.40, -6.56, 0.18, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-26.00, -6.05, 0.10, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-25.60, -5.91, 0.09, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-25.20, -5.82, 0.08, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.80, -5.64, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.40, -5.64, 0.06, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.00, -5.50, 0.06, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-23.60, -5.58, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-23.20, -5.38, 0.10, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-22.80, -5.69, 0.45, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    
def plot_obsQLF_z1pt3(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_bongiorno07(-20., 5e-6, 7e-6, 2e-6, ax, markersize=markersize, marker='o',label='Bongiorno et al. 2007; z:1-1.55',hubble_h=hubble_h)
        plot_croom09(-29.25, -9.31, 0.77, 0.52, ax, markersize=markersize, marker='.',label='Croom et al. 2009; z:1.06-1.44',hubble_h=hubble_h)
        plot_magneville13(-26.40, -6.60, 0.19, ax, markersize=markersize, marker='v',label='Palanque-Delabrouille et al. 2013; z:1.06-1.44',hubble_h=hubble_h)
    else:
        plot_bongiorno07(-20., 5e-6, 7e-6, 2e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
        plot_croom09(-29.25, -9.31, 0.77, 0.52, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
        plot_magneville13(-26.40, -6.60, 0.19, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_bongiorno07(-21., 1.5e-5, 1.8e-5, 1e-5, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-22., 8e-6,1e-5,5e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-23., 5e-6,6e-6,3e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-24., 2e-6,3e-6,1e-6, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_bongiorno07(-26., 6e-7,1e-6,7e-8, ax, markersize=markersize, marker='o',hubble_h=hubble_h)
    plot_croom09(-28.25, -8.65, 0.20, 0.19, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-27.75, -7.96, 0.08, 0.07, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-27.25, -7.65, 0.06, 0.05, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-26.75, -7.11, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-26.25, -6.73, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.75, -6.32, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.25, -6.08, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.75, -5.82, 0.03, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.25, -5.65, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.75, -5.53, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.25, -5.46, 0.02, 0.02, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-22.75, -5.52, 0.04, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-22.25, -4.99, 0.77, 0.52, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_magneville13(-26.00, -6.44, 0.16, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-25.60, -6.99, 0.10, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-25.20, -5.87, 0.09, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.80, -5.79, 0.08, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.40, -5.57, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-24.00, -5.64, 0.08, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-23.60, -5.56, 0.08, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-23.20, -5.30, 0.07, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-22.80, -5.29, 0.12, ax, markersize=markersize, marker='v',hubble_h=hubble_h)
    plot_magneville13(-22.40, -5.22, 0.25, ax, markersize=markersize, marker='v',hubble_h=hubble_h)

def plot_obsQLF_z0pt5(ax,markersize=10,legend=True,hubble_h=0.678):
    if legend:
        plot_giallongo12(_M14502MB(-21.75), -5.82, -5.79, -5.85, ax, markersize=markersize, marker='D', hubble_h=hubble_h,label='Richards et al. 2005')
        plot_croom09(-26.75, -8.28, 0.24, 0.22, ax, markersize=markersize, marker='.',label='Croom et al. 2009; z:0.4-0.68',hubble_h=hubble_h)
    else:
        plot_giallongo12(_M14502MB(-21.75), -5.82, -5.79, -5.85, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
        plot_croom09(-26.75, -8.28, 0.24, 0.22, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-22.05), -5.97, -5.95, -6.01, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-22.35), -6.08, -6.06, -6.10, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-22.65), -6.15, -6.12, -6.17, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-22.95), -6.39, -6.37, -6.43, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-23.25), -6.58, -6.54, -6.61, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-23.55), -6.81, -6.76, -6.85, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-23.85), -7.05, -6.99, -7.12, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-24.15), -7.25, -7.18, -7.34, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-24.45), -7.43, -7.36, -7.54, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_giallongo12(_M14502MB(-24.75), -7.71, -7.60, -7.87, ax, markersize=markersize, marker='D', hubble_h=hubble_h)
    plot_croom09(-26.25, -8.10, 0.18, 0.17, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.75, -7.37, 0.07, 0.06, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-25.25, -7.09, 0.05, 0.05, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.75, -6.71, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-24.25, -6.44, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.75, -6.16, 0.03, 0.03, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-23.25, -6.04, 0.05, 0.05, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-22.75, -5.69, 0.04, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-22.25, -5.59, 0.04, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-21.75, -5.58, 0.04, 0.04, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-21.25, -5.52, 0.05, 0.05, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-20.75, -5.34, 0.08, 0.07, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    plot_croom09(-20.25, -5.11, 0.20, 0.19, ax, markersize=markersize, marker='.',hubble_h=hubble_h)
    
def plot_obsQLF_legend(ax,markersize=10):
    plot_willott10(-22.2, 5e-8, 1e-7, 1e-11, ax, markersize=markersize, marker='v',label='Willott et al. 2010')
    plot_giallongo15(-19,-4.1,-3.80,-4.50,ax, markersize=markersize, marker='s',label='Giallongo et al. 2015')
    plot_dijkstra06(20,1.2,ax, markersize=markersize, marker='^',label='Dijkstra & Wyithe 2006; z~6.6')
    plot_dijkstra06(14,3.9,ax, markersize=markersize, marker='>',label='Dijkstra & Wyithe 2006; z~5.7')
    plot_dijkstra06(20,1.6,ax, markersize=markersize, marker='<',label='Dijkstra & Wyithe 2006; z~4.5')
    plot_mcgreer13(-28.05,-9.45,0.21,ax, markersize=markersize, marker='D',label='Shen & Kelly 2012')
    plot_mcgreer13(-27.00,-8.40,2.81,ax, markersize=markersize, marker='^',label='McGreer et al. 2013')
    plot_wolf03(-24, -6.02, ax, markersize=markersize, marker='*',label='Wolf et al. 2003')
    plot_gilkman11(-28.45,0.008,0.005,0.011,ax, markersize=markersize, marker='s',label='Gilkman et al. 2011; SDSS')
    plot_gilkman11(-25.37,24,9,13,ax, markersize=markersize, marker='o',label='Gilkman et al. 2011; NDWFS&DLS')
    plot_masters12(-25.5,-23.5, 1.9,0.8, ax, markersize=markersize, marker='x',label='Masters et al. 2012')
    plot_bongiorno07(-20.25,4.5e-6, 8e-6, 1.2e-6, ax, markersize=markersize, marker='o',label='Bongiorno et al. 2007')
    plot_hunt04(-20.5, 1.5e-6, 3.0e-6,4.5e-6, ax, markersize=markersize, marker='h',label='Hunt et al. 2004')
    plot_croom09(-29.25, -9.08, 0.34, 0.29, ax, markersize=markersize, marker='.',label='Croom et al. 2009')
    plot_magneville13(-28.40, -7.40, 0.43, ax, markersize=markersize, marker='v',label='Palanque-Delabrouille et al. 2013')
    plot_giallongo12(_M14502MB(-21.75), -5.82, -5.79, -5.85, ax, markersize=markersize, marker='D',label='Richards et al. 2005')

def plot_obsQLF(ax,z,markersize=10,legend=True,hubble_h=0.678,silent=False):
    _plot_obsQLF = {6.0: plot_obsQLF_z6pt0,\
                    5.0: plot_obsQLF_z5pt0,\
                    4.0: plot_obsQLF_z4pt0,\
                    3.0: plot_obsQLF_z3pt0,\
                    2.0: plot_obsQLF_z2pt0,\
                    1.5: plot_obsQLF_z1pt5,\
                    1.3: plot_obsQLF_z1pt3,\
                    0.5: plot_obsQLF_z0pt5}
    try:
        (_plot_obsQLF[z])(ax, markersize=markersize,legend=legend,hubble_h=hubble_h)
    except:    
        if not silent:
            print 'available QLF at redshift:', _plot_obsQLF.keys()
        return -1
