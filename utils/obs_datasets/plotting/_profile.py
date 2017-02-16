#!/usr/bin/env python
import numpy as np
from scipy import constants
import random


#hubbe (km/s)*Mpc**-1
G = 4.302e-6 #Kpc*Msol**-1.0*(km/s)**2

def _delta_c(c):
    return 66.66666667*(c**3)/(np.log(1+c)-c/(1+c))

def _rho_c(h,omegaL0,omegaM0,omegaB0,z):
    return (h*100*1e-3)**2*(omegaM0*((1+z)**3)+omegaB0*((1+z)**2)+omegaL0)*3.0/(8.0*constants.pi*G)

def _nfw(rho_c, r_200, c, xrange=(0.5,2.0)):
    delta_c = _delta_c(c)
    r_s = r_200/c
    x = np.linspace(xrange[0],xrange[1],1000)
    y = np.log10(rho_c*delta_c)-x+np.log10(r_s)-2.0*np.log10(1.0+(10.0**x)/r_s)
    return x,y

def _nfw_mass(rho_c, r_200, c, xrange=(0.5,2.0)):
    delta_c = _delta_c(c)
    r_s = r_200/c
    x = np.linspace(xrange[0],xrange[1],1000)
    y = np.log10(4*constants.pi*rho_c*delta_c)+3.0*np.log10(r_s)+np.log10(np.log((r_s+(10.0**x))/r_s)-(10.0**x)/(r_s+(10.0**x)))
    return x,y

def _sis(sigma_v,xrange=(0.5,2.0)):
    x = np.linspace(xrange[0],xrange[1],1000)
    y = 2.0*np.log10(sigma_v)-2.0*x-np.log10(2.0*constants.pi*G)
    return x,y

def _grazian(xrange=(9.0,12.0),z=5,return_error=0):
    if((z>=3.5) & (z<4.5)):
        alpha = -1.63
        delta_alpha = 0.05
        log_M = 10.96
        delta_log_M = 0.13
        log_phi = -3.94
        delta_log_phi = 0.16
    elif((z>=4.5) & (z<5.5)):
        alpha = -1.63
        delta_alpha = 0.09
        log_M = 10.78
        delta_log_M = 0.23
        log_phi = -4.18
        delta_log_phi = 0.29
    elif((z>=5.5) & (z<6.5)):
        alpha = -1.55
        delta_alpha = 0.19
        log_M = 10.49
        delta_log_M = 0.32
        log_phi = -4.16
        delta_log_phi = 0.47
    elif((z>=6.5) & (z<7.5)):
        alpha = -1.88
        delta_alpha = 0.36
        log_M = 10.69
        delta_log_M = 1.58
        log_phi = -5.24
        delta_log_phi = 2.02
    else:
        print 'Given redshift is beyond the limits of Grazian 15'
        return -1
    x = np.linspace(xrange[0],xrange[1],1000)
    y = log_phi + (alpha+1)*(x-log_M) - 0.4343*10**(x-log_M)
    if return_error:
        delta_y = delta_log_phi + abs(x-log_M)*delta_alpha + abs(-alpha-1 + 10**(x-log_M))*delta_log_M
        return x,y,delta_y
    else:
        return x,y

from _function import _sigmoid

def _Sawala(xrange=[7.5,11],M_t=10**11.6,a=0.69,b=0.98,w=0.79):
    x = np.linspace(xrange[0],xrange[1],1000)
    y = a + (b-a) * _sigmoid(10**x/M_t,w)
    return x,y

def _Schaller(xrange=[7.5,11]):
    M_12 = 10**11.33
    M_23 = 10**13.19
    x = np.linspace(xrange[0],xrange[1],1000)
    y = 0.7309 + 0.1123 * _sigmoid(10**x/M_12,1.721) + 0.1625* _sigmoid(10**x/M_23,2.377)
    return x,y
