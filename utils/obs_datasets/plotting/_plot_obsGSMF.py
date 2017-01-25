#!/usr/bin/env python
import numpy as np

__author__ = ('Yuxiang Qin')
__email__ = ("Yuxiang.L.Qin@Gmail.com")
__all__  = ('plot_obsGSMF', 'plot_obsGSMF_legend',)

eta = 0.06
alpha_q = 1.57 
alpha_q_op = 0.44
Cha2Sal = 1.65
Sal2Cha = 1./Cha2Sal
Kro2Sal = 1.5

def _Schechter(x, logM_star, alpha1, phi_star1, alpha2=0.0, phi_star2=0.0):
    return np.log(10)*(phi_star1*(10**((x-logM_star)*(1.+alpha1)))+phi_star2*(10**((x-logM_star)*(1.+alpha2))))*np.exp(-10**(x-logM_star))

def _h_convertor_x(hubble_h_obs,hubble_h=0.678):
    return -1. * np.log10(hubble_h/hubble_h_obs)

def _h_convertor_x2(hubble_h_obs,hubble_h=0.678):
    return -1. * np.log10(hubble_h/hubble_h_obs**2.)

def _h_convertor_y(hubble_h_obs,hubble_h=0.678):
    return 3. * np.log10(hubble_h/hubble_h_obs)

def plot_obsGSMF_z7pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Gonzalez 11
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Gonzalez11_z7_smf.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]+_h_convertor_y(0.69999999,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='s',color='gray',label='Gonzalez et al. 2011')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]+_h_convertor_y(0.69999999,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='s',color='gray')
    
    #Duncan 14
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Duncan14_MF_z7.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[8:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h)+np.log10(Cha2Sal),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                          np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='o',color='gray',label='Duncan et al. 2014')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h)+np.log10(Cha2Sal),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                          np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='o',color='gray')
    
    #Song 15
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Song2015_z7_smf.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[19:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='D',color='gray',label='Song et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='D',color='gray')
    
    #Grazian 15
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFgrazian2015_z7.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                          np.log10(data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='^',color='gray',label='Grazian et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                          np.log10(data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='^',color='gray')


def plot_obsGSMF_z5pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Gonzalez 11
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Gonzalez11_z5_smf.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]+_h_convertor_y(0.69999999,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='s',color='gray',label='Gonzalez et al. 2011')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]+_h_convertor_y(0.69999999,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='s',color='gray')
    
    #Duncan 14
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Duncan14_MF_z5.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[8:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h)+np.log10(Cha2Sal),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                          np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='o',color='gray',label='Duncan et al. 2014')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h)+np.log10(Cha2Sal),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                          np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='o',color='gray')
    
    
    #Song 15
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Song2015_z5_smf.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[19:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='D',color='gray',label='Song et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='D',color='gray')
    
    #Grazian 15
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFgrazian2015_z5.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                          np.log10(data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='^',color='gray',label='Grazian15')
    else: 
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                          np.log10(data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='^',color='gray')
    
    #Stefanon 16
    data = []
    file = '/home/yqin/dragons/data/smf_Stefanon2016_z5.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    data[:,0] +=np.log10(Cha2Sal)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])-5.+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,3]),\
                          np.log10(data[:,1]+data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='p',color='gray',mfc='None',\
                    label='Stefanon et al. 2016')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])-5.+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,3]),\
                          np.log10(data[:,1]+data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='p',color='gray',mfc='None')
    
    # Davidzon 17 
    #z: 4.5-5.5
    logM_star = 11.30
    alpha = -2.11
    phi_star = 0.003e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)

def plot_obsGSMF_z4pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Duncan 14
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Duncan14_MF_z4.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[8:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h)+np.log10(Cha2Sal),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                          np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='o',color='gray',label='Duncan et al. 2014')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h)+np.log10(Cha2Sal),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                          np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='o',color='gray')
    
    #Grazian 15
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFgrazian2015_z4.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                          np.log10(data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='^',color='gray',label='Grazian15')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                          np.log10(data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='^',color='gray')
    
    #Stefanon 16
    data = []
    file = '/home/yqin/dragons/data/smf_Stefanon2016_z4.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    data[:,0] +=np.log10(Cha2Sal)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])-5.+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,3]),\
                          np.log10(data[:,1]+data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='p',color='gray',mfc='None',\
                    label='Stefanon et al. 2016')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,1])-5.+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,3]),\
                          np.log10(data[:,1]+data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='p',color='gray',mfc='None')
    
    # Davidzon 17 
    #z: 3.5-4.5
    logM_star = 11.10
    alpha = -1.98
    phi_star = 0.016e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)

    
def plot_obsGSMF_z2pt0(ax,hubble_h=0.678,markersize=10,legend=True):
        
    #Marchesini 09
    data = []
    file = '/home/yqin/dragons/data/smf_Marchesini2009_z2pt5.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    data[:,2] +=np.log10(Kro2Sal)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,3],data[:,2]],\
                    markersize=markersize,fmt='*',color='gray',label='Marchesini et al. 2009')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,3],data[:,2]],\
                    markersize=markersize,fmt='*',color='gray')
    
    # Mortlock 11
    data = []
    file = '/home/yqin/dragons/data/smf_Mortlock2011_z2pt25.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray',label='Mortlock et al. 2011')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray')
    
    #Khostovan 16
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFkhostovan16_z2.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='>',color='gray',label='Khostovan16_HbetaOIII')
    else:
        ax.errorbar(data[:,0]+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='>',color='gray')
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFkhostovan16OII_z2.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='>',color='gray',label='Khostovan16_OII')
    else:
        ax.errorbar(data[:,0]+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='>',color='gray')
    
    # Pozzetti 07
    # K selected complex
    # z~1.6-2.5
    alpha = -1.17
    phi_star = 1.25e-3
    logM_star = 10.97
    M_star = 10**logM_star
    x = np.linspace(10.0,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',\
                lw=2,label='Pozzetti et al. 2007',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',lw=2,alpha=1)
    
    # Santini 12 #z:1.8-2.5
    alpha = -1.58
    logM_star = 11.29
    logphi_star = -3.52
    M_star = 10**logM_star
    phi_star = 10**logphi_star
    x = np.linspace(8,13)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=5,label='Santini et al. 2012',alpha=1)
    else:
        ax.plot(x+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=5,alpha=1)
    
    # Ilbert 13 #z: 2-2.5
    alpha1 = -0.22
    phi_star1 = 0.62e-3
    alpha2 = -1.6
    phi_star2 = 0.15e-3
    logM_star = 10.74
    M_star = 10**logM_star
    x = np.linspace(10.04,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',\
                lw=2,label='Ilbert et al. 2013',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=2,alpha=1)
    
    # Muzzin 13 #z: 2-2.5
    alpha = -0.55
    logM_star = 10.81
    M_star = 10**logM_star
    phi_star = 4.79e-4
    x = np.linspace(10.54,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',\
                lw=5,label='Muzzin et al. 2013',alpha=1)
    else:
       ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',lw=5,alpha=1)
    
    # Tomczak 14 #z:2.0-2.5
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFtomczak_z2pt25.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,2],data[:,3]],\
                    markersize=markersize,fmt='v',color='gray',label='Tomczak et al. 2014')
    else:
        ax.errorbar(data[:,0]+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    markersize=markersize,fmt='v',color='gray')
   
    # Huertas-Company 16 #z:2-2.5
    alpha = -1.19
    phi_star = 0.41e-3
    logM_star = 11.05
    M_star = 10**logM_star
    x = np.linspace(10.21,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=5,label='Huertas-Company et al. 2016',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=5,alpha=1)
    
    # Davidzon 17 
    #z: 2.0-2.5
    logM_star = 10.60
    alpha1 = -1.57
    phi_star1 = 0.295e-3
    alpha2 = 0.07
    phi_star2 = 0.45e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)


def plot_obsGSMF_z1pt75(ax,hubble_h=0.678,markersize=10,legend=True):
    
    # Marchesini et al. 2009
    data = []
    file = '/home/yqin/dragons/data/smf_Marchesini2009_z1pt65.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    data[:,2] +=np.log10(Kro2Sal)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,3],data[:,2]],\
                    markersize=markersize,fmt='*',color='gray',label='Marchesini et al. 2009')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,3],data[:,2]],\
                    markersize=markersize,fmt='*',color='gray')
    
    # Mortlock 11
    data = []
    file = '/home/yqin/dragons/data/smf_Mortlock2011_z1pt75.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray',label='Mortlock et al. 2011')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray')
    
    # Santini 12 
    #z:1.4-1.8
    alpha = -1.42
    logM_star = 11.32
    logphi_star = -3.41
    M_star = 10**logM_star
    phi_star = 10**logphi_star
    x = np.linspace(8,13)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
       ax.plot(x+_h_convertor_x2(0.7,hubble_h),\
               np.log10(y)+_h_convertor_y(0.7,hubble_h),\
               color='gray',linestyle='--',lw=5,label='Santini et al. 2012',alpha=1)
    else:
       ax.plot(x+_h_convertor_x2(0.7,hubble_h),\
               np.log10(y)+_h_convertor_y(0.7,hubble_h),\
               color='gray',linestyle='--',lw=5,alpha=1)
    
    # Ilbert 13 
    #z: 1.5-2
    alpha1 = -0.24
    phi_star1 = 0.88e-3
    alpha2 = -1.6
    phi_star2 = 0.33e-3
    logM_star = 10.74
    M_star = 10**logM_star
    x = np.linspace(9.5,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',\
                lw=2,label='Ilbert et al. 2013',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=2,alpha=1)

    # Muzzin 13 
    #z: 1.5-2
    alpha = -0.86
    logM_star = 10.81
    M_star = 10**logM_star
    phi_star = 10.13e-4
    x = np.linspace(10.03,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',\
                lw=5,label='Muzzin et al. 2013',alpha=1)
    else:
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',lw=5,alpha=1)
    
    # Huertas-Company 16
    #z~1.5-2
    alpha = -0.88
    phi_star = 1.22e-3
    logM_star = 10.89
    M_star = 10**logM_star
    x = np.linspace(10.02,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=5,label='Huertas-Company et al. 2016',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=5,alpha=1)

    # Davidzon 17 
    #z: 1.5-2.0
    logM_star = 10.51
    alpha1 = -1.28
    phi_star1 = 0.969e-3
    alpha2 = 0.82
    phi_star2 = 0.64e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)


def plot_obsGSMF_z1pt3(ax,hubble_h=0.678,markersize=10,legend=True):

    # Pozzetti 07
    # K selected complex
    # z~1.2-1.6
    alpha = -1.17
    phi_star = 1.39e-3
    logM_star = 10.93
    M_star = 10**logM_star
    x = np.linspace(9.5,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:                                                      
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',\
                lw=2,label='Pozzetti et al. 2007',alpha=1)          
    else:                                                           
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',lw=2,alpha=1)            
    
    # Mortlock 11
    data = []
    file = '/home/yqin/dragons/data/smf_Mortlock2011_z1pt25.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    if legend:                                                                              
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray',label='Mortlock et al. 2011')
    else:                                                                                   
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray')                             
    
    # Santini 12 
    #z:1.0-1.4
    alpha = -1.46
    logM_star = 11.50
    logphi_star = -3.49
    M_star = 10**logM_star
    phi_star = 10**logphi_star
    x = np.linspace(8,13)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:                                                                       
        ax.plot(x+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=5,label='Santini et al. 2012',alpha=1)
    else:                                                                            
        ax.plot(x+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=5,alpha=1)                            
    
    # Ilbert 13 
    #z: 1.1-1.5
    alpha1 = -0.08
    phi_star1 = 1.35e-3
    alpha2 = -1.46
    phi_star2 = 0.67e-3
    logM_star = 10.71
    M_star = 10**logM_star
    x = np.linspace(9.5,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',\
                lw=2,label='Ilbert et al. 2013',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=2,alpha=1)
    
    # Muzzin 13 
    #z: 1.-1.5
    alpha = -1.02
    logM_star = 10.87
    M_star = 10**logM_star
    phi_star = 13.91e-4
    x = np.linspace(10.03,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:                                                      
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',\
                lw=5,label='Muzzin et al. 2013',alpha=1)            
    else:                                                           
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',lw=5,alpha=1)           
    
    # Huertas-Company 16
    #z~1.1-1.5
    alpha = -1.2
    phi_star = 0.94e-3
    logM_star = 10.97
    M_star = 10**logM_star
    x = np.linspace(9.6,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:                                                     
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=5,label='Huertas-Company et al. 2016',alpha=1)  
    else:                                                          
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=5,alpha=1)           
    
    # Davidzon 17 
    #z: 1.1-1.5
    logM_star = 10.62
    alpha1 = -1.28
    phi_star1 = 1.069e-3
    alpha2 = 0.29
    phi_star2 = 1.21e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)

        
def plot_obsGSMF_z0pt95(ax,hubble_h=0.678,markersize=10,legend=True):
    
    # Pozzetti 07
    # K selected complex
    # z~0.9-1.2
    alpha = -1.2
    phi_star = 1.34e-3
    logM_star = 11.07
    M_star = 10**logM_star
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',\
                lw=2,label='Pozzetti et al. 2007',alpha=1)           
    else:                                                            
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',lw=2,alpha=1)             
    
    # Drory 09
    data = []
    file = '/home/yqin/dragons/data/smf_Drory2009.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[39:48],dtype=np.float)
    data[:,2] +=np.log10(Cha2Sal)
    if legend:
        ax.errorbar(data[:,2]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,7]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,8],\
                    markersize=markersize,fmt='x',color='gray',label='Drory et al. 2009')
    else:
        ax.errorbar(data[:,2]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,7]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,8],\
                    markersize=markersize,fmt='x',color='gray')
    
    # Mortlock 11
    data = []
    file = '/home/yqin/dragons/data/smf_Mortlock2011_z1pt25.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray',label='Mortlock et al. 2011') 
    else:                                                                                    
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray')                              
    
    # Ilbert 13 
    #z: 0.8-1.1
    alpha1 = -0.52
    phi_star1 = 2.03e-3
    alpha2 = -1.62
    phi_star2 = 0.29e-3
    logM_star = 10.87
    M_star = 10**logM_star
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',\
                lw=2,label='Ilbert et al. 2013',alpha=1)            
    else:                                                           
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=2,alpha=1)           
    
    # Muzzin 13 
    #z: 0.5-1
    alpha = -1.17
    logM_star = 11.00
    M_star = 10**logM_star
    phi_star = 16.25e-4
    x = np.linspace(10.03,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',\
                lw=5,label='Muzzin et al. 2013',alpha=1)            
    else:                                                           
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',lw=5,alpha=1)           
        
    # Huertas-Company 16
    #z~0.8-1.1
    alpha = -1.25
    phi_star = 1.33e-3
    logM_star = 11.03
    M_star = 10**logM_star
    x = np.linspace(9.2,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:                                                    
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=5,label='Huertas-Company et al. 2016',alpha=1) 
    else:                                                         
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=5,alpha=1)          
        
    # Davidzon 17 
    #z: 0.8-1.1
    logM_star = 10.56
    alpha1 = -1.31
    phi_star1 = 1.428e-3
    alpha2 = 0.51
    phi_star2 = 2.19e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)

def plot_obsGSMF_z0pt55(ax,hubble_h=0.678,markersize=10,legend=True):
        
    # Pozzetti 07
    # K selected complex
    # z~0.4-0.7
    alpha = -1.16
    phi_star = 1.74e-3
    logM_star = 10.98
    M_star = 10**logM_star
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:                                                      
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',\
                lw=2,label='Pozzetti et al. 2007',alpha=1)          
    else:                                                           
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-',lw=2,alpha=1)            
    
    # Drory 09
    data = []
    file = '/home/yqin/dragons/data/smf_Drory2009.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[16:27],dtype=np.float)
    data[:,2] +=np.log10(Cha2Sal)
    if legend:
        ax.errorbar(data[:,2]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,7]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,8],\
                    markersize=markersize,fmt='x',color='gray',label='Drory et al. 2009') 
    else:                                                                                 
        ax.errorbar(data[:,2]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,7]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,8],\
                    markersize=markersize,fmt='x',color='gray')                           
    
    # Ilbert 13 
    #z: 0.5-0.8
    alpha1 = -1.0
    phi_star1 = 1.22e-3
    alpha2 = -1.64
    phi_star2 = 0.16e-3
    logM_star = 11.03
    M_star = 10**logM_star
    x = np.linspace(8.5,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',\
                lw=2,label='Ilbert et al. 2013',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='--',lw=2,alpha=1)

    # Muzzin 13 
    #z: 0.5-1
    alpha = -1.17
    logM_star = 11.00
    M_star = 10**logM_star
    phi_star = 16.25e-4
    x = np.linspace(10.03,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    if legend:
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',\
                lw=5,label='Muzzin et al. 2013',alpha=1)
    else:
        ax.plot(x+np.log10(Kro2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle='-.',lw=5,alpha=1)
    
    # Huertas-Company 16
    #z~0.5-0.8
    alpha1 = -0.82
    phi_star1 = 2.22e-3
    alpha2 = -1.6
    phi_star2 = 0.45e-3
    logM_star = 10.86
    x = np.linspace(8.9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=5,label='Huertas-Company et al. 2016',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=5,alpha=1)

    # Davidzon 17 
    #z: 0.5-0.8
    logM_star = 10.77
    alpha1 = -1.36
    phi_star1 = 1.070e-3
    alpha2 = 0.03
    phi_star2 = 1.68e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    if legend:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',\
                lw=2,label='Davidzon et al. 2017',alpha=1)
    else:
        ax.plot(x+np.log10(Cha2Sal)+_h_convertor_x2(0.7,hubble_h),\
                np.log10(y)+_h_convertor_y(0.7,hubble_h),\
                color='gray',linestyle=':',lw=2,alpha=1)

    

def plot_obsGSMF_z0pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    # Drory 09
    data = []
    file = '/home/yqin/dragons/data/smf_Drory2009.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[2:15],dtype=np.float)
    data[:,2] +=np.log10(Cha2Sal)
    if legend:                                                                            
        ax.errorbar(data[:,2]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,7]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,8],\
                    markersize=markersize,fmt='x',color='gray',label='Drory et al. 2009') 
    else:                                                                                 
        ax.errorbar(data[:,2]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,7]+_h_convertor_y(0.7,hubble_h),\
                    yerr=data[:,8],\
                    markersize=markersize,fmt='x',color='gray')                           
    
    # Cole 01
    data = []
    file = '/home/yqin/dragons/data/smf_Cole2001.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[5:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(1.0,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(1.0,hubble_h),\
                    yerr=np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                    markersize=markersize,fmt='v',color='black',mfc='None',label='Cole et al. 2001')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(1.0,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(1.0,hubble_h),\
                    yerr=np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                    markersize=markersize,fmt='v',color='black',mfc='None')
    
    # Bell 03
    data = []
    file = '/home/yqin/dragons/data/smf_Bell2003.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[6:],dtype=np.float)
    data[:,0] -= np.log10(0.7)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(1.0,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(1.0,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,2]),\
                          np.log10(data[:,3])-np.log10(data[:,1])],
                    markersize=markersize,fmt='^',color='black',mfc='None',label='Bell et al. 2003')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(1.0,hubble_h),\
                    np.log10(data[:,1])+_h_convertor_y(1.0,hubble_h),\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,2]),\
                          np.log10(data[:,3])-np.log10(data[:,1])],
                    markersize=markersize,fmt='^',color='black',mfc='None')
    
    # Baldry 08
    data = []
    file = '/home/yqin/dragons/data/smf_Baldry2008.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[13:],dtype=np.float)
    data[:,0] -= np.log10(0.7)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,2])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,2])-np.log10(data[:,4]),\
                          np.log10(data[:,5])-np.log10(data[:,2])],
                    markersize=markersize,fmt='s',color='black',mfc='None',label='Baldry et al. 2008')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    np.log10(data[:,2])+_h_convertor_y(0.7,hubble_h),\
                    yerr=[np.log10(data[:,2])-np.log10(data[:,4]),\
                          np.log10(data[:,5])-np.log10(data[:,2])],
                    markersize=markersize,fmt='s',color='black',mfc='None')
    
    # Thanjavur 15
    data = []
    file = '/home/yqin/dragons/data/smf_Thanjavur15.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[2:],dtype=np.float)
    data[:,0] +=np.log10(Cha2Sal)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    yerr=[data[:,3],data[:,2]],\
                    markersize=markersize,fmt='v',color='black',mfc='None',label='Thanjavur et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x2(0.7,hubble_h),\
                    data[:,1]+_h_convertor_y(0.7,hubble_h),\
                    markersize=markersize,fmt='v',color='black',mfc='None')


def plot_obsGSMF_legend(ax,markersize=10):
        
    #Gonzalez 11
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Gonzalez11_z7_smf.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    ax.errorbar(data[:,0],data[:,1],yerr=[data[:,2],data[:,3]],\
                markersize=markersize,fmt='s',color='gray',label='Gonzalez et al. 2011')
    
    #Duncan 14
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Duncan14_MF_z7.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[8:],dtype=np.float)
    ax.errorbar(data[:,0]+np.log10(Cha2Sal),np.log10(data[:,1]),\
                yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,2]),\
                      np.log10(data[:,1]+data[:,3])-np.log10(data[:,1])],\
                markersize=markersize,fmt='o',color='gray',label='Duncan et al. 2014')
    
    #Song 15
    data = []
    file = '/home/yqin/bitbucket/meraxes/meraxes/utils/obs_datasets/smf/Song2015_z7_smf.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[19:],dtype=np.float)
    ax.errorbar(data[:,0],data[:,1],yerr=[data[:,2],data[:,3]],\
                markersize=markersize,fmt='D',color='gray',label='Song et al. 2015')
    
    #Grazian 15
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFgrazian2015_z7.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0],np.log10(data[:,1]),\
                yerr=[np.log10(data[:,1])-np.log10(data[:,3]),\
                      np.log10(data[:,2])-np.log10(data[:,1])],\
                markersize=markersize,fmt='^',color='gray',label='Grazian et al. 2015')
    
    #Marchesini 09
    data = []
    file = '/home/yqin/dragons/data/smf_Marchesini2009_z2pt5.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    data[:,2] +=np.log10(Kro2Sal)
    ax.errorbar(data[:,0],data[:,1],yerr=[data[:,3],data[:,2]],\
                    markersize=markersize,fmt='*',color='gray',label='Marchesini et al. 2009')
    
    # Mortlock 11
    data = []
    file = '/home/yqin/dragons/data/smf_Mortlock2011_z2pt25.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    ax.errorbar(data[:,0],data[:,1],yerr=data[:,2],\
                    markersize=markersize,fmt='.',color='gray',label='Mortlock et al. 2011')
    
    #Khostovan 16
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFkhostovan16_z2.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+np.log10(Cha2Sal),data[:,1],\
                 yerr=[data[:,2],data[:,3]],\
                 markersize=markersize,fmt='>',color='gray',label='Khostovan16_HbetaOIII')
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFkhostovan16OII_z2.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+np.log10(Cha2Sal),data[:,1],\
                 yerr=[data[:,2],data[:,3]],\
                 markersize=markersize,fmt='<',color='gray',label='Khostovan16_OII')
    
    # Pozzetti 07
    # K selected complex
    # z~1.6-2.5
    alpha = -1.17
    phi_star = 1.25e-3
    logM_star = 10.97
    M_star = 10**logM_star
    x = np.linspace(10.0,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    ax.plot(x+np.log10(Cha2Sal),np.log10(y),color='gray',linestyle='-',\
            lw=2,label='Pozzetti et al. 2007',alpha=1)
    
    # Santini 12 #z:1.8-2.5
    alpha = -1.58
    logM_star = 11.29
    logphi_star = -3.52
    M_star = 10**logM_star
    phi_star = 10**logphi_star
    x = np.linspace(8,13)
    y = _Schechter(x, logM_star, alpha, phi_star)
    ax.plot(x,np.log10(y),color='gray',linestyle='--',lw=5,label='Santini et al. 2012',alpha=1)
    
    # Ilbert 13 #z: 1.5-2
    alpha1 = -0.24
    phi_star1 = 0.88e-3
    alpha2 = -1.6
    phi_star2 = 0.33e-3
    logM_star = 10.74
    M_star = 10**logM_star
    x = np.linspace(9.67,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    ax.plot(x+np.log10(Cha2Sal),np.log10(y),color='gray',linestyle='--',\
            lw=2,label='Ilbert et al. 2013',alpha=1)
    
    # Muzzin 13 #z: 1.5-2
    alpha = -0.86
    logM_star = 10.81
    M_star = 10**logM_star
    phi_star = 10.13e-4
    x = np.linspace(10.03,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    ax.plot(x+np.log10(Kro2Sal),np.log10(y),color='gray',linestyle='-.',\
            lw=5,label='Muzzin et al. 2013',alpha=1)
    
    # Tomczak 14 #z:1.5-2.0
    data = []
    file = '/home/yqin/dragons/programs/hydro/MFtomczak_z1pt75.dat'
    for line in open(file,'r'): data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+np.log10(Cha2Sal),data[:,1],\
                 yerr=[data[:,2],data[:,3]],\
                 markersize=markersize,fmt='v',color='gray',label='Tomczak et al. 2014')
    
    # Huertas-Company 16 #z:1.5-2
    alpha = -0.88
    phi_star = 1.22e-3
    logM_star = 10.89
    M_star = 10**logM_star
    x = np.linspace(10.02,12)
    y = _Schechter(x, logM_star, alpha, phi_star)
    ax.plot(x+np.log10(Cha2Sal),np.log10(y),color='gray',linestyle=':',\
                lw=5,label='Huertas-Company et al. 2016',alpha=1)
    
    # Drory 09
    data = []
    file = '/home/yqin/dragons/data/smf_Drory2009.txt'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[39:48],dtype=np.float)
    data[:,2] +=np.log10(Cha2Sal)
    ax.errorbar(data[:,2],data[:,7],yerr=data[:,8],\
                markersize=markersize,fmt='x',color='gray',label='Drory et al. 2009')

    #Stefanon 16
    data = []
    file = '/home/yqin/dragons/data/smf_Stefanon2016_z5.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[3:],dtype=np.float)
    data[:,0] +=np.log10(Cha2Sal)
    ax.errorbar(data[:,0],np.log10(data[:,1])-5.,\
                    yerr=[np.log10(data[:,1])-np.log10(data[:,1]-data[:,3]),\
                          np.log10(data[:,1]+data[:,2])-np.log10(data[:,1])],\
                    markersize=markersize,fmt='p',color='gray',mfc='None',\
                    label='Stefanon et al. 2016')

    # Davidzon 17 
    logM_star = 10.77
    alpha1 = -1.36
    phi_star1 = 1.070e-3
    alpha2 = 0.03
    phi_star2 = 1.68e-3
    x = np.linspace(9,12)
    y = _Schechter(x, logM_star, alpha1, phi_star1, alpha2, phi_star2)
    ax.plot(x+np.log10(Cha2Sal),np.log10(y),color='gray',linestyle=':',\
            lw=2,label='Davidzon et al. 2017',alpha=1)

def plot_obsGSMF(ax,z,hubble_h=0.678,markersize=10,legend=True,silent=False):
    _plot_obsGSMF = {7.0:  plot_obsGSMF_z7pt0,\
                     5.0:  plot_obsGSMF_z5pt0,\
                     4.0:  plot_obsGSMF_z4pt0,\
                     2.0:  plot_obsGSMF_z2pt0,\
                     1.75: plot_obsGSMF_z1pt75,\
                     1.3:  plot_obsGSMF_z1pt3,\
                     0.95: plot_obsGSMF_z0pt95,\
                     0.55: plot_obsGSMF_z0pt55,\
                     0.0:  plot_obsGSMF_z0pt0}

    try:
        (_plot_obsGSMF[z])(ax, hubble_h=hubble_h,markersize=markersize,legend=legend)
    except:
        if not silent:
            print 'available GSMF at redshift:', _plot_obsGSMF.keys()
        return -1
