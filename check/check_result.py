#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import math

legendtitle = ''
color = ['k', 'g', 'r', 'b', 'y', 'c']
workroot = os.getenv('YICHAO_WORK_PS')

#FIGNAME = 'auto_power_%s'%maps
FIGNAME = 'auto_power_chechreslt'
ymax = 3.e-4
ymin = 3.e-12

def getpower(p_root, k_root, v_root=None):
    p = np.load(p_root)
    k = np.load(k_root)
    if v_root!=None:
        v = np.load(v_root)

    dk = math.sqrt(k[1]/k[0])

    p[p<0] = 0
    non0 = p.nonzero()

    p = p.take(non0)[0]
    k = k.take(non0)[0]

#   print p
#   print k
#   print non0

    kl= k
    k = k*dk
    ku= k*dk
    ke= np.ndarray(shape=(2, len(non0[0])))
    ke[0]= k-kl
    ke[1]= ku-k

    if v_root!=None:
        e = np.ndarray(shape=(2, len(non0[0])))
        e[0] = v.take(non0)[0]
        e[1] = v.take(non0)[0]
        for ii in range(len(e[0])):
            if e[0][ii] >= p[ii]:
                e[0][ii] =  p[ii]-1.e-15

        return k, p, e, ke
    else:
        return k, p

def plotsim(modes):
    
    p_root = workroot + '%s_%d/simmaps_p_combined.npy'%(maps_sim, modes)
    k_root = workroot + '%s_%d/simmaps_k_combined.npy'%(maps_sim, modes)
    v_root = workroot + '%s_%d/simmaps_p_var_combined.npy'%(maps_sim, modes)

    k, p, e, ke = getpower(p_root, k_root, v_root)
    for i in range(len(p)):
        print "%e\t%e\t%e\t\n"%(k[i], p[i], e[0][i])
    plt.errorbar(k, p, e, ke, fmt='o', c='0.5',
        label='simulations 512*128*64',
        capsize=4.5, elinewidth=2)

def plotth(resultroot):
    p_root = resultroot[0] + 'pk_camb.npy'
    p = np.load(p_root)

    #OmegaHI = 0.4 * 1.e-3
    #OmegaHI = 1.e-3
    OmegaHI = 1.e-3
    Omegam = 0.27
    OmegaL = 0.73
    z = 0.8
    def get_distortion(b):
        f = (Omegam*(1+z)**3)**0.55
        t = (1.+(2./3.)*(f/b)+(1./5.)*(f/b)**2)
        return t

    b_opt = 1.2
    t_opt = get_distortion(b_opt)
    t_opt = 1
    b_gbt = 1.35
    t_gbt = get_distortion(b_gbt)
    t_gbt = 1.

    a3 = (1+z)**(-3)
    Tb = 0.39*(OmegaHI)*((Omegam + a3*OmegaL)/0.29)**(-0.5)*((1.+z)/2.5)**0.5

    p[1] = p[1]*Tb**2
    #p[1] = p[1]*b_gbt*b_gbt*t_gbt
    p[1] = p[1]*p[0]*p[0]*p[0]/2./3.1415926/3.1415926

    plt.plot(p[0], p[1], label="theoretical power spectrum", 
        c='0.65', linewidth=3)

def plot(resultroot, resultfile, label=None, combined=True):
    if label == None:
        label = resultfile
    plt.figure(figsize=(8,7))

    #plotsim(15)
    #plotth()
    
    for i in range(len(resultfile)):
        if combined:
            p_root = resultroot[i] + resultfile[i] + '_p_combined.npy'
            k_root = resultroot[i] + resultfile[i] + '_k_combined.npy'
            v_root = resultroot[i] + resultfile[i] + '_p_var_combined.npy'
            k, p, e, ke = getpower(p_root, k_root, v_root)

            plt.errorbar(k, p, e, ke, fmt='o', c=color[i],
                label=resultfile[i],
                capsize=4.5, elinewidth=1)
        else:
            p_root = resultroot[i] + resultfile[i] + '_p.npy'
            k_root = resultroot[i] + resultfile[i] + '_k.npy'
            k, p = getpower(p_root, k_root)

            plt.scatter(k, p, marker='o', c=color[i], label=label[i], 
                s=40, alpha=0.5)

    
    plt.loglog()
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.xlim(xmin=0.025, xmax=1.5)
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('$\Delta^2$ [$K^2$]')
    #plt.ylabel('$\Delta$ [$K$]')
    plt.legend(loc=0, scatterpoints=1, frameon=False, title=legendtitle)
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    
    plt.savefig('./png/'+FIGNAME+'.png', format='png')
    #plt.savefig('./eps/'+FIGNAME+'.eps', format='eps')
    #plt.savefig(FIGNAME+'.png', format='png', dpi=130)
    #plt.savefig(FIGNAME+'.eps', format='eps', dpi=130)
    #plt.show()

if __name__=="__main__":
    resultroot = [
        #'/mnt/raid-project/gmrt/ycli/ps_result/bias/auto_GBT_15hr_map_oldcal_legendre_modes_5gwj_15/',
        '/mnt/raid-project/gmrt/ycli/ps_result/bias/auto_GBT_15hr_map_oldcal_legendre_modes_0gwj_15/',
        '/mnt/raid-project/gmrt/ycli/ps_result/bias/auto_GBT_15hr_map_oldcal_legendre_modes_0gwj_15/',
        '/mnt/raid-project/gmrt/ycli/ps_result/reference_auto_sim_15hr_oldmap_ideal_15/',
        #'/mnt/raid-project/gmrt/ycli/ps_result/bias/auto_GBT_15hr_map_oldcal_legendre_modes_5gwj_15/',
        ]

    resultfile = ['simmaps1','simmaps0','simmaps']

    label = ['P[clean_{sim+map}(sim+map)]', 'P[clean_{map}(map)]', 'P[sim]']
    
    plot(resultroot, resultfile, label=label, combined=False)
