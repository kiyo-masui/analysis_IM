#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import math


def plot2dpower(fileroot, filename, savename, cros=False):
    print 'using map:'
    print fileroot+filename
    if cros:
        p2 = np.load(fileroot + filename + '_p2.npy')
    else:
        p2 = np.load(fileroot + filename + '_p2_combined.npy')
    dp2 = np.load(fileroot + filename + '_p2_var_combined.npy')
    k2 = np.load(fileroot + filename + '_k2_combined.npy')

    #p, modenum, dp, bin_cent, bin_edgs = twod2oned(p2, k2, dp2)

    p2 = np.log10(p2)
    k2 = np.log10(k2)


    plt.figure(figsize=(10,8))
    plt.imshow(p2, origin='lower', extent=(k2[0][0], k2[0][-1], k2[1][0], k2[1][-1]),
        interpolation='none', aspect=1)
    plt.colorbar(aspect=30).set_label('2D $\Delta^2$ [$K^2$]')
    plt.xlabel('log(k_v) [log(h/Mpc)]')
    plt.ylabel('log(k_p) [log(h/Mpc)]')
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.title(filename)
    plt.savefig('./png/'+savename+'_2d.png', format='png')


def plot1dpower(fileroot, filename, ax1, ax2=None, savename=None, weights=None, cros=False):
    print 'using map:'
    print fileroot+filename
    if cros:
        p2 = np.load(fileroot + filename + '_p2.npy')
    else:
        p2 = np.load(fileroot + filename + '_p2_combined.npy')
    dp2 = np.load(fileroot + filename + '_p2_var_combined.npy')
    k2 = np.load(fileroot + filename + '_k2_combined.npy')
 
    if os.path.exists('/mnt/raid-project/gmrt/ycli/ps_result/bias/' + filename):
        b = np.load('/mnt/raid-project/gmrt/ycli/ps_result/bias/'
                    + filename + '/b2_bias.npy')
        if weights=='count':
            kn = np.load(fileroot + filename + '_n2_combined.npy')
            weight = getcountweight(kn, b)
        elif weights=='noise':
            pn = np.load(fileroot + 
                         filename.replace('auto', 'nois') + 
                         '_p2_var_combined.npy')
            weight = getnoiseweight(pn, b)
        else:
            weight = np.ones(p2.shape)
    else:
        print 'no bias'
        if weights=='count':
            kn = np.load(fileroot + filename + '_n2_combined.npy')
            weight = getcountweight(kn)
        elif weights=='noise':
            pn = np.load(fileroot + 
                         filename.replace('auto', 'nois') + 
                         '_p2_var_combined.npy')
            weight = getnoiseweight(pn)
        else:
            weight = np.ones(p2.shape)


    p, modenum, dp, bin_cent, bin_edgs = twod2oned(p2, k2, dp2, weight)

    if savename!=None:
        plt.figure(figsize=(8,7))
    ke = np.ndarray(shape=k2.shape)
    ke[0] = bin_cent - bin_edgs[:-1]
    ke[1] = bin_edgs[1:] - bin_cent
    pe = np.ndarray(shape=k2.shape)
    pe[0] = dp
    pe[1] = dp
    pe[0].put(np.where(pe[0]>=p), p[pe[0]>=p]-1.e-15)
    if not weights==None:
        label = filename + ' ' + weights
    else:
        label = filename
    #f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    #plt.errorbar(bin_cent, p, pe, ke, fmt='o',  
    ax1.errorbar(bin_cent, p, pe, ke, fmt='o',  
        label=label , capsize=4.5, elinewidth=1)
    if not ax2 is None:
        #p[p==0] = np.inf
        ax2.plot(bin_cent, pe[1], 'o', label=label)

    if savename!=None:
        plt.loglog()
        plt.ylim(ymin=3.e-12, ymax=3.e-4)
        plt.xlim(xmin=0.025, xmax=1.5)
        plt.xlabel('k [h/Mpc]')
        plt.ylabel('$\Delta^2$ [$K^2$]')
        #plt.ylabel('$\Delta$ [$K$]')
        plt.legend(loc=0, scatterpoints=1, frameon=False )
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        
        plt.savefig('./png/'+savename+'_1d.png', format='png')

def plotth(fileroot, maptype, ax):
    p_root = fileroot + 'pk_camb.npy'
    p = np.load(p_root)

    #OmegaHI = 0.4 * 1.e-3
    OmegaHI = 0.7*1.e-3
    #OmegaHI = 1.e-3
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

    p[1] = p[1]*Tb
    if maptype == 'auto':
        p[1] = p[1]*Tb
    #p[1] = p[1]*b_gbt*b_gbt*t_gbt
    p[1] = p[1]*p[0]*p[0]*p[0]/2./3.1415926/3.1415926

    #plt.plot(p[0], p[1], label="theoretical power spectrum", 
    ax.plot(p[0], p[1], label="theoretical power spectrum", 
        c='0.65', linewidth=3)

def twod2oned(p2, k2, dp2, weight=None):
    goodlist = np.isnan(p2).__invert__()
    if not weight==None:
        normal = weight[goodlist].sum()
        p2 = p2*weight#/normal
        dp2 = dp2*weight#/normal
    kp= k2[0].reshape([k2.shape[1], 1])
    kv= k2[1]
    k = np.zeros(shape=(2,)+p2.shape)
    k[0,:,:] = kp[:]
    k[1,:,:] = kv[:]
    kmode = np.sqrt(k[0]*k[0] + k[1]*k[1])
    dk = k2[0][1]/k2[0][0]
    bins = k2[0].tolist() +  [k2[0][-1]*dk,]
    modenum, bin_edges = np.histogram(kmode[goodlist], 
                                      bins=bins, 
                                      weights=weight[goodlist])
    p, bin_edges = np.histogram(kmode[goodlist], bins=bins, weights=p2[goodlist])
    dp2 = dp2**2
    dp, bin_edges = np.histogram(kmode[goodlist], bins=bins, weights=dp2[goodlist])
    dp = np.sqrt(dp)
    modenum[modenum==0] = np.inf
    p = p/modenum
    dp= dp/modenum
    bin_cent = bin_edges[:-1]*math.sqrt(dk)
    return p, modenum, dp, bin_cent, bin_edges

def getcountweight(count, transfer_f=None):
    if not transfer_f==None:
        transfer_f[transfer_f==0] = np.inf
        weight = count/transfer_f/transfer_f
    else:
        weight = count
    weight[np.isnan(weight)] = 0.
    return weight
def getnoiseweight(noise, transfer_f=None):
    if not transfer_f==None:
        noise = noise*transfer_f
    noise[noise==0] = np.inf
    weight = 1./noise/noise
    return weight


if __name__=='__main__':

    power_2d = False
    power_1d = True

    workroot = '/mnt/raid-project/gmrt/ycli/ps_result/'

    #maplist = ['cros_GBT_1hr_map_oldcal'] #['cros_GBT_15hr_map_oldcal']
    #maplist = ['cros_IQ', 'cros_IU', 'cros_GBT_15hr_41-90_fdgp_RM']
    maplist = ['cros_II', 'cros_II_extend']
    #maplist = ['cros_IU',]# 'cros_GBT_15hr_map_oldcal']
    maptype = 'cros'#'auto' 
    modelist = [20, ]
    svdnlist = [0, ]

    weightlist = ['count',] #[None, 'count', 'noise']

    #savename = 'auto_15hr_oldmap_str_25gwj_80'
    #savename = 'cros_15hr_oldmap_str_0gwj_20'
    #savename = 'cros_1hr_oldmap_str_0gwj_50'
    #savename = 'cros_15hr_IQ_first_0gwj_1gwj3'
    savename = 'cros_15hr_II_extend_0gwj_20'
    #savename = 'cros_15hr_IU_II_0gwj_20'

    if power_2d:
        #fileroot = workroot + 'simulation_auto_sim_15hr_oldmap_str_50/'
        #filename = 'simmaps'
        #savename = filename + '_2dpower'
        #plot2dpower(fileroot, filename, savename)

        #fileroot = workroot + 'simulation_cros_sim_15hr_oldmap_ideal_20/'
        #filename = 'simmaps'
        #savename = filename + '_2dpower'
        #plot2dpower(fileroot, filename, savename)

        #fileroot = workroot + 'simulation_cros_sim_1hr_oldmap_str_50/'
        #filename = 'simmaps'
        #savename = filename + '_2dpower'
        #plot2dpower(fileroot, filename, savename)

        #exit()

        for map in maplist:
            for mode in modelist:
                for svdn in svdnlist:
                    fileroot = workroot
                    fileroot += 'power_%s_legendre_modes_%dgwj_%d/'%(map, svdn, mode)
                    filename = '%s_legendre_modes_%dgwj_%d'%(map, svdn, mode)
                    #fileroot += 'power_%s_cleaned_%d/'%(map, mode)
                    #filename = '%s_cleaned_%d'%(map, mode)
                    if not os.path.exists(fileroot):
                        print 'file not found: %s%s'%(fileroot, filename)
                        continue
                    #savename = filename + '_2dpower'
                    plot2dpower(fileroot, filename, savename, cros=True)


    if power_1d:
        #plt.figure(figsize=(10,7))
        #plt.subplots_adjust(hspace=0.0)
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        f.set_figheight(10)
        f.set_figwidth(8)
        f.subplots_adjust(hspace=0.0)
        
        for map in maplist:
            for mode in modelist:
                for svdn in svdnlist:
                    fileroot = workroot
                    fileroot += 'power_%s_legendre_modes_%dgwj_%d/'\
                                 %(map,svdn,mode)
                    filename = '%s_legendre_modes_%dgwj_%d'%(map, svdn, mode)
                    
                    #fileroot += 'power_%s_first_legendre_modes_%dgwj_%d/'\
                    #            %(map, svdn, mode)
                    #filename = '%s_first_legendre_modes_%dgwj_%d'%(map, svdn, mode)

                    #fileroot += 'power_%s_legendre_modes_%dgwj_%d/'%(map,svdn,mode)
                    #filename = '%s_legendre_modes_%dgwj_%d'%(map, svdn, mode)

                    #fileroot += 'power_%s_extend_legendre_modes_%dgwj_%d/'\
                    #             %(map,svdn,mode)
                    #filename = '%s_extend_legendre_modes_%dgwj_%d'%(map, svdn, mode)

                    #fileroot += 'power_%s_legendre_modes_%dgwj_conv_%d/'\
                    #             %(map,svdn,mode)
                    #filename = '%s_legendre_modes_%dgwj_conv_%d'%(map, svdn, mode)

                    #fileroot += 'power_%s_cleaned_%d/'%(map, mode)
                    #filename = '%s_cleaned_%d'%(map, mode)
                    if not os.path.exists(fileroot): 
                        print 'No such file %s' %(fileroot)
                        continue
                    for weights in weightlist:
                        plot1dpower(fileroot, filename, ax1, ax2, 
                                    weights=weights, cros=True)
                        #plot1dpower(fileroot, filename, weights=weights, cros=True)
                        #plot1dpower(fileroot, filename, weights=weights)

        if maptype == 'auto':
            # -- compare weight --
            #fileroot = workroot + 'simulation_auto_sim_15hr_oldmap_str_80/'
            #filename = 'simmaps'
            #plot1dpower(fileroot, filename)

            #fileroot = workroot + 'simulation_auto_sim_15hr_oldmap_str_80/'
            #filename = 'simmaps'
            #plot1dpower(fileroot, filename, weights='count')

            fileroot = workroot + 'simulation_auto_sim_15hr_oldmap_str_25/'
            filename = 'simmaps'
            plot1dpower(fileroot, filename, ax1, ax2, weights='noise')

            plotth(fileroot, maptype, ax1)

            plt.ylim(ymin=3.e-12, ymax=3.e-4)
            plt.xlim(xmin=0.025, xmax=1.5)
            plt.xlabel('k [h/Mpc]')
            plt.ylabel('$\Delta^2$ [$K^2$]')

        if maptype == 'cros':
            # -- compare weight --
            fileroot = workroot + 'simulation_cros_15hr_20/'
            filename = 'simmaps'
            #plot1dpower(fileroot, filename, ax1, )

            #fileroot = workroot + 'simulation_cros_sim_15hr_oldmap_str_20/'
            #filename = 'simmaps'
            #plot1dpower(fileroot, filename)

            #fileroot = workroot + 'simulation_cros_sim_15hr_oldmap_str_20/'
            #filename = 'simmaps'
            #plot1dpower(fileroot, filename, weights='count')

            #fileroot = workroot + 'simulation_cros_sim_15hr_oldmap_str_80/'
            #filename = 'simmaps'
            #plot1dpower(fileroot, filename, weights='noise')

            #fileroot = workroot + 'simulation_cros_sim_1hr_oldmap_str_50/'
            #filename = 'simmaps'
            #plot1dpower(fileroot, filename, weights='count')

            plotth(fileroot, maptype, ax1)

            ax1.set_xlim(xmin=0.025, xmax=1.5)
            #ax1.set_ylim(ymin=3.e-9, ymax=3.e3)
            ax1.set_ylim(ymin=3.e-9, ymax=3.e-1)
            ax1.set_ylabel('$\Delta^2$ [$K$]')
            ax1.loglog()
            ax1.legend(loc=0, scatterpoints=1, frameon=False)
            #ax1.ylim(ymin=3.e-8, ymax=3.e-1)
            #ax1.xlim(xmin=0.025, xmax=1.5)
            #ax1.xlabel('k [h/Mpc]')
            #ax1.ylabel('$\Delta^2$ [$K$]')

            ax2.set_xlim(xmin=0.025, xmax=1.5)
            #ax2.set_ylim(ymin=3.e-8, ymax=3.e-1)
            ax2.set_ylabel('Error')
            ax2.legend(loc=0, scatterpoints=1, frameon=False)
            ax2.loglog()
            #ax2.ylim(ymin=3.e-8, ymax=3.e-1)
            #ax2.xlim(xmin=0.025, xmax=1.5)
            #ax2.xlabel('k [h/Mpc]')
            #ax2.ylabel('Error')

            #plt.ylim(ymin=3.e-8, ymax=3.e-1)
            #plt.xlim(xmin=0.025, xmax=1.5)
            plt.xlabel('k [h/Mpc]')
            #3plt.ylabel('$\Delta^2$ [$K$]')

        savename = savename

        #plt.loglog()
        #plt.ylabel('$\Delta$ [$K$]')
        #plt.legend(loc=0, scatterpoints=1, frameon=False )
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        
        plt.savefig('./png/'+savename+'_1d.png', format='png')
