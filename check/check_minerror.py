#! /usr/bin/env python 

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import check_2dpower
from mpl_toolkits.axes_grid1 import ImageGrid

np.set_printoptions(threshold=np.nan)

def plot2dpower_minerror(file_root, file_name, mode_list, kind='2d'):

    power = []
    power_error = []
    power_weight = []

    power_k = np.load(file_root%mode_list[0] + 
                         file_name%mode_list[0] + 
                         '_k2_combined.npy')

    for mode in mode_list:
        root = file_root%mode
        name = file_name%mode

        p2 = np.load(root + name + '_p2_combined.npy')
        dp2= np.load(root + name + '_p2_var_combined.npy')

        b2 = check_2dpower.getbias(name, from_cross=True)

        w2 = check_2dpower.getweight(root, name, 'noise', b = b2)

        sim= np.load(file_root.replace('power', 'reference')%20 + 
                     'simmaps_p2_combined.npy')
        sim[sim==0] = np.inf
        dp2 /= sim
        p2  /= sim

        p2 *= b2
        dp2*= b2

        power.append(p2)
        power_error.append(dp2)
        power_weight.append(w2)

    power = np.array(power)
    power_error = np.array(power_error)
    power_weight = np.array(power_weight)

    min_error = np.min(power_error, axis=0)

    power_modes = np.array(mode_list)[:,None,None] * np.ones(min_error.shape)

    print min_error.shape
    power[power_error!=min_error] = ma.masked
    power_error[power_error!=min_error] = ma.masked
    power_weight[power_error!=min_error] = ma.masked
    power_modes[power==0] = ma.masked
    power_modes[power_error!=min_error] = ma.masked

    power = np.sum(power, axis=0)
    power_error = np.sum(power_error, axis=0)
    power_weight = np.sum(power_weight, axis=0)
    power_modes = np.sum(power_modes, axis=0)

    power[power_error>10] = ma.masked
    power_weight[power_error>10] = ma.masked
    power_modes[power_error>10] = ma.masked
    power_error[power_error>10] = ma.masked

    kp = power_k[0]
    kv = power_k[1]

    p1, modenum, dp, bc, be = check_2dpower.twod2oned(power, power_k, 
                                                      power_error, power_weight)
    #p1, modenum, dp, bc, be, power, power_error, power_weight, kp, kv = \
    #    check_2dpower.twod2oned_reshape(power, power_k, power_error,
    #    k_range = [0.04, 0, 0.08, 0.2], weight=power_weight)
    if kind=='2draw':
        # plot 2d power
        p2 = ma.array(power)
        dp2 = ma.array(power_error)

        k = power_k[0]

        p2_negative = -p2
        dp2 = dp2**2

        p2[p2<=0] = ma.masked
        p2_negative[p2_negative<=0] = ma.masked
        dp2[dp2<=0] = ma.masked

        dp2[dp2==0] = np.inf
        p2 /= dp2
        p2_negative /= dp2

        p2 = np.ma.log10(p2)
        p2_negative = np.ma.log10(p2_negative)
        dp2= np.ma.log10(dp2)

        fig = plt.figure(figsize=(15,9))
        ax = ImageGrid(fig, 111,
                       nrows_ncols = (1, 3),
                       direction = "row",
                       axes_pad = 0.3,
                       add_all = True,
                       label_mode = "L",
                       share_all = True,
                       cbar_location = "right",
                       cbar_mode = "each",
                       cbar_size = "5%",
                       cbar_pad = 0.00,
                       )

        label = file_name.replace("legendre_modes_0gwj_", "")
        label = label.replace("_%d", "")
        cmax = 0 
        cmin = -4

        im0 = ax[0].pcolormesh(kv,kp,p2)
        im0.set_clim(cmin, cmax)
        ax[0].set_xlim(kp.min(), kp.max())
        ax[0].set_ylim(kp.min(), kp.max())
        ax[0].loglog()
        ax[0].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[0].set_ylabel('log(k_p) [log(h/Mpc)]')
        ax[0].set_title(label + '\n2d power \n$log_{10}(P(k)/\delta P(k)^{2})$')

        ax[0].cax.colorbar(im0)

        im1 = ax[1].pcolormesh(kv,kp,p2_negative)
        im1.set_clim(cmin, cmax)
        ax[1].set_xlim(kp.min(), kp.max())
        ax[1].set_ylim(kp.min(), kp.max())
        ax[1].loglog()
        ax[1].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[1].set_title(label + '\n2d power negative \n$log_{10}(-P(k)/\delta P(k)^{2})$')

        ax[1].cax.colorbar(im1)

        im2 = ax[2].pcolormesh(kv,kp,dp2)
        #im2.set_clim(2*cmin, 2*cmax)
        ax[2].set_xlim(kp.min(), kp.max())
        ax[2].set_ylim(kp.min(), kp.max())
        ax[2].loglog()
        ax[2].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[2].set_title(label + '\n2d error \n$log_{10}((\delta P(k)/P_{sim}(k))^2)$')

        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        #plt.title(label)
        plt.savefig('./png/minerror_2d2d_power_%s'%file_name[:-3])


    if kind=='2d':
        # plot 2d power
        p2 = ma.array(power)
        dp2 = ma.array(power_error)

        k = power_k[0]

        p2_negative = -p2
        dp2 = dp2**2

        p2[p2<=0] = ma.masked
        dp2[dp2<=0] = ma.masked

        p2 = np.ma.log10(p2)
        p2_negative = np.ma.log10(p2_negative)
        dp2= np.ma.log10(dp2)

        fig = plt.figure(figsize=(20,9))
        ax = ImageGrid(fig, 111,
                       nrows_ncols = (1, 4),
                       direction = "row",
                       axes_pad = 0.3,
                       add_all = True,
                       label_mode = "L",
                       share_all = True,
                       cbar_location = "right",
                       cbar_mode = "each",
                       cbar_size = "5%",
                       cbar_pad = 0.00,
                       )

        label = file_name.replace("legendre_modes_0gwj_", "")
        label = label.replace("_%d", "")
        cmax = 3
        cmin = -1

        im0 = ax[0].pcolormesh(kv,kp,p2)
        im0.set_clim(cmin, cmax)
        ax[0].set_xlim(kp.min(), kp.max())
        ax[0].set_ylim(kp.min(), kp.max())
        ax[0].loglog()
        ax[0].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[0].set_ylabel('log(k_p) [log(h/Mpc)]')
        ax[0].set_title(label + '\n2d power \n$log_{10}(P(k)/P_{sim}(k))$')

        ax[0].cax.colorbar(im0)

        im1 = ax[1].pcolormesh(kv,kp,p2_negative)
        im1.set_clim(cmin, cmax)
        ax[1].set_xlim(kp.min(), kp.max())
        ax[1].set_ylim(kp.min(), kp.max())
        ax[1].loglog()
        ax[1].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[1].set_title(label + '\n2d power negative \n$log_{10}(-P(k)/P_{sim}(k))$')

        ax[1].cax.colorbar(im1)

        im2 = ax[2].pcolormesh(kv,kp,dp2)
        #im2.set_clim(2*cmin, 2*cmax)
        ax[2].set_xlim(kp.min(), kp.max())
        ax[2].set_ylim(kp.min(), kp.max())
        ax[2].loglog()
        ax[2].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[2].set_title(label + '\n2d error \n$log_{10}((\delta P(k)/P_{sim}(k))^2)$')

        ax[2].cax.colorbar(im2)

        power_modes = ma.array(power_modes)
        power_modes[power_modes<min(mode_list)] = ma.masked
        im3 = ax[3].pcolormesh(kv,kp,power_modes)
        #im2.set_clim(2*cmin, 2*cmax)
        ax[3].set_xlim(kp.min(), kp.max())
        ax[3].set_ylim(kp.min(), kp.max())
        ax[3].loglog()
        ax[3].set_xlabel('log(k_v) [log(h/Mpc)]')
        ax[3].set_title(label + '\nmodes')

        ax[3].cax.colorbar(im3, ticks=mode_list)

        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        #plt.title(label)
        plt.savefig('./png/minerror_2d2d_power_%s'%file_name[:-3])

    if kind == '1d':
        # plot 1d power
        p = p1[p1>0]
        p_negative = -p1[p1<0]

        pe = np.ndarray(shape=(2, p.shape[0]))
        pe[0] = dp[p1>0]
        pe[1] = dp[p1>0]
        pe[0].put(np.where(pe[0]>=p), p[pe[0]>=p]-1.e-12)

        pe_negative = np.ndarray(shape=(2, p_negative.shape[0]))
        pe_negative[0] = dp[p1<0]
        pe_negative[1] = dp[p1<0]
        pe_negative[0].put(np.where(pe_negative[0]>=p_negative), 
                         p_negative[pe_negative[0]>=p_negative]-1.e-12)

        fig = plt.figure(figsize=(8,7))
        ax = fig.add_subplot(111)
        label = file_name.replace("legendre_modes_0gwj_", "")
        label = label.replace("_%d", "")
        ax.errorbar(bc[p1>0], p, pe, fmt='ro', label=label, capsize=4.5, elinewidth=1)
        ax.errorbar(bc[p1<0], p_negative, pe_negative, fmt='rs', mec='r', mfc='None', capsize=4.5, elinewidth=0.5, label=label + ' negative')

        #fileroot = file_root.replace('power', 'reference')%20
        #filename = 'simmaps'
        #check_2dpower.plot1dpower(fileroot, filename, ax, weights='noise', 
        #                          binedgs=True, color='0.6')
        #filename = 'simmaps_beam'
        #check_2dpower.plot1dpower(fileroot, filename, ax, weights='noise', 
        #                          binedgs=True, color='0.3')

        ymin = 1.e-1
        ymax = 1.e4
        plt.ylim(ymin=ymin, ymax=ymax)
        plt.xlim(xmin=0.025, xmax=1.5)
        plt.xlabel('k [h/Mpc]')
        plt.ylabel('$\Delta^2$ $P(k)/P_{sim}(k)$')
        plt.loglog()
        plt.legend(loc=2, scatterpoints=1, frameon=False)
        
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)


        plt.loglog()
        plt.savefig('./png/minerror_2d1d_power_%s'%file_name[:-3])

def plot1dpower_minerror(file_root, file_name, mode_list):

    power = []
    power_error_upper= []
    power_error_lower= []
    modenum = []
    bin_cent = []
    bin_edgs = []

    for mode in mode_list:
        p, mn, pe, bc, be =\
            check_2dpower.get1dpower(file_root%mode, file_name%mode, weights='noise')
        power.append(p)
        power_error_upper.append(pe[1,:])
        power_error_lower.append(pe[0,:])
        modenum.append(mn)
        bin_cent.append(bc)
        bin_edgs.append(be)


    power = ma.array(power)
    power_error_upper = ma.array(power_error_upper)
    power_error_lower = ma.array(power_error_lower)
    modenum = ma.array(modenum)
    bin_cent = ma.array(bin_cent)

    power[power<=0] = ma.masked
    power_error_upper[power<=0] = ma.masked
    power_error_lower[power<=0] = ma.masked
    modenum[power<=0] = ma.masked
    bin_cent[power<=0] = ma.masked

    #min_error = power_error_upper.min(axis=0)
    min_error = ma.min(power_error_upper, axis=0)

    print min_error
    power_error_upper[power_error_upper!=min_error] = ma.masked
    power_error_lower[power_error_upper!=min_error] = ma.masked
    power[power_error_upper!=min_error] = ma.masked
    modenum[power_error_upper!=min_error] = ma.masked
    bin_cent[power_error_upper!=min_error] = ma.masked

    fig = plt.figure(figsize=(8,7))
    ax = fig.add_subplot(111)
    label = file_name.replace("legendre_modes_0gwj_", "")
    for i in range(len(mode_list)):
        if power.mask[i].all(): continue
        p = np.array(power[i][~power.mask[i]])
        pe= np.array([power_error_lower[i][~power.mask[i]], 
                      power_error_upper[i][~power.mask[i]]])
        bc= np.array(bin_cent[i][~power.mask[i]])

        #p[p==0] = ma.masked
        #pe[0][p==0] = ma.masked
        #pe[1][p==0] = ma.masked
        #bc[p==0] = ma.masked
        
        ax.errorbar(bc, p, pe, fmt='o', label=label%mode_list[i], 
                    capsize=4.5, elinewidth=1, )

    #fileroot = work_root + 'reference_auto_1hr_IE_legendre_modes_0gwj_14conv_20/'
    fileroot = file_root.replace('power', 'reference')%20
    filename = 'simmaps'
    check_2dpower.plot1dpower(fileroot, filename, ax, weights='noise', 
                              binedgs=True, color='0.6')
    filename = 'simmaps_beam'
    check_2dpower.plot1dpower(fileroot, filename, ax, weights='noise', 
                              binedgs=True, color='0.3')

    plt.ylim(ymin=ymin, ymax=ymax)
    plt.xlim(xmin=0.025, xmax=1.5)
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('$\Delta^2$ [$K^2$]')
    plt.loglog()
    plt.legend(loc=2, scatterpoints=1, frameon=False)
    
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)


    plt.loglog()
    plt.savefig('./png/minerror_1d_power_%s'%file_name[:-3])

    #print power_error_upper
    #print power_error_lower
    #print modenum


if __name__=="__main__":

    ymin = 3.e-11
    ymax = 3.e-4

    conv = '14conv'
    hour = '1hr'

    mode_list = [10, 15, 20, 25, 30, 35, 40]
    
    work_root = "/mnt/raid-project/gmrt/ycli/ps_result/"

    file_root = work_root + "power_auto_" + hour +\
                "_IE_legendre_modes_0gwj_" + conv + "_%d/"
    file_name = "auto_" + hour + "_IE_legendre_modes_0gwj_" + conv + "_%d"

    #plot1dpower_minerror(file_root, file_name, mode_list)
    #plot2dpower_minerror(file_root, file_name, mode_list, kind='2draw')
    plot2dpower_minerror(file_root, file_name, mode_list, kind='2d')
    plot2dpower_minerror(file_root, file_name, mode_list, kind='1d')

