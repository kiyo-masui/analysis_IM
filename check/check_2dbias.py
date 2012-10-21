#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import math

def plot2dtransfer_3(fileroot, filename, savename):
    k  = np.load(fileroot + filename + '/k_bias.npy')

    b2 = np.load(fileroot + filename + '/b2_bias.npy')
    b2[b2==0.] = np.inf
    b2 = 1./b2
    b2[np.isnan(b2)] = 0.
    #b2 = np.ma.masked_equal(b2, np.inf)
    b2 = np.ma.masked_equal(b2, 0)

    b2_beam = np.load(fileroot + filename + '/b2_beam_bias.npy')
    b2_beam[b2_beam==0.] = np.inf
    b2_beam = 1./b2_beam
    b2_beam[np.isnan(b2_beam)] = 0.
    #b2_beam = np.ma.masked_equal(b2_beam, np.inf)
    b2_beam = np.ma.masked_equal(b2_beam, 0.)

    b2_lose = np.load(fileroot + filename + '/b2_lose_bias.npy')
    b2_lose[b2_lose==0.] = np.inf
    b2_lose = 1./b2_lose
    b2_lose[np.isnan(b2_lose)] = 0.
    #b2_lose = np.ma.masked_equal(b2_lose, np.inf)
    b2_lose = np.ma.masked_equal(b2_lose, 0.)

    #print b2
    #print b2_beam
    #print b2_lose

    cmax = 1.2

    f = plt.figure(figsize=(27,11))
    #f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(23, 7))
    #plt.subplots_adjust(wspace=0)
    ax = ImageGrid(f, 111,
                   nrows_ncols = (1, 3),
                   direction = "row",
                   axes_pad = 0.05,
                   add_all = True,
                   label_mode = "L",
                   share_all = True,
                   cbar_location = "right",
                   cbar_mode = "single",
                   cbar_size = "5%",
                   cbar_pad = 0.05,
                   )
    im0 = ax[0].pcolormesh(k, k, b2)
    im0.set_clim(0, cmax)
    ax[0].set_xlim(k.min(), k.max())
    ax[0].set_ylim(k.min(), k.max())
    ax[0].loglog()
    ax[0].set_xlabel('log(k_v) [log(h/Mpc)]')
    ax[0].set_ylabel('log(k_p) [log(h/Mpc)]')
    ax[0].set_title(filename)

    im1 = ax[1].pcolormesh(k, k, b2_beam)
    im1.set_clim(0, cmax)
    ax[1].set_xlim(k.min(), k.max())
    ax[1].loglog()
    ax[1].set_xlabel('log(k_v) [log(h/Mpc)]')
    ax[1].set_title(filename + ' beam')

    im2 = ax[2].pcolormesh(k, k, b2_lose)
    im2.set_clim(0, cmax)
    ax[2].set_xlim(k.min(), k.max())
    ax[2].loglog()
    ax[2].set_xlabel('log(k_v) [log(h/Mpc)]')
    ax[2].set_title(filename + ' lost')

    ax[2].cax.colorbar(im2)

    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.savefig('./png/'+savename + '3.png', format='png')


def plot2dtransfer(fileroot, filename, savename, middle=''):
    b2 = np.load(fileroot + filename + '/b2%s_bias.npy'%middle)
    k  = np.load(fileroot + filename + '/k%s_bias.npy'%middle)

    b2 = 1./b2
    b2 = np.ma.masked_equal(b2, np.inf)
    print b2
    #b2 = np.log10(b2)

    #k = np.log10(k)

    #fig = plt.figure(figsize=(10,8))
    #ax = fig.add_subplot(111)
    #cax = ax.imshow(b2, origin='lower', extent=(k[0], k[-1], k[0], k[-1]),
    #                interpolation='none', aspect=1)
    #cb = fig.colorbar(cax, aspect=30)
    #cb.set_label('2d transfer')
    #cb.set_clim(vmin=0.5, vmax=0.6)
    plt.figure(figsize=(10,8))
    plt.pcolormesh(k, k, b2)
    plt.clim(0., 1.5)
    #plt.clabel('2d transfer')
    plt.colorbar(aspect=30).set_label('2d transfer')
    plt.loglog()
    plt.xlim(k.min(), k.max())
    plt.ylim(k.min(), k.max())
    plt.xlabel('log(k_v) [log(h/Mpc)]')
    plt.ylabel('log(k_p) [log(h/Mpc)]')
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.title(filename)
    plt.savefig('./png/'+savename + middle +'.png', format='png')

if __name__=='__main__':
    workroot = '/mnt/raid-project/gmrt/ycli/ps_result/'

    #fileroot = workroot + 'bias/'
    #filename = 'auto_GBT_15hr_map_oldcal_legendre_modes_60gwj_80'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_GBT_15hr_map_oldcal_legendre_modes_20gwj_50'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_GBT_15hr_map_oldcal_legendre_modes_5gwj_50'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_GBT_15hr_map_oldcal_legendre_modes_40gwj_50'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_GBT_1hr_map_oldcal_legendre_modes_0gwj_50'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IQ_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IU_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_II_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IQ_legendre_modes_0gwj_10'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IU_legendre_modes_0gwj_10'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_II_legendre_modes_0gwj_10'
    #savename = filename + '_transfer'
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IQ_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
    #plot2dtransfer(fileroot, filename, savename, middle='_beam')
    #plot2dtransfer(fileroot, filename, savename, middle='_lose')
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IU_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
    #plot2dtransfer(fileroot, filename, savename, middle='_beam')
    #plot2dtransfer(fileroot, filename, savename, middle='_lose')
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_II_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
    #plot2dtransfer(fileroot, filename, savename, middle='_beam')
    #plot2dtransfer(fileroot, filename, savename, middle='_lose')
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_II_extend_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_1hr_II_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_1hr_IQ_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_1hr_IU_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_1hr_II_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_1hr_IQ_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_1hr_IU_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_15hr_IE_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
     
    #fileroot = workroot + 'bias/'
    #filename = 'cros_15hr_EI_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_15hr_IE_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_15hr_EI_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_1hr_IE_legendre_modes_0gwj_conv_10'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_1hr_IE_legendre_modes_0gwj_conv_5'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_1hr_IE_legendre_modes_0gwj_conv_15'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_1hr_IE_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    fileroot = workroot + 'bias/'
    filename = 'auto_1hr_ABCD_legendre_modes_0gwj_conv_10'
    savename = filename + '_transfer'
    plot2dtransfer_3(fileroot, filename, savename)

    fileroot = workroot + 'bias/'
    filename = 'auto_1hr_ABCD_legendre_modes_0gwj_conv_20'
    savename = filename + '_transfer'
    plot2dtransfer_3(fileroot, filename, savename)

    fileroot = workroot + 'bias/'
    filename = 'auto_1hr_ABCD_legendre_modes_0gwj_conv_25'
    savename = filename + '_transfer'
    plot2dtransfer_3(fileroot, filename, savename)

    fileroot = workroot + 'bias/'
    filename = 'auto_1hr_ABCD_legendre_modes_0gwj_25'
    savename = filename + '_transfer'
    plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'auto_1hr_IE_legendre_modes_0gwj_conv_15'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IQ_first_legendre_modes_0gwj_3'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
    #plot2dtransfer(fileroot, filename, savename, middle='_beam')
    #plot2dtransfer(fileroot, filename, savename, middle='_lose')
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IU_first_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
    #plot2dtransfer(fileroot, filename, savename, middle='_beam')
    #plot2dtransfer(fileroot, filename, savename, middle='_lose')
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_II_first_legendre_modes_0gwj_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
    #plot2dtransfer(fileroot, filename, savename, middle='_beam')
    #plot2dtransfer(fileroot, filename, savename, middle='_lose')
    #plot2dtransfer(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_II_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IQ_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)

    #fileroot = workroot + 'bias/'
    #filename = 'cros_IU_legendre_modes_0gwj_conv_20'
    #savename = filename + '_transfer'
    #plot2dtransfer_3(fileroot, filename, savename)
