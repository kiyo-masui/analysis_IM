#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import math


def plot2dtransfer(fileroot, filename, savename):
    b2 = np.load(fileroot + filename + '/b2_bias.npy')
    k  = np.load(fileroot + filename + '/k_bias.npy')

    b2 = 1./b2
    print b2
    #b2 = np.log10(b2)

    k = np.log10(k)

    plt.figure(figsize=(10,8))
    plt.imshow(b2, origin='lower', extent=(k[0], k[-1], k[0], k[-1]),
        interpolation='none', aspect=1)
    plt.colorbar(aspect=30).set_label('2d transfer')
    plt.xlabel('log(k_v) [log(h/Mpc)]')
    plt.ylabel('log(k_p) [log(h/Mpc)]')
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.title(filename)
    plt.savefig('./png/'+savename+'.png', format='png')

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

    fileroot = workroot + 'bias/'
    filename = 'cros_GBT_1hr_map_oldcal_legendre_modes_0gwj_50'
    savename = filename + '_transfer'
    plot2dtransfer(fileroot, filename, savename)

