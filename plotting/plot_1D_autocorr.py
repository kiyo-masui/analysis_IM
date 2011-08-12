import scipy as sp
import scipy.signal as sig
import scipy.stats.stats as st
import numpy as np
import numpy.ma as ma
import numpy.random as rn
import numpy.fft as fft
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib import cm
import copy
import os
import cPickle
from core import fitsGBT
from core import algebra
from time_stream import flag_data
from time_stream import rotate_pol
#import window_thing as gw
from time_stream import flag_data as fd
from time_stream import cal_scale as cs
from correlate import freq_slices as fs
from correlate.freq_slices import *

def save_1D_corr(mode_number, high, save):
    '''very hardcoded'''
    # Load pickle file.
#    file_name = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest/73_ABCD_all_%d_modes_real3map/New_Slices_object.pkl" % (mode_number)
    file_name = "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/73_ABCD_all_15_modes_realmap_ra_fixed/New_Slices_object.pkl"
    f = open(file_name, "r")
    F = cPickle.load(f)
    f.close()

    # Setting axis info after pickling.
    map_file = F.params["input_root"] + "sec_A_15hr_41-73_clean_map_I.npy"
    exMap = algebra.make_vect(algebra.load(map_file))
    for Pair in F.Pairs:
        Pair.Map1.info = exMap.info
        Pair.Map2.info = exMap.info
        Pair.Noise_inv1.info = exMap.info
        Pair.Noise_inv2.info = exMap.info

    ########## Getting 1D thing. #####################################
    out_list=[]
    d1_list=[]
    for i in range(0,6):
        # The corr to use.
        corr = F.Pairs[i].corr
        # The lags used
        lags = sp.array(F.params['lags'])
        real_lags = copy.deepcopy(lags)
        real_lags[0] = 0
        real_lags[1:] -= sp.diff(lags)/2.0
        # The range selected in ini file.
        frange = F.params['freq']
        # The corresponding real frequencies for that range.
        realrange = [F.Pairs[i].Map1.get_axis('freq')[f] for f in frange]
        # The 2D correlation.
        out = fs.rebin_corr_freq_lag(corr, realrange, nfbins=200, weights=F.Pairs[i].counts,
                                     return_fbins=True)
        out_list.append(out[0])
        #plt.figure()
        #plt.imshow(out[0])
        #plt.colorbar()
        # the 1D correlation.
        d1 = fs.collapse_correlation_1D(out[0], out[2], real_lags, out[1])
        d1_list.append(copy.deepcopy(d1[0]))
    #    plt.figure()
    #    plt.plot(d1[2], d1[0],'.')
        x_axis = d1[2][1]

    matrixx=[]
    for d in d1_list:
        matrixx.append(d.tolist())

    matrixx = sp.array(matrixx)
    print matrixx

    vals=[]
    std=[]
    for i in range(0, matrixx.shape[1]):
        # Get the sqrt to get mK.
        vals.append(sp.mean(sp.sign(matrixx[:,i])*sp.sqrt(abs(matrixx[:,i]))))
        std.append(sp.std(sp.sign(matrixx[:,i])*sp.sqrt(abs(matrixx[:,i]))))

    vals = sp.array(vals)
    std = sp.array(std)

    print
    print
    for i in range(0, matrixx.shape[0]):
        print sp.sign(matrixx[i,:])*sp.sqrt(abs(matrixx[i,:]))

    # Set plot up for log log axes.    
    plt.figure()
    ax = plt.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    # Plot positives.
    t_inds = vals >= 0
    n_inds = vals < 0 
    plt.plot(x_axis[t_inds],vals[t_inds]*1000,'b.') # take out *1000 for simmaps. 
    plt.plot(x_axis[t_inds],(vals[t_inds]+std[t_inds])*1000,'g_')
    plt.plot(x_axis[t_inds],(vals[t_inds]-std[t_inds])*1000,'g_')
    # Plot negatives.
    neg_vals = -1*vals
    plt.plot(x_axis[n_inds],neg_vals[n_inds]*1000,'r.') # take out *1000 for simmaps. 
    plt.plot(x_axis[n_inds],(neg_vals[n_inds]+std[n_inds])*1000,'g_')
    plt.plot(x_axis[n_inds],(neg_vals[n_inds]-std[n_inds])*1000,'g_')
    #plt.axis([1, 100, 0.01, 500.0])
    plt.xlabel('lag (Mpc/h)')
    plt.ylabel('correlation (mK)')

    # Plot the model(?).
    t_lags = sp.arange(0.1, 100, 0.1)
    r0 = 5.5
    rb = 7.0
    t = (sp.sqrt(((rb + t_lags) / r0)**(-1.8)))
    t = t * 0.15 / t[0]
    f = plt.plot(t_lags, t, marker='None', color='k', linestyle='-')

    if high:
       plt.axis([0.9, 100, 0.0001, 1.0])
    else:
       plt.axis([0.9, 100, 0.001, 1.0])
    
    if save:
        name = "/cita/h/home-2/calinliv/Desktop/figures/check/1D_corr_map_%d_ra_fixed.png" % (mode_number)
        #f = open(name, "w")
        plt.savefig(name)

    ######## end 1D thing. ########################################




