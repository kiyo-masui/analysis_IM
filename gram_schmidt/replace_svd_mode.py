#! /usr/bin/env python

import gs
import scipy as sp
import numpy as np
import cPickle
import matplotlib.pyplot as plt
import data_paths
import gs
import file_tools as ft

path_key = 'GBT_15hr_map_oldcal_cleaned_noconv'
def get_svd_modes(path_key, pairitem):
    datapath_db = data_paths.DataPath()
    root = datapath_db.fetch(path_key, silent=True)

    filename_svd = root[1]['%s;SVD'%pairitem]
    svd_info = ft.load_pickle(filename_svd)

    return svd_info
    


def process_mode_files(path_key):
    datapath_db = data_paths.DataPath()
    root = datapath_db.fetch(path_key, silent=True)
    #print root[1]['A_with_B;SVD']

    pairs = ["A_with_B", "A_with_C", "A_with_D",
             "B_with_C", "B_with_D", "C_with_D"]

    # open one file to get the dimensions
    modepkl = cPickle.load(open(root[1]['%s;SVD'%pairs[0]], "r"))
    (amp_ind, mode0_ind, mode1_ind) = (0,1,2)
    num_modes = len(modepkl[amp_ind])
    num_freq = len(modepkl[mode0_ind][0])
    num_pairs = len(pairs)

    amplitudes = np.zeros((num_modes, num_pairs))
    mode_functions_l = np.zeros((num_modes, num_freq, num_pairs))
    mode_functions_r = np.zeros((num_modes, num_freq, num_pairs))

    for pair, pairind in zip(pairs, range(num_pairs)):
        modepkl = cPickle.load(open(root[1]['%s;SVD'%pair], "r"))
        amplitudes[:,pairind] = modepkl[amp_ind]
        for modeind in range(num_modes):
            mode_functions_l[modeind, :, pairind] = modepkl[mode0_ind][modeind]
            mode_functions_r[modeind, :, pairind] = modepkl[mode1_ind][modeind]


    return (amplitudes, mode_functions_l, mode_functions_r)

def get_mode_functions(path_key):
    (amplitudes, mode_functions_l, mode_functions_r) = process_mode_files(path_key)
    mode_functions_l_avg = np.mean(mode_functions_l, axis=2)
    return mode_functions_l_avg

def plot_mode_functions(path_key, n=6):
    functions = get_mode_functions(path_key)
    (n_mode, n_freq) = functions.shape
    x = range(n_freq)
    if n>n_mode: n=n_mode
    plt.figure(figsize=(8,7))
    for i in range(n):
        plt.plot(x, functions[i])
    plt.savefig('./png/%s_%dmode.png'%(path_key, n))
    plt.show()
    
# ---------------------------

def replace_modes(mode_functions, n, m=None):
    mode_functions = np.array(mode_functions)
    #print mode_functions.shape
    (n_mode, n_freq) = mode_functions.shape
    mode_functions_new = mode_functions[:n, :]
    #print mode_functions_new.shape
    #print n_freq
    if m==None: m = n_mode
    for i in range(n, m):
        mode_functions_new = gs.gs_array(mode_functions_new)
        #mode_functions_new[-1] *= 0.1
    #print mode_functions_new.shape
    return mode_functions_new

if __name__=='__main__':
    
    svd_info = get_svd_modes(path_key, 'A_with_B')
    svd_info_new = (svd_info[0], 
                    replace_modes(svd_info[1], 4, m=50),
                    replace_modes(svd_info[2], 4, m=50))
    freq_n = len(svd_info[1][0])
    x = range(freq_n)

    plt.figure(figsize=(8,8))
    plt.subplot(611)
    plt.plot(x, svd_info[1][0])
    plt.plot(x, svd_info_new[1][0])
    plt.subplot(612)
    plt.plot(x, svd_info[1][1])
    plt.plot(x, svd_info_new[1][1])
    plt.subplot(613)
    plt.plot(x, svd_info[1][2])
    plt.plot(x, svd_info_new[1][2])
    plt.subplot(614)
    plt.plot(x, svd_info[1][3])
    plt.plot(x, svd_info_new[1][3])
    plt.subplot(615)
    plt.plot(x, svd_info[1][4])
    plt.plot(x, svd_info_new[1][4])
    plt.subplot(616)
    plt.plot(x, svd_info[1][5])
    plt.plot(x, svd_info_new[1][5])
    plt.savefig('./png/%s.png'%(path_key))
    plt.show()

