#! /usr/bin/env python

import gs
import scipy as sp
import numpy as np
import cPickle
import matplotlib.pyplot as plt
import data_paths
import file_tools as ft
from core import algebra

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

def replace_modes(mode_functions, n, m=None, weight=None):
    mode_functions = np.array(mode_functions)
    #print mode_functions.shape
    (n_mode, n_freq) = mode_functions.shape
    mode_functions_new = mode_functions[:n, :]
    #print mode_functions_new.shape
    #print n_freq
    if m==None: m = n_mode
    for i in range(n, m):
        mode_functions_new = gs.gs_array(mode_functions_new, 
            legendre_order = i-n+1, weight=weight)  # note: start from P1
        #mode_functions_new[-1] *= 0.1
    #print mode_functions_new.shape
    return mode_functions_new

def replace_mode_newleg(mode_functions, n, m=None, weight=None):
    mode_functions = np.array(mode_functions)
    (n_mode, n_freq) = mode_functions.shape
    mode_functions_new = mode_functions[:n, :]
    if m==None: m = n_mode
    mode_functions_new = gs.gs_array_newleg(mode_functions_new, m, n_freq, weight=weight)
    return mode_functions_new

if __name__=='__main__':
    freq_list = range(50,60)
    freq_n_all = 256
    if True :
      freq_list = range(freq_n_all)
      freq_list.remove(6)
      freq_list.remove(7)
      freq_list.remove(8)
      freq_list.remove(15)
      freq_list.remove(16)
      freq_list.remove(18)
      freq_list.remove(19)
      freq_list.remove(20)
      freq_list.remove(21)
      freq_list.remove(22)
      freq_list.remove(37)
      freq_list.remove(80)
      freq_list.remove(103)
      freq_list.remove(104)
      freq_list.remove(105)
      freq_list.remove(106)
      freq_list.remove(107)
      freq_list.remove(108)
      freq_list.remove(130)
      freq_list.remove(131)
      freq_list.remove(132)
      freq_list.remove(133)
      freq_list.remove(134)
      freq_list.remove(171)
      freq_list.remove(175)
      freq_list.remove(177)
      freq_list.remove(179)
      freq_list.remove(182)
      freq_list.remove(183)
      freq_list.remove(187)
      freq_list.remove(189)
      freq_list.remove(192)
      freq_list.remove(193)
      freq_list.remove(194)
      freq_list.remove(195)
      freq_list.remove(196)
      freq_list.remove(197)
      freq_list.remove(198)
      freq_list.remove(201)
      freq_list.remove(204)
      freq_list.remove(208)
      freq_list.remove(209)
      freq_list.remove(212)
      freq_list.remove(213)
      freq_list.remove(218)
      freq_list.remove(219)
      freq_list.remove(229)
      freq_list.remove(233)
      freq_list.remove(237)
      freq_list.remove(244)
      freq_list.remove(254)
      freq_list.remove(255)
    freq_list = tuple(freq_list)

    svd_info = get_svd_modes(path_key, 'A_with_B')
    mode_n = len(svd_info[1])
    freq_n = len(svd_info[1][0])
    svd_info_all = np.zeros(shape=(2, mode_n, freq_n_all))
    for i in range(mode_n):
        np.put(svd_info_all[0][i], freq_list, svd_info[1][i])
        np.put(svd_info_all[1][i], freq_list, svd_info[2][i])
    svd_info_new_all = (svd_info[0], 
                        replace_modes(svd_info_all[0], 4, m=50),
                        replace_modes(svd_info_all[1], 4, m=50))
    svd_info_new = (svd_info[0],
                    np.take(svd_info_new_all[1], freq_list, axis=1),
                    np.take(svd_info_new_all[2], freq_list, axis=1),)

    svd_info_new_all_weighted = (svd_info[0], 
               replace_modes(svd_info_all[0], 4, m=50, weight=svd_info_all[0][0]),
               replace_modes(svd_info_all[1], 4, m=50, weight=svd_info_all[1][0]))
    svd_info_new_weighted = (svd_info[0],
                    np.take(svd_info_new_all_weighted[1], freq_list, axis=1),
                    np.take(svd_info_new_all_weighted[2], freq_list, axis=1),)

    svd_info_new_nojump = (svd_info[0],
                        replace_modes(svd_info[1], 4, m=50),
                        replace_modes(svd_info[2], 4, m=50))

    #x = range(freq_n_all)
    x = range(freq_n)

    plt.figure(figsize=(8,9))
    fig_n = 6
    for i in range(fig_n):
        plt.subplot(int('%d%d%d'%(fig_n, 1, i+1)))
        m = i + 5
        #plt.plot(x, svd_info_all[0][m], label='svd mode %d'%m)
        plt.plot(x, svd_info[1][m], label='svd mode %d'%m)
        #plt.plot(x, svd_info_new[1][m], label='leg mode %d (jump)'%m)
        #plt.plot(x, svd_info_new_nojump[1][m], label='leg mode %d (no jump)'%m)
        plt.plot(x, svd_info_new_weighted[1][m], label='leg mode %d (jump weighted)'%m)
        plt.legend(frameon=False)
    #plt.subplot(611)
    #plt.plot(x, svd_info_all[0][0])
    #plt.plot(x, svd_info_new[1][0])
    #plt.subplot(612)
    #plt.plot(x, svd_info_all[0][1])
    #plt.plot(x, svd_info_new[1][1])
    #plt.subplot(613)
    #plt.plot(x, svd_info_all[0][2])
    #plt.plot(x, svd_info_new[1][2])
    #plt.subplot(614)
    #plt.plot(x, svd_info_all[0][3])
    #plt.plot(x, svd_info_new[1][3])
    #plt.subplot(615)
    #plt.plot(x, svd_info_all[0][4])
    #plt.plot(x, svd_info_new[1][4])
    #plt.subplot(616)
    #plt.plot(x, svd_info_all[0][5])
    #plt.plot(x, svd_info_new[1][5])
    plt.savefig('./png/%s.png'%(path_key))
    plt.show()

