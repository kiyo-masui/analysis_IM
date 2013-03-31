#! /usr/bin/env python 

import numpy as np
from mkpower import ps_summary
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid


rf_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_rf_20mode_2dpow'
tr_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_tr_20mode_2dpow'
ns_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_ns_20mode_2dpow'
ps_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_ps_20mode_2dpow'

def plot_1d_power_spectrum(rf_root, tr_root, ns_root, ps_root):
    power_1d, power_1d_err, k_1d_centre =\
        ps_summary.convert_2dps_to_1dps(ps_root, ns_root, tr_root, rf_root)
    print power_1d
    print power_1d_err

    power_1d_err_2 = np.zeros(shape=(2, power_1d_err.shape[0]))
    power_1d_err_2[0] = power_1d_err
    power_1d_err_2[1] = power_1d_err
    power_1d_err_2[0].put(np.where(power_1d_err_2[0] >= power_1d), 
                          power_1d[power_1d_err_2[0] >= power_1d] - 1.e-15)

    power_1d_sim, power_1d_err_sim, k_1d_centre_sim =\
        ps_summary.convert_2dps_to_1dps_sim(rf_root)
    print power_1d_sim
    print power_1d_err_sim

    power_1d_err_2_sim = np.zeros(shape=(2, power_1d_err_sim.shape[0]))
    power_1d_err_2_sim[0] = power_1d_err_sim
    power_1d_err_2_sim[1] = power_1d_err_sim
    power_1d_err_2_sim[0].put(np.where(power_1d_err_2_sim[0] >= power_1d_sim), 
                          power_1d_sim[power_1d_err_2_sim[0] >= power_1d_sim] - 1.e-15)

    fig = plt.figure(figsize=(8, 8))
    plt.errorbar(k_1d_centre, power_1d, power_1d_err_2, 
                 fmt='ro', mec='r', capsize=4.5, elinewidth=1)
    plt.errorbar(k_1d_centre_sim, power_1d_sim, power_1d_err_2_sim, 
                 fmt='ko', mec='k', mfc='none', capsize=4.5, elinewidth=1)
    plt.loglog()
    plt.ylim(ymin=1.e-12, ymax=1.e-1)
    plt.xlim(xmin=0.025, xmax=1.5)
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.show()

def plot_2d_everthing(rf_root, tr_root, ns_root, ps_root):

    plot_list = []
    x_list = []
    y_list = []
    title_list = []
    cmin = -9
    cmax = -3
    
    transfer_function, k_p_edges, k_v_edges =\
        ps_summary.load_transfer_function(rf_root, tr_root)
    transfer_function = np.ma.array(transfer_function)
    transfer_function[transfer_function==0] = np.ma.masked
    transfer_function = np.log10(transfer_function)
    plot_list.append(transfer_function)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('transfer function')
    
    weight, k_p_edges, k_v_edges =\
        ps_summary.load_weight(ns_root, transfer_function)
    weight = np.ma.array(weight)
    weight[weight==0] = np.ma.masked
    weight = np.log10(weight)
    plot_list.append(weight)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('noise weight')
    
    power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
        ps_summary.load_power_spectrum(ps_root)
    k_mode_number = np.ma.array(k_mode_number)
    k_mode_number[k_mode_number==0] = np.ma.masked
    k_mode_number = np.ma.log(k_mode_number)
    plot_list.append(k_mode_number)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('2d k mode number')
    
    power_spectrum = np.ma.array(power_spectrum)
    power_spectrum[power_spectrum==0] = np.ma.masked
    power_spectrum = np.log10(power_spectrum)
    plot_list.append(power_spectrum)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('2d power spectrum')
    
    power_spectrum_compensated = power_spectrum + transfer_function
    plot_list.append(power_spectrum_compensated)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('2d power spectrum compensated')
    
    n_row = 2 
    n_col = int(np.ceil(float(len(plot_list))/float(n_row)))
    
    fig = plt.figure(figsize=(5*n_col,5*n_row))
    ax = ImageGrid(fig, 111,
                   nrows_ncols = (n_row, n_col),
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
    
    for i in range(len(plot_list)):
        x = x_list[i]
        y = y_list[i]
        image = plot_list[i]
        im = ax[i].pcolormesh(x, y, image)
        #im0.set_clim(cmin, cmax)
        ax[i].set_xlabel('k vertical [h/Mpc]')
        ax[i].set_ylabel('k parallel [h/Mpc]')
        ax[i].set_xlim(xmin=x.min(), xmax=x.max())
        ax[i].set_ylim(ymin=y.min(), ymax=y.max())
        ax[i].set_title(title_list[i])
        ax[i].loglog()
        ax[i].cax.colorbar(im)
    
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.show()

if __name__=='__main__':

    result_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/'
    rf_root = result_root + 'auto_rf_20mode_2dpow'
    tr_root = result_root + 'auto_tr_20mode_2dpow'
    ns_root = result_root + 'auto_ns_20mode_2dpow'
    ps_root = result_root + 'auto_ps_20mode_2dpow'
    
    #plot_2d_everthing(rf_root, tr_root, ns_root, ps_root)

    plot_1d_power_spectrum(rf_root, tr_root, ns_root, ps_root)
