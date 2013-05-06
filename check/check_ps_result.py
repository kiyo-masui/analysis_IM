#! /usr/bin/env python 

import numpy as np
from mkpower import ps_summary
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid


rf_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_rf_40mode_2dpow'
tr_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_tr_40mode_2dpow'
ns_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_ns_40mode_2dpow'
ps_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_ps_40mode_2dpow'
si_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_si_40mode_2dpow'

def plot_1d_power_spectrum(rf_root, tr_root, ns_root, ps_root, si_root, filename, 
                           truncate_range = None):
    power_1d, power_1d_err, k_1d_centre =\
        ps_summary.convert_2dps_to_1dps(ps_root, ns_root, tr_root, 
                                        rf_root, truncate_range)
    #print power_1d
    #print power_1d_err

    k_1d_centre_positive = k_1d_centre[power_1d>0]
    power_1d_positive = power_1d[power_1d>0]
    power_1d_positive_err = power_1d_err[power_1d>0]
    power_1d_positive_err_2 = np.zeros(shape=(2, power_1d_positive_err.shape[0]))
    power_1d_positive_err_2[0] = power_1d_positive_err
    power_1d_positive_err_2[1] = power_1d_positive_err
    power_1d_positive_err_2[0].put(
        np.where(power_1d_positive_err_2[0] >= power_1d_positive), 
        power_1d_positive[power_1d_positive_err_2[0] >= power_1d_positive] - 1.e-15)

    k_1d_centre_negative = k_1d_centre[power_1d<0]
    power_1d_negative = -power_1d[power_1d<0]
    power_1d_negative_err = power_1d_err[power_1d<0]
    power_1d_negative_err_2 = np.zeros(shape=(2, power_1d_negative_err.shape[0]))
    power_1d_negative_err_2[0] = power_1d_negative_err
    power_1d_negative_err_2[1] = power_1d_negative_err
    power_1d_negative_err_2[0].put(
        np.where(power_1d_negative_err_2[0] >= power_1d_negative), 
        power_1d_negative[power_1d_negative_err_2[0] >= power_1d_negative] - 1.e-15)

    print 'positive power'
    for i in range(len(k_1d_centre_positive)):
        print k_1d_centre_positive[i], power_1d_positive[i], power_1d_positive_err[i]
    print 'negative power'
    for i in range(len(k_1d_centre_negative)):
        print k_1d_centre_negative[i], power_1d_negative[i], power_1d_negative_err[i]

    power_1d_sim, power_1d_err_sim, k_1d_centre_sim =\
        ps_summary.convert_2dps_to_1dps_sim(si_root, ns_root, 
        truncate_range = truncate_range)

    #print power_1d_sim
    #print power_1d_err_sim

    power_1d_err_2_sim = np.zeros(shape=(2, power_1d_err_sim.shape[0]))
    power_1d_err_2_sim[0] = power_1d_err_sim
    power_1d_err_2_sim[1] = power_1d_err_sim
    power_1d_err_2_sim[0].put(np.where(power_1d_err_2_sim[0] >= power_1d_sim), 
                          power_1d_sim[power_1d_err_2_sim[0] >= power_1d_sim] - 1.e-15)

    fig = plt.figure(figsize=(8, 8))
    label = filename.replace('_', ' ')
    if truncate_range != None:
        title = 'truncated k_p[%5.3f %5.3f] k_v[%5.3f %5.3f]'%tuple(truncate_range)
    else:
        title = 'no truncating'
    plt.errorbar(k_1d_centre_positive, power_1d_positive, power_1d_positive_err_2, 
        fmt='ro', mec='r', capsize=4.5, elinewidth=1, label=label + ' positive')
    plt.errorbar(k_1d_centre_negative, power_1d_negative, power_1d_negative_err_2, 
        fmt='ro', mec='r', mfc='none', capsize=4.5, 
        elinewidth=1, label=label + ' negative')
    plt.errorbar(k_1d_centre_sim, power_1d_sim, power_1d_err_2_sim, 
        fmt='ko', mec='k', mfc='none', capsize=4.5, elinewidth=1, label='simulation')
    plt.loglog()
    plt.ylim(ymin=1.e-12, ymax=1.e-1)
    plt.xlim(xmin=0.025, xmax=1.5)
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('$\Delta^2$ [$K^2$]')
    plt.legend(frameon=False, title=title)
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    if truncate_range != None:
        plt.savefig('./png/' + filename + '_truncated.png', format='png')
    else:
        plt.savefig('./png/' + filename + '.png', format='png')
    plt.show()

def plot_2d_power_spectrum(rf_root, tr_root, ps_root, filename, truncate_range=None):
    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
        ps_summary.load_power_spectrum(ps_root)

    power_spectrum_err, k_p_edges, k_v_edges =\
        ps_summary.load_power_spectrum_err(ps_root)

    transfer_function, k_p_edges, k_v_edges =\
        ps_summary.load_transfer_function(rf_root, tr_root)

    power_spectrum *= transfer_function
    power_spectrum_err *= transfer_function

    power_spectrum = np.ma.array(power_spectrum)
    power_spectrum[power_spectrum==0] = np.ma.masked
    power_spectrum = np.log10(power_spectrum)
    plot_list.append(power_spectrum)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('%s\n 2d power spectrum'%filename)
    
    power_spectrum_err = np.ma.array(power_spectrum_err)
    power_spectrum_err[power_spectrum_err==0] = np.ma.masked
    power_spectrum_err = np.log10(power_spectrum_err)
    plot_list.append(power_spectrum_err)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('%s\n 2d power spectrum error'%filename)
    
    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        power_spectrum, kpe, kve = ps_summary.truncate_2dps(k_p_range, 
                                k_v_range, k_p_edges, k_v_edges, power_spectrum)
        plot_list.append(power_spectrum)
        x_list.append(kve)
        y_list.append(kpe)
        title_list.append('%s\n 2d power spectrum truncated'%filename)

        power_spectrum_err, kpe, kve = ps_summary.truncate_2dps(k_p_range, 
                                k_v_range, k_p_edges, k_v_edges, power_spectrum_err)
        plot_list.append(power_spectrum_err)
        x_list.append(kve)
        y_list.append(kpe)
        title_list.append('%s\n 2d power spectrum error truncated'%filename)

    image_box_2d(x_list, y_list, plot_list, title_list=title_list, n_row=1,
                 xlim=[min(k_v_edges), max(k_v_edges)],
                 ylim=[min(k_p_edges), max(k_p_edges)],
                 clim=[-7, -1])
    plt.savefig('./png/%s_2d_power_spectrum.png'%filename, format='png')


def plot_2d_weight(rf_root, tr_root, ns_root, filename, truncate_range=None):

    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    transfer_function, k_p_edges, k_v_edges =\
        ps_summary.load_transfer_function(rf_root, tr_root)
    
    weight, k_p_edges, k_v_edges =\
        ps_summary.load_weight(ns_root, transfer_function)
    weight = np.ma.array(weight)
    weight[weight==0] = np.ma.masked
    weight = np.log10(weight)
    plot_list.append(weight)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('%s\n noise weight'%filename)
    
    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        weight_truncated, kpe, kve = ps_summary.truncate_2dps(k_p_range, 
                                k_v_range, k_p_edges, k_v_edges, weight)
        plot_list.append(weight_truncated)
        x_list.append(kve)
        y_list.append(kpe)
        title_list.append('%s\n noise weight truncated'%filename)

    image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                 xlim=[min(k_v_edges), max(k_v_edges)],
                 ylim=[min(k_p_edges), max(k_p_edges)],
                 clim=[3, 15])
    plt.savefig('./png/%s_noise_weight.png'%filename, format='png')

def plot_2d_transfer_function(rf_root, tr_root, filename, truncate_range=None):

    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    transfer_function, k_p_edges, k_v_edges =\
        ps_summary.load_transfer_function(rf_root, tr_root)
    transfer_function = np.ma.array(transfer_function)
    transfer_function[transfer_function==0] = np.ma.masked
    transfer_function = np.ma.log10(transfer_function)
    plot_list.append(transfer_function)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)
    title_list.append('%s\n transfer function'%filename)

    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        transfer_function_truncated, kpe, kve = ps_summary.truncate_2dps(k_p_range, 
                                k_v_range, k_p_edges, k_v_edges, transfer_function)
        plot_list.append(transfer_function_truncated)
        x_list.append(kve)
        y_list.append(kpe)
        title_list.append('%s\n transfer function truncated'%filename)

    image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                 xlim=[min(k_v_edges), max(k_v_edges)],
                 ylim=[min(k_p_edges), max(k_p_edges)],
                 clim=[3, 9])
    plt.savefig('./png/%s_transfer_function.png'%filename, format='png')

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

    image_box_2d(x_list, y_list, plot_list, n_row=2, title_list=title_list)
    plt.savefig('./png/everything.png', format='png')
    
def image_box_2d(x_list, y_list, plot_list, n_row=1, n_col=None, title_list=None, 
                 show=False, xlim=None, ylim=None, clim=None, ):
    if n_col == None:
        n_col = int(np.ceil(float(len(plot_list))/float(n_row)))
    
    fig = plt.figure(figsize=(5*n_col+2,5*n_row+3))
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
        if clim != None:
            im.set_clim(clim[0], clim[1])
        ax[i].set_xlabel('k vertical [h/Mpc]')
        ax[i].set_ylabel('k parallel [h/Mpc]')
        if xlim == None:
            xlim = [x.min(), x.max()]
        if ylim == None:
            ylim = [y.min(), y.max()]
        ax[i].set_xlim(xlim[0], xlim[1])
        ax[i].set_ylim(ylim[0], ylim[1])
        if title_list != None:
            ax[i].set_title(title_list[i])
        ax[i].loglog()
        ax[i].cax.colorbar(im)
    
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    if show:
        plt.show()

if __name__=='__main__':

    mode = 25 
    #mode = 10 
    #filename = '15hr_II_14conv_auto_ps_15hour_%dmode'%mode
    filename = '15hr_IE_14conv_auto_ps_15hour_%dmode'%mode
    #filename = '1hr_IE_14conv_auto_ps_01hour_%dmode'%mode

    result_root = '/Users/ycli/DATA/ps_result/ps_result_subreal/'\
                  + filename + '/'
    #result_root = '/Users/ycli/DATA/ps_result/' + filename + '/'
    rf_root = result_root + 'auto_rf_%dmode_2dpow'%mode
    tr_root = result_root + 'auto_tr_%dmode_2dpow'%mode
    ns_root = result_root + 'auto_ns_%dmode_2dpow'%mode
    ps_root = result_root + 'auto_ps_%dmode_2dpow'%mode
    si_root = result_root + 'auto_si_%dmode_2dpow'%mode
    
    truncate_range = [0.035,1.,0.08,0.3]
    #truncate_range = [0.035,1.,0.04,0.3]

    #plot_2d_everthing(rf_root, tr_root, ns_root, ps_root)
    #plot_2d_transfer_function(rf_root, tr_root, filename)
    #plot_2d_transfer_function(rf_root, tr_root, filename, truncate_range)
    #plot_2d_weight(rf_root, tr_root, ns_root, filename, truncate_range)
    #plot_2d_power_spectrum(rf_root, tr_root, ps_root, filename, truncate_range)

    plot_1d_power_spectrum(rf_root, tr_root, ns_root, ps_root, si_root, filename)
    plot_1d_power_spectrum(rf_root, tr_root, ns_root, ps_root, si_root, filename, 
                           truncate_range=truncate_range)
