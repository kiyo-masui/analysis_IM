#! /usr/bin/env python 

import numpy as np
from mkpower import ps_summary
from mkpower import ps_fitting
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import os
import copy


rf_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_rf_40mode_2dpow'
tr_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_tr_40mode_2dpow'
ns_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_ns_40mode_2dpow'
ps_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_ps_40mode_2dpow'
si_root = '/Users/ycli/DATA/ps_result/new_ps_estimation_test/auto_si_40mode_2dpow'

def seperate_positive_and_negative_power(power, power_err, k_centre, lower=1.e-15):
    '''
        This function used to separate the positive and negative power,
        and change the lower bound of the error bar to a tiny positive value, 
        if it was negative
    '''

    k_1d_centre_positive = k_centre[power>0]
    power_1d_positive = power[power>0]
    power_1d_positive_err = power_err[power>0]
    power_1d_positive_err_2 = np.zeros(shape=(2, power_1d_positive_err.shape[0]))
    power_1d_positive_err_2[0] = power_1d_positive_err
    power_1d_positive_err_2[1] = power_1d_positive_err
    power_1d_positive_err_2[0].put(
        np.where(power_1d_positive_err_2[0] >= power_1d_positive), 
        power_1d_positive[power_1d_positive_err_2[0] >= power_1d_positive] - lower)

    k_1d_centre_negative = k_centre[power<0]
    power_1d_negative = -power[power<0]
    power_1d_negative_err = power_err[power<0]
    power_1d_negative_err_2 = np.zeros(shape=(2, power_1d_negative_err.shape[0]))
    power_1d_negative_err_2[0] = power_1d_negative_err
    power_1d_negative_err_2[1] = power_1d_negative_err
    power_1d_negative_err_2[0].put(
        np.where(power_1d_negative_err_2[0] >= power_1d_negative), 
        power_1d_negative[power_1d_negative_err_2[0] >= power_1d_negative] - lower)

    return power_1d_positive, power_1d_positive_err_2, k_1d_centre_positive,\
           power_1d_negative, power_1d_negative_err_2, k_1d_centre_negative

def plot_1d_power_spectrum_opt_multi(file_root, file_name, ps_type, 
                                     save_name='power_opt',
                                     positive_only=False, plot_error=False, 
                                     from_sec=True, sec=['A', 'B', 'C', 'D'],
                                     average=True, truncate_range=None):
    '''
        This function used to plot a batch of optical power spectrum
    '''
    
    color = "bgrcmykw"
    #color = "rcmykw"
    color_index = 0

    if isinstance(file_name, np.ndarray):
        file_name = file_name.tolist()
    if not isinstance(file_name, list):
        file_name = [file_name]

    fig = plt.figure(figsize=(8, 8))

    #k_paper_cole, power_paper_cole, power_err_paper_cole = \
    #    ps_summary.load_2df_ps_from_paper('Cole')
    #k_paper, power_paper, power_err_paper = ps_summary.load_2df_ps_from_paper()
    #plt.errorbar(k_paper, power_paper*1.26, power_err_paper, 
    #    fmt='bo', mec='b', capsize=2.5, elinewidth=1, markersize=5,
    #    label='Percival et. al 2001 x 1.26')
    #plt.errorbar(k_paper_cole, power_paper_cole*1.26, power_err_paper_cole, 
    #    fmt='go', mec='g', capsize=2.5, elinewidth=1, markersize=5,
    #    label='Cole et. al 2005 x 1.26')

    # load eric's result
    #
    #power_eric = np.loadtxt('/Users/ycli/DATA/ps_result_eswitzer/gbt_pwrspec_plots_wsvd/GBT_15hr_avg_fdgp_new_x_GBT_15hr_avg_fdgp_new_onesided_ixi_conv1p4_clipnoise_bfwindow_svdweighted_noiseweight/pk_20modes.dat')
    #power_eric = np.loadtxt('/Users/ycli/DATA/ps_result_eswitzer/gbt_pwrspec_plots/GBT_15hr_avg_fdgp_new_x_GBT_15hr_avg_fdgp_new_onesided_ixi_conv1p4_clipnoise_bfwindow_svdweighted_noiseweight_notrans/pk_20modes.dat')
    #power_eric = np.loadtxt('/Users/ycli/DATA/ps_result_eswitzer/gbt_pwrspec_plots/GBT_15hr_avg_fdgp_new_x_GBT_15hr_avg_fdgp_new_onesided_ixi_conv1p4_clipnoise_bfwindow_svdweighted_self/pk_20modes.dat')

    #power_eric_positive, power_eric_positive_err_2, k_eric_centre_positive,\
    #power_eric_negative, power_eric_negative_err_2, k_eric_centre_negative\
    #= seperate_positive_and_negative_power( power_eric[:,1], 
    #                                        power_eric[:,2], 
    #                                        power_eric[:,0])
    #plt.errorbar(k_eric_centre_positive, power_eric_positive, 
    #        power_eric_positive_err_2, fmt='ks', mec='k', capsize=4.5, 
    #        elinewidth=1, label='Eric positive 20 modes notrans')
    #plt.errorbar(k_eric_centre_negative, power_eric_negative, 
    #        power_eric_negative_err_2, fmt='ks', mec='k', mfc='none', 
    #        capsize=4.5, elinewidth=1,)
    #        #label=label + ' negative')

    #power_eric = np.loadtxt('/Users/ycli/DATA/ps_result_eswitzer/gbt_pwrspec_plots/GBT_15hr_avg_fdgp_new_x_GBT_15hr_avg_fdgp_new_onesided_ixi_conv1p4_clipnoise_bfwindow_svdweighted_noiseweight_beamtrans/pk_20modes.dat')
    #power_eric_positive, power_eric_positive_err_2, k_eric_centre_positive,\
    #power_eric_negative, power_eric_negative_err_2, k_eric_centre_negative\
    #= seperate_positive_and_negative_power( power_eric[:,1], 
    #                                        power_eric[:,2], 
    #                                        power_eric[:,0])
    #plt.errorbar(k_eric_centre_positive, power_eric_positive, 
    #        power_eric_positive_err_2, fmt='gs', mec='g', capsize=4.5, 
    #        elinewidth=1, label='Eric positive 20 modes notrans')
    #plt.errorbar(k_eric_centre_negative, power_eric_negative, 
    #        power_eric_negative_err_2, fmt='gs', mec='g', mfc='none', 
    #        capsize=4.5, elinewidth=1,)
    #        #label=label + ' negative')


    power_th_k = np.logspace(np.log10(0.01), np.log10(5.), 100)
    power_th = ps_summary.load_theory_ps(power_th_k, redshift=0.08, 
                                         cross=ps_type=='cros', in_Jy=False)
    #plt.plot(power_th_k, power_th, 'k', linewidth=2)

    offsets = np.arange(len(file_name)) - 0.5*float(len(file_name) - 1.)
    #offsets = np.arange(len(file_name) + 1) - 0.8*float(len(file_name))

    if average:
        ps_power_list = []
        ps_error_list = []
        ps_kbins_list = []

    for name in file_name:
        if ps_type=='cros':
            result_root = file_root + name
            ps_root = result_root%(ps_type, 'ps')
            sn_root = result_root%(ps_type, 'sn')
            ps_power, ps_error, ps_kbins = \
                plot_1d_power_spectrum_opt(ps_root, sn_root, name.split('/')[0], 
                                           fig=fig, color = color[color_index], 
                                           positive_only = positive_only,
                                           plot_error = plot_error, 
                                           offset = offsets[color_index],
                                           convert_to_K = True)
            if average:
                ps_power_list.append(ps_power)
                ps_error_list.append(ps_error)
                ps_kbins_list.append(ps_kbins)
            #bias = ps_fitting.chisquare_fitting(ps_power, ps_error, ps_kbins, 
            #                            redshift=0.08, cross=ps_type=='cros')
            #plt.plot(power_th_k, bias[0]*power_th, '%s'%color[color_index], 
            #         linewidth=2, linestyle='--',
            #         label=r'$\Omega_{HI}b_{HI}b_{opt} = %3.2f\times10^{-3}$'\
            #                 %bias[0])

        elif ps_type=='auto':
            result_root = file_root + name
            rf_root = result_root%(ps_type, 'rf')
            tr_root = result_root%(ps_type, 'tr')
            if from_sec:
                ps_root = result_root%(ps_type, 'ps') + '_%s'
                ne_root = result_root%(ps_type, 'ne') + '_%s'
                ns_root = result_root%(ps_type, 'ns') + '_%s'
                print ns_root
            else:
                ps_root = result_root%(ps_type, 'ps')
                ns_root = result_root%(ps_type, 'ns')
                ne_root = None

            ps_power, ps_error, ps_kbins = \
                plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, 
                                           name.split('/')[0], 
                                           ne_root=ne_root, sec=sec, fig=fig, 
                                           filename=None, color=color[color_index], 
                                           plot_error=plot_error,
                                           truncate_range = truncate_range,
                                           offset = offsets[color_index], 
                                           from_sec = from_sec,
                                           convert_to_K = False)
            if average:
                ps_power_list.append(ps_power)
                ps_error_list.append(ps_error)
                ps_kbins_list.append(ps_kbins)
            #bias = ps_fitting.chisquare_fitting(ps_power, ps_error, ps_kbins, 
            #                            redshift=0.08, cross=ps_type=='cros')
            #bias = ps_fitting.chisquare_estimation([ps_power], [ps_error], 
            #        [ps_kbins], redshift=0.8, cross=ps_type=='cros')
            #plt.plot(power_th_k, bias[0]*power_th, '%s'%color[color_index], 
            #         linewidth=2, label=r'$\Omega_{HI}b_{HI} = ' + \
            #         r'%3.2f\pm%3.2f\times10^{-3}$'%(bias[0],bias[1]))

        color_index += 1

    if average:
        bias = ps_fitting.chisquare_fitting(ps_power_list, ps_error_list, 
                                            ps_kbins_list, bias=[0.2,0.6,1000],
                                            redshift=0.8, cross=ps_type=='cros',
                                            plot=False)
        print bias
        bias = ps_fitting.chisquare_estimation(ps_power_list, ps_error_list, 
                                               ps_kbins_list, 
                                               redshift=0.8, cross=ps_type=='cros')
        print bias
        if ps_type == 'auto':
            plt.plot(power_th_k, bias[0]**2*power_th, 'k-', linewidth=2,
                     label=r'$\Omega_{HI}b_{HI} = ' + \
                     r'%3.2f\pm%3.2f\times10^{-3}$'%(bias[0],bias[1]))
                     #r'%3.2f^{+%3.2f}_{-%3.2f}\times10^{-3}$'\
                     #%(bias[0], bias[1], bias[2]))
        else:
            plt.plot(power_th_k, bias[0]*power_th, 'k-', linewidth=2,
                     label=r'$\Omega_{HI}b_{HI} = ' + \
                     r'%3.2f\pm%3.2f\times10^{-3}$'%(bias[0],bias[1]))
                     #r'%3.2f^{+%3.2f}_{-%3.2f}\times10^{-3}$'\
                     #%(bias[0], bias[1], bias[2]))
        #plt.plot(power_th_k, bias[1]*power_th, 'k.', linewidth=1)
        #plt.plot(power_th_k, bias[2]*power_th, 'k.', linewidth=1)

        ps_power_list = np.ma.array(ps_power_list)
        ps_error_list = np.ma.array(ps_error_list)
        ps_kbins_list = np.ma.array(ps_kbins_list)

        #mask = ps_power_list < 0
        #ps_power_list[mask] = np.ma.masked
        #ps_error_list[mask] = np.ma.masked
        #ps_kbins_list[mask] = np.ma.masked

        ps_power_average = np.ma.mean(ps_power_list, axis=0)
        ps_error_average = np.ma.sqrt(np.ma.mean(ps_error_list**2, axis=0))
        ps_kbins_average = np.ma.mean(ps_kbins_list, axis=0)

        ps_power_positive, ps_power_positive_err_2, ps_kbins_positive,\
        ps_power_negative, ps_power_negative_err_2, ps_kbins_negative\
            = seperate_positive_and_negative_power(ps_power_average, 
                                                   ps_error_average, 
                                                   ps_kbins_average)

        plt.errorbar(ps_kbins_positive, ps_power_positive, ps_power_positive_err_2, 
            fmt='ks', mec='k', capsize=4.5, elinewidth=1, label='average positive')
        plt.errorbar(ps_kbins_negative, ps_power_negative, ps_power_negative_err_2, 
            fmt='ks', mec='k', mfc='none', capsize=4.5, elinewidth=1,)
    else:
        bias = 0.35*100./75. * 0.63
        #bias = 0.43
        if ps_type == 'auto':
            plt.plot(power_th_k, bias**2*power_th, 'k-', linewidth=2,
                     label=r'$\Omega_{HI}b_{HI} = ' + \
                     r'%3.2f\times10^{-3}$'%bias)
        else:
            plt.plot(power_th_k, bias*power_th, 'k-', linewidth=2,
                     label=r'$\Omega_{HI}b_{HI} = ' + \
                     r'%3.2f\times10^{-3}$'%bias)


    plt.loglog()
    if ps_type=='cros':
        plt.ylim(ymin=1.e-9, ymax=4.e-1)
    elif ps_type=='auto':
        #plt.ylim(ymin=1.e-10, ymax=5.e-5)
        plt.ylim(ymin=1.e-9, ymax=4.e-3) # paper
    #plt.xlim(xmin=0.05, xmax=2.5)
    plt.xlim(xmin=0.09, xmax=0.7) # paper
    plt.xlabel('$k [h/Mpc]$')
    #plt.ylabel('$\Delta^2$ [$K^2$]')
    plt.ylabel('$\Delta^2$ [$P(k)k^3/2\pi^2$]')
    #plt.ylabel('$P(k)$ [$h^{-3}Mpc^3$]')
    plt.legend(frameon=False, loc=2, ncol=1)
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.savefig('./png/' + save_name + '.png', format='png')
    #plt.show()

def plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, label, 
                               ne_root=None, sec = ['A', 'B', 'C', 'D'],
                               si_root=None, fig=None, filename=None, color='r',
                               plot_error=False, truncate_range = None, from_sec=True,
                               offset=0., convert_to_K=False):
    if from_sec:
        power_1d, power_1d_err, k_1d_centre =\
            ps_summary.convert_2dps_to_1dps_sec(ps_root, ns_root, tr_root, rf_root, 
                                                sec, truncate_range)
    else:
        power_1d, power_1d_err, k_1d_centre =\
            ps_summary.convert_2dps_to_1dps(ps_root, ns_root, tr_root, rf_root,
                                            ne_root, sec, truncate_range)

    print power_1d
    print k_1d_centre

    # exclude some k bin for test
    #power_1d[-10] = 0.
    #power_1d_err[-10] = 0.
    #power_1d[-11] = 0.
    #power_1d_err[-11] = 0.

    if convert_to_K:
        power_1d /= 1.1**2
        power_1d_err /= 1.1**2

    power_1d_positive, power_1d_positive_err_2, k_1d_centre_positive,\
    power_1d_negative, power_1d_negative_err_2, k_1d_centre_negative\
        = seperate_positive_and_negative_power(power_1d, power_1d_err, k_1d_centre)

    dk = k_1d_centre_positive[1]/k_1d_centre_positive[0]
    offset = dk**(0.05*offset)
    k_1d_centre_positive *= offset
    k_1d_centre_negative *= offset


    #print 'positive power'
    #for i in range(len(k_1d_centre_positive)):
    #    print k_1d_centre_positive[i],power_1d_positive[i],power_1d_positive_err_2[0][i]
    #print 'negative power'
    #for i in range(len(k_1d_centre_negative)):
    #    print k_1d_centre_negative[i],power_1d_negative[i],power_1d_negative_err_2[0][i]

    #power_1d_sim, power_1d_err_sim, k_1d_centre_sim =\
    #    ps_summary.convert_2dps_to_1dps_sim(si_root, ns_root, 
    #    truncate_range = truncate_range)

    #print power_1d_sim
    #print power_1d_err_sim

    #power_1d_err_2_sim = np.zeros(shape=(2, power_1d_err_sim.shape[0]))
    #power_1d_err_2_sim[0] = power_1d_err_sim
    #power_1d_err_2_sim[1] = power_1d_err_sim
    #power_1d_err_2_sim[0].put(np.where(power_1d_err_2_sim[0] >= power_1d_sim), 
    #                      power_1d_sim[power_1d_err_2_sim[0] >= power_1d_sim] - 1.e-15)

    #power_th = ps_summary.load_theory_ps(k_1d_centre_sim)
    #power_th /= 0.72**3
    #power_th, power_th_k = ps_summary.load_theory_ps(k_1d_centre_sim)

    if fig == None:
        fig = plt.figure(figsize=(8, 8))
        filename = label
    else:
        filename = None

    fmt=color + 'o'
    mec=color

    label = label.replace('_', ' ')
    if truncate_range != None:
        title = 'truncated k_p[%5.3f %5.3f] k_v[%5.3f %5.3f]'%tuple(truncate_range)
    else:
        title = 'no truncating'
    plt.errorbar(k_1d_centre_positive, power_1d_positive, power_1d_positive_err_2, 
        fmt=fmt, mec=mec, capsize=4.5, elinewidth=1, label='%s positive'%label)
    plt.errorbar(k_1d_centre_negative, power_1d_negative, power_1d_negative_err_2, 
        fmt=fmt, mec=mec, mfc='none', capsize=4.5, elinewidth=1,)
        #label=label + ' negative')

    if plot_error:
        plt.step(k_1d_centre, power_1d_err, where='mid', c=color, linestyle='--',
                 label='error of %s'%label)
    #plt.errorbar(k_1d_centre_sim, power_1d_sim, power_1d_err_2_sim, 
    #    fmt='ko', mec='k', mfc='none', capsize=4.5, elinewidth=1, label='simulation')
    #plt.plot(k_1d_centre_sim, power_th, c='r', label='theory')
    #plt.plot(power_th_k, power_th, c='r', label='theory')

    if filename != None:
        plt.loglog()
        plt.ylim(ymin=1.e-12, ymax=1.e-1)
        plt.xlim(xmin=0.025, xmax=1.5)
        plt.xlabel('$k [h/Mpc]$')
        plt.ylabel('$\Delta^2$ [$K^2$]')
        plt.legend(frameon=False, title=title)
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        if truncate_range != None:
            plt.savefig('./png/' + filename + '_truncated.png', format='png')
        else:
            plt.savefig('./png/' + filename + '.png', format='png')
        #plt.show()
    return power_1d, power_1d_err, k_1d_centre

def plot_1d_power_spectrum_opt(ps_root, sn_root, label, si_root=None, from_1d=False,
                               truncate_range = None, fig=None, color='r', 
                               positive_only=False, plotshortnoise=True, 
                               plot_error=False, offset=0.,
                               convert_to_K=False):
    '''
        This function used to plot optical power spectrum
    '''
    if from_1d:
        power_1d, power_1d_err, power_1d_kmn, k_1d_centre =\
            ps_summary.load_power_spectrum_opt_1d(ps_root, sn_root)
    else:
        power_1d, power_1d_err, k_1d_centre, shortnoise =\
            ps_summary.convert_2dps_to_1dps_opt(ps_root, sn_root, truncate_range)

    if convert_to_K:
        power_1d /= 1.1
        power_1d_err /= 1.1

    power_1d_positive, power_1d_positive_err_2, k_1d_centre_positive,\
    power_1d_negative, power_1d_negative_err_2, k_1d_centre_negative\
        = seperate_positive_and_negative_power(power_1d, power_1d_err, k_1d_centre)

    dk = k_1d_centre_positive[1]/k_1d_centre_positive[0]
    offset = dk**(0.04*offset)
    k_1d_centre_positive *= offset
    k_1d_centre_negative *= offset


    if fig == None:
        fig = plt.figure(figsize=(8, 8))
        filename = label
    else:
        filename = None
    label = label.replace('_', ' ')

    fmt = color + 'o'
    mec = color
    plt.errorbar(k_1d_centre_positive, power_1d_positive, power_1d_positive_err_2, 
        fmt=fmt, mec=mec, capsize=4.5, elinewidth=1, label='%s positive'%label)
    if not positive_only:
        plt.errorbar(k_1d_centre_negative, power_1d_negative, 
                     power_1d_negative_err_2, 
                     fmt=fmt, mec=mec, mfc='none', capsize=4.5, 
                     elinewidth=1,)
                     #label=label + ' negative')
    if plot_error:
        plt.step(k_1d_centre, power_1d_err, where='mid', c=color, linestyle='--',
                 label='error of %s'%label)

        #plt.step(k_1d_centre, shortnoise, where='mid', c=color, 
        #         label='short noise of ' + label)

    if filename != None:
        plt.loglog()
        plt.ylim(ymin=1.e-3, ymax=1.e2)
        plt.xlim(xmin=0.01, xmax=2.5)
        plt.xlabel('$k [h/Mpc]$')
        #plt.ylabel('$\Delta^2$ [$K^2$]')
        plt.ylabel('$\Delta^2$ [$P(k)k^3/2\pi^2$]')
        #plt.ylabel('$P(k)$ [$h^{-3}Mpc^3$]')
        plt.legend(frameon=False)
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        plt.savefig('./png/' + filename + '.png', format='png')
        plt.show()
    return power_1d, power_1d_err, k_1d_centre

def plot_1d_power_spectrum_raw(ps_root, filename): 
    power_1d, k_1d_centre = ps_summary.load_power_spectrum_1d(ps_root, raw=True)

    fig = plt.figure(figsize=(8, 8))
    label = filename.replace('_', ' ')

    for i in range(power_1d.shape[0]):
        print power_1d[i]
        #print k_1d_centre
        power_positive, power_positive_err, k_positive,\
        power_negative, power_negative_err, k_negative = \
        seperate_positive_and_negative_power(power_1d[i], power_1d[i], k_1d_centre)
        plt.scatter(k_positive, power_positive, c='g', marker='o', alpha=0.4, 
                    edgecolors = 'g')
        plt.scatter(k_negative, power_negative, c='g', marker='o', alpha=0.4, 
                    edgecolors = 'g', facecolors = 'none')

    plt.loglog()
    plt.ylim(ymin=1.e-12, ymax=1.e-1)
    plt.xlim(xmin=0.025, xmax=1.5)
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('$\Delta^2$ [$K^2$]')
    #plt.legend(frameon=False)
    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    plt.savefig('./png/' + filename + '.png', format='png')
    plt.show()

def plot_1d_power_spectrum_sim(si_root, filename=None):

    power_1d_sim, power_1d_err_sim, k_1d_centre_sim =\
        ps_summary.convert_2dps_to_1dps_sim(si_root)

    power_1d_sim_positive, power_1d_sim_positive_err_2, k_1d_centre_sim_positive,\
    power_1d_sim_negative, power_1d_sim_negative_err_2, k_1d_centre_sim_negative\
        = seperate_positive_and_negative_power(power_1d_sim, 
                                               power_1d_err_sim, 
                                               k_1d_centre_sim)

    power_th, power_th_k = ps_summary.load_theory_ps(k_1d_centre_sim)

    if filename != None:
        fig = plt.figure(figsize=(8, 8))
        plt.errorbar(k_1d_centre_sim_positive, power_1d_sim_positive, 
            power_1d_sim_positive_err_2, fmt='ko', mec='k', mfc='none', 
            capsize=4.5, elinewidth=1, label='simulation')
        #plt.plot(k_1d_centre_sim_positive, power_th, c='r', label='theory')
        plt.plot(power_th_k, power_th, c='r', label='theory')
        plt.loglog()
        plt.ylim(ymin=1.e-12, ymax=1.e-1)
        plt.xlim(xmin=0.025, xmax=1.5)
        plt.xlabel('k [h/Mpc]')
        plt.ylabel('$\Delta^2$ [$K^2$]')
        plt.legend(frameon=False,)
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        plt.savefig('./png/test.png', format='png')
        plt.show()
    return power_1d_sim, power_1d_err_sim, k_1d_centre_sim

def plot_2d_power_spectrum(rf_root, tr_root, ps_root, filename, truncate_range=None):
    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
        ps_summary.load_power_spectrum(ps_root)

    power_spectrum_err, k_p_edges, k_v_edges =\
        ps_summary.load_power_spectrum_err(ps_root)

    if os.path.exists(rf_root) and os.path.exists(tr_root):
        transfer_function, k_p_edges, k_v_edges =\
            ps_summary.load_transfer_function(rf_root, tr_root)
    else:
        print 'No transfer function'
        transfer_function = np.ones(power_spectrum.shape)

    power_spectrum *= transfer_function
    power_spectrum_err *= transfer_function

    power_spectrum_positive = np.ma.array(copy.deepcopy(power_spectrum))
    power_spectrum_positive[power_spectrum<=0] = np.ma.masked
    power_spectrum_positive = np.log10(power_spectrum_positive)
    plot_list.append(power_spectrum_positive)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)

    title_list.append('%s\n 2d power spectrum positive'%filename)

    power_spectrum_negative = np.ma.array(-copy.deepcopy(power_spectrum))
    power_spectrum_negative[power_spectrum>=0] = np.ma.masked
    power_spectrum_negative = np.log10(power_spectrum_negative)
    plot_list.append(power_spectrum_negative)
    x_list.append(k_v_edges)
    y_list.append(k_p_edges)

    title_list.append('%s\n 2d power spectrum negative'%filename)
    
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
                 clim=[-8, -3])
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
    plt.show()

def plot_2d_gauss_error(ne_root, filename, sec=None):

    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    if sec == None:
        sec = ['A', 'B', 'C', 'D']
    for s1 in sec:
        for s2 in sec:
            if s1 == s2:
                continue
            tab = "%s%s"%(s1, s2)
            ne_root_sec = ne_root%tab
            noise_error, k_p_edges, k_v_edges = ps_summary.load_noise_err(ne_root_sec)
            noise_error = np.ma.array(noise_error)
            noise_error[noise_error==0] = np.ma.masked
            noise_error = np.ma.log10(noise_error)
            plot_list.append(noise_error)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)
            title_list.append('%s\n gauss noise %s'\
                             %(filename.replace('_', ' '), tab))
    image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                 xlim=[min(k_v_edges), max(k_v_edges)],
                 ylim=[min(k_p_edges), max(k_p_edges)],
                 n_row = len(sec),
                 clim=[-8, -3],
                 )
    plt.savefig('./png/%s_gauss_power.png'%filename, format='png')

def plot_2d_noise_error(ne_root, filename, sec=None):

    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    if sec == None:
        sec = ['A', 'B', 'C', 'D']
    for s1 in sec:
        for s2 in sec:
            if s1 == s2:
                continue
            tab = "%s%s"%(s1, s2)
            ne_root_sec = ne_root%tab
            noise_error, k_p_edges, k_v_edges = ps_summary.load_noise_err(ne_root_sec)
            noise_error = np.ma.array(noise_error)
            noise_error[noise_error==0] = np.ma.masked
            noise_error = np.ma.log10(noise_error)
            plot_list.append(noise_error)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)
            title_list.append('%s\n noise diag power %s'\
                             %(filename.replace('_', ' '), tab))
    image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                 xlim=[min(k_v_edges), max(k_v_edges)],
                 ylim=[min(k_p_edges), max(k_p_edges)],
                 n_row = len(sec),
                 clim=[-1, 2],
                 )
    plt.savefig('./png/%s_noise_diag_power.png'%filename, format='png')

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
    
    si_root = "/Users/ycli/DATA/ps_result/test_ideal_2df_ps/2df_si_2dpow"
    #plot_1d_power_spectrum_sim(si_root)
    #exit()

    # -----------------------------------------------------------------------
    file_root = "/Users/ycli/DATA/ps_result/"
    file_name_list = [#"FULL_2df_ps_selection_nosubmean", 
                      #"FULL_2df_ps_selection_but_separableweight",
                      #"FULL_2df_ps_selection_submean",
                      #"FULL_2df_ps_separable_submean",
                      #"FULL_2df_ps_separable_100mock",
                      #"FULL_2df_ps_selection_100mock",
                      #"FULL_2df_ps_selection_1000mock",
                      #"FULL_2df_ps_selection_10000mock",
                      #"FULL_PHY_2df_ps",
                      "FULL_2df_ps",
                      #"FULL_2df_ps",
                      #"FULL_2df_ps_selection_submean_oneseed",
                      #"29RA_2df_ps_selection_submean",
                      #"29RA_2df_ps_separable_submean",
                      ]
    #plot_1d_power_spectrum_opt_multi(file_root, file_name_list, 
    #                                 ps_type = '2df',
    #                                 save_name ='2df', 
    #                                 positive_only = True)

    #exit()

    # -----------------------------------------------------------------------

    #parkes_field = 'n03n30_10by7_0627pixel_'
    #parkes_field = '07n30_10by7_0627pixel_'
    #parkes_field = '17n30_10by7_0627pixel_'
    #parkes_field = '27n30_10by7_'
    parkes_field = 'p3500n3000_parkes_2010_10_24-28_'
    #parkes_field = '15hr_41-80_pointcorr_ABCD_'
    #parkes_field = '15hr_41-80_avg_fdgp_new_ABCD_'

    file_root = "/Users/ycli/DATA/ps_result/"
    #file_root = "/Users/ycli/DATA/ps_result/PKS_p3500n3000_parkes_2010_10_24-28_old/"
    ps_type = 'cros'
    #ps_type = 'auto'
    file_name_list = [
        #"PKS_cros_ps/%s_%s_2dpow",

        #"PKS_cros_ps_00mode/%s_%s_0mode_2dpow",
        #"PKS_cros_ps_01mode/%s_%s_1mode_2dpow",
        #"PKS_cros_ps_02mode/%s_%s_2mode_2dpow",
        #"PKS_cros_ps_03mode/%s_%s_3mode_2dpow",
        #"PKS_cros_ps_04mode/%s_%s_4mode_2dpow",

        #"PKS_cros_ps_00mode/%s_%s_0mode_2dpow",
        #"PKS_cros_ps_01mode/%s_%s_1mode_2dpow",
        #"PKS_cros_ps_02mode/%s_%s_2mode_2dpow",
        #"PKS_cros_ps_03mode/%s_%s_3mode_2dpow",

        #"PKS_auto_ps_00mode/%s_%s_0mode_2dpow",
        #"PKS_auto_ps_01mode/%s_%s_1mode_2dpow",
        #"PKS_auto_ps_02mode/%s_%s_2mode_2dpow",
        #"PKS_auto_ps_03mode/%s_%s_3mode_2dpow",

        #"PKS_1pt1_cov_auto_ps_00mode/%s_%s_0mode_2dpow",
        #"PKS_1pt1_cov_auto_ps_01mode/%s_%s_1mode_2dpow",
        #"PKS_1pt1_cov_auto_ps_02mode/%s_%s_2mode_2dpow",
        #"PKS_1pt1_cov_auto_ps_03mode/%s_%s_3mode_2dpow",

        #"PKS_%s1pt1_cov_%s_ps_00mode/"%(parkes_field, ps_type) + "%s_%s_0mode_2dpow",
        #"PKS_%s1pt1_cov_%s_ps_01mode/"%(parkes_field, ps_type) + "%s_%s_1mode_2dpow",
        #"PKS_%s1pt1_cov_%s_ps_02mode/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_%s_ps_03mode/"%(parkes_field, ps_type) + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_%s_ps_04mode/"%(parkes_field, ps_type) + "%s_%s_4mode_2dpow",
        #"PKS_%s1pt1_cov_%s_ps_05mode/"%(parkes_field, ps_type) + "%s_%s_5mode_2dpow",
        #"PKS_%s1pt1_cov_%s_ps_08mode/"%(parkes_field, ps_type) + "%s_%s_8mode_2dpow",


        #"PKS_%s1pt1_cov_bgcal_prior2_%s_ps_02mode/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_prior2En1_%s_ps_02mode/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",

        #"PKS_%s1pt1_cov_bgcal_%s_ps_00mode/"%(parkes_field, ps_type) + "%s_%s_0mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_01mode/"%(parkes_field, ps_type) + "%s_%s_1mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_02mode/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_small_%s_ps_02mode/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_03mode/"%(parkes_field, ps_type) + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_02mode_2dfmaglim19pt35/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_03mode/"%(parkes_field, ps_type) + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%(parkes_field, ps_type) + "%s_%s_5mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_j2k_%s_ps_03mode/"%(parkes_field, ps_type) + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%(parkes_field, ps_type) + "%s_%s_5mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_03mode/"%(parkes_field, ps_type) + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_04mode/"%(parkes_field, ps_type) + "%s_%s_4mode_2dpow",
        #"PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%(parkes_field, ps_type) + "%s_%s_5mode_2dpow",

        #"PKS_%s1pt1_cov_cal_%s_ps_02mode/"%(parkes_field, ps_type) + "%s_%s_2mode_2dpow",

        #"PKS_%s1pt1_cov_fewfreqcut_cros_ps_00mode/"%parkes_field + "%s_%s_0mode_2dpow",
        #"PKS_%s1pt1_cov_fewfreqcut_cros_ps_01mode/"%parkes_field + "%s_%s_1mode_2dpow",
        #"PKS_%s1pt1_cov_fewfreqcut_cros_ps_02mode/"%parkes_field + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_fewfreqcut_cros_ps_03mode/"%parkes_field + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_fewfreqcut_cros_ps_04mode/"%parkes_field + "%s_%s_4mode_2dpow",
        #"GBT_%s1pt4_cov_%s_ps_20mode/"%(parkes_field, ps_type) + "%s_%s_20mode_2dpow",
        #"GBT_%s1pt4_cov_wsvd_%s_ps_15mode/"%(parkes_field, ps_type) + "%s_%s_15mode_2dpow",
        #"GBT_15hr_41-80_pointcorr_ABCD_1pt4_cov_wsvd_%s_ps_15mode/"%(ps_type) + "%s_%s_15mode_2dpow",
        #"GBT_15hr_41-80_pointcorr_ABCD_1pt4_cov_%s_ps_20mode/"%(ps_type) + "%s_%s_20mode_2dpow",
        #"GBT_15hr_41-80_pointcorr_ABCD_1pt4_cov_%s_ps_05mode/"%(ps_type) + "%s_%s_5mode_2dpow",
        #"GBT_15hr_41-80_pointcorr_ABCD_1pt4_cov_%s_ps_10mode/"%(ps_type) + "%s_%s_10mode_2dpow",
        #"GBT_15hr_41-80_pointcorr_ABCD_1pt4_cov_%s_ps_20mode/"%(ps_type) + "%s_%s_20mode_2dpow",
        #"GBT_15hr_41-80_avg_fdgp_new_ABCD_1pt4_cov_%s_ps_20mode/"%(ps_type) + "%s_%s_20mode_2dpow",
        #"GBT_%s1pt4_cov_wsvd_%s_ps_25mode/"%(parkes_field, ps_type) + "%s_%s_25mode_2dpow",

        ]

    #truncate_range = [0.035,1.,0.04,0.3]
    #truncate_range = [0.035,1.,0.08,0.3]
    truncate_range = None

    #plot_1d_power_spectrum_opt_multi(file_root, file_name_list, 
    #        ps_type = ps_type,
    #        save_name ='pks_%s%s_1pt1_cov'%(parkes_field, 
    #                                        ps_type), 
    #        #save_name ='gbt_%s%s_1pt1_cov_sec_compensated'%(
    #        #    parkes_field, ps_type), 
    #        #save_name ='gbt_%s_1pt4_cov_sec_compensated'%(
    #        #    ps_type), 
    #        positive_only = False,
    #        plot_error = True, 
    #        from_sec = True,
    #        average=False,
    #        truncate_range = truncate_range,
    #        sec = ['A', 'B', 'C', 'D'])
    #exit()


    file_name_list = [
            "PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%('07n30_10by7_0627pixel_', ps_type) + "%s_%s_5mode_2dpow",
            "shifttest_PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%('07n30_10by7_0627pixel_', ps_type) + "%s_%s_5mode_2dpow",
            #"PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%('17n30_10by7_0627pixel_', ps_type) + "%s_%s_5mode_2dpow",
            #"PKS_%s1pt1_cov_cal_%s_ps_05mode/"%('27n30_10by7_', ps_type) + "%s_%s_5mode_2dpow",
            #"PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_5mode_2dpow",
            #"nonetest_PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_5mode_2dpow",
            #"PKS_%s1pt1_cov_bgcal_%s_ps_03mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_3mode_2dpow",
            #"GBT_test/%s_%s_20mode_2dpow",
            #"PKS_%s1pt1_cov_%s_ps_05mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_5mode_2dpow",

            #"PKS_%s1pt1_cov_bgcal_%s_ps_03mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_3mode_2dpow",
            #"PKS_%s1pt1_cov_bgcal_%s_ps_01mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_1mode_2dpow",
            #"PKS_%s1pt1_cov_bgcal_%s_ps_05mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_5mode_2dpow",

        #"PKS_%s1pt1_cov_cal_%s_ps_03mode/"%('27n30_10by7_', ps_type) + "%s_%s_3mode_2dpow",
        #"PKS_%s1pt1_cov_cal_j2000test_%s_ps_03mode/"%('27n30_10by7_', ps_type) + "%s_%s_3mode_2dpow",

        #"PKS_%s1pt1_cov_cal_%s_ps_02mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_2mode_2dpow",
        #"PKS_%s1pt1_cov_cal_j2000test_%s_ps_02mode/"%('p3500n3000_parkes_2010_10_24-28_', ps_type) + "%s_%s_2mode_2dpow",
        ]

    plot_1d_power_spectrum_opt_multi(file_root, file_name_list, 
                                     ps_type = ps_type,
                                     save_name ='pks_average_%s_1pt1_cov'%(ps_type), 
                                     positive_only = False,
                                     plot_error = False, 
                                     from_sec = False,
                                     average=False,)
                                     #sec = ['A', 'B', 'C', 'D', 'E'])
    exit()

    # -----------------------------------------------------------------------

    mode = 5
    #mode = 10 
    #filename = '15hr_II_14conv_auto_ps_15hour_%dmode'%mode
    #filename = '15hr_IE_14conv_auto_ps_15hour_%dmode'%mode
    #filename = '1hr_IE_14conv_auto_ps_01hour_%dmode'%mode
    #filename = '15hr_ABCD_14conv_auto_ps_%dmode'%mode
    #filename = 'PKS_auto_ps_%02dmode'%mode
    #filename = 'PKS_1pt1_cov_auto_ps_%02dmode'%mode
    #filename = '15hr_ABCD_14conv_test_me_auto_ps_%dmode'%mode
    #filename = '15hr_ABCD_14conv_test_eric_auto_ps_%dmode'%mode

    #result_root = '/Users/ycli/DATA/ps_result/'\
    #              + filename + '/'

    filename = "PKS_%s1pt1_cov_bgcal_%s_ps_%02dmode"%(parkes_field, ps_type, mode)
    result_root = '/Users/ycli/DATA/ps_result/' + filename + '/'

    #result_root = '/Users/ycli/DATA/ps_result/' + filename + '/'
    rf_root = result_root + '%s_rf_%dmode_2dpow'%(ps_type, mode)
    tr_root = result_root + '%s_tr_%dmode_2dpow'%(ps_type, mode)
    ns_root = result_root + '%s_ns_%dmode_2dpow'%(ps_type, mode)
    ps_root = result_root + '%s_ps_%dmode_2dpow'%(ps_type, mode)
    si_root = result_root + '%s_si_%dmode_2dpow'%(ps_type, mode)
    
    #truncate_range = [0.035,1.,0.08,0.3]
    #truncate_range = [0.035,1.,0.04,0.3]

    filename = "PKS_%s1pt1_cov_bgcal_%s_ps_%02dmode"%('p3500n3000_', ps_type, mode)
    #plot_2d_everthing(rf_root, tr_root, ns_root, ps_root)
    #plot_2d_transfer_function(rf_root, tr_root, filename)
    #plot_2d_transfer_function(rf_root, tr_root, filename, truncate_range)
    #plot_2d_weight(rf_root, tr_root, ns_root, filename, truncate_range)
    #plot_2d_power_spectrum(rf_root, tr_root, ps_root, filename, truncate_range)
    plot_2d_power_spectrum(rf_root, tr_root, ps_root, filename)

    #plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, rf_root, filename)
    #plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, tr_root, filename)
    #plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, si_root, filename)
    #plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, si_root, 
    #                           filename, truncate_range=truncate_range)

    #si_root = result_root + 'sim_raw_x_beammeansub/auto_si_%dmode_2dpow'%mode
    #si_root = result_root + 'sim_beammeansub_x_beammeansub/auto_si_%dmode_2dpow'%mode
    #plot_1d_power_spectrum_gbt(rf_root, tr_root, ns_root, ps_root, si_root, filename)

    #si_root = result_root + 'auto_tr_%dmode_1draw'%mode
    #si_root = result_root + 'auto_si_%dmode_1draw'%mode
    #si_root = result_root + 'sim_raw_x_beammeansub/auto_si_%dmode_1draw'%mode
    #si_root = result_root + 'sim_beammeansub_x_beammeansub/auto_si_%dmode_1draw'%mode
    #plot_1d_power_spectrum_raw(si_root, filename + '_1draw')


    # -----------------------------------------------------------------------

    mode = 0
    filename = "PKS_%s1pt1_cov_bgcal_%s_ps_%02dmode"%(parkes_field, ps_type, mode)
    result_root = '/Users/ycli/DATA/ps_result/' + filename + '/'

    ne_root = result_root + 'auto_ne_%dmode_'%mode + '2dpow_%s'
    ns_root = result_root + 'auto_ns_%dmode_'%mode + '2dpow_%s'

    filename = "PKS_%s1pt1_cov_bgcal_%s_ps_%02dmode"%('p3500n3000_', ps_type, mode)
    #plot_2d_noise_error(ne_root, filename, sec=['A', 'B', 'C', 'D', 'E'])

    #3plot_2d_gauss_error(ns_root, filename, sec=['A', 'B', 'C', 'D', 'E'])






