#! /usr/bin/env python 

import numpy as np
import ps_summary as pss
import os
import copy
from kiyopy import parse_ini

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
    return fig, ax

def plot_2d_power_spectrum(power_spectrum_list, k_mode_number_list, k_p_edges, 
        k_v_edges, filename, label_list=None, output='./'):

    plot_list = []
    x_list = []
    y_list = []
    title_list = []

    if not isinstance(power_spectrum_list, list):
        power_spectrum_list = [power_spectrum_list,]
        k_mode_number_list  = [k_mode_number_list, ]
        label_list = [label_list,]

    for i in range(len(power_spectrum_list)):
        power_spectrum = power_spectrum_list[i]
        k_mode_number  = k_mode_number_list[i]
        label = label_list[i]

        power_spectrum_positive = np.ma.array(copy.deepcopy(power_spectrum))
        power_spectrum_positive[power_spectrum<=0] = np.ma.masked
        power_spectrum_positive = np.log10(power_spectrum_positive)
        plot_list.append(power_spectrum_positive)
        x_list.append(k_v_edges)
        y_list.append(k_p_edges)

        title_list.append('%s\n%s positive'%(filename, label))

        if np.any(power_spectrum<0):
            power_spectrum_negative = np.ma.array(-copy.deepcopy(power_spectrum))
            power_spectrum_negative[power_spectrum>=0] = np.ma.masked
            power_spectrum_negative = np.log10(power_spectrum_negative)
            plot_list.append(power_spectrum_negative)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)

            title_list.append('%s\n%s negative'%(filename, label))

        if k_mode_number != None:
            k_mode_number = np.ma.array(copy.deepcopy(k_mode_number))
            k_mode_number[k_mode_number==0] = np.ma.masked
            k_mode_number = np.log10(k_mode_number)
            plot_list.append(k_mode_number)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)
            title_list.append('%s\n%s k mode'%(filename, label))
    
    fig, ax = image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                           n_row=int(round(float(len(power_spectrum_list))/2.)),
                           xlim=[min(k_v_edges), max(k_v_edges)],
                           ylim=[min(k_p_edges), max(k_p_edges)],
                           clim=[-8, -3])
    plt.savefig('%s/%s_2dpow.png'%(output, filename), format='png')
    fig.clf()
    plt.close()




autops_init = {
        "inputroot" : "./",
        "outputroot" : "./",
        "sections" : ["A", "B", "C", "D"],
        "truncate" : [],
        "mode" : 20,
        }
autopsprefix = 'autops_'

class GBTAutoPowerSpectrum_Analysis(object):
    r""" Get the power spectrum
    """

    def __init__(self, parameter_file=None, param_dict=None, feedback=0):

        self.params = params_dict

        if parameter_file:
            self.params = parse_init.parse(parameter_file, 
                                           autops_init, 
                                           prefix=autopsprefix)

    def execute(self, processes):

        pre = ''
        inputroot = self.params['inputroot']
        mode = self.params['mode']
        rf_root = inputroot + 'auto_rf_%dmode_2dpow'%mode
        tr_root = inputroot + 'auto_tr_%dmode_2dpow'%mode

        ps_root = inputroot + 'auto_ps_%dmode_2dpow'%mode + '_%s'
        ne_root = inputroot + 'auto_ne_%dmode_2dpow'%mode + '_%s'
        ns_root = inputroot + 'auto_ns_%dmode_2dpow'%mode + '_%s'

        # load the power spectrum of each pair
        power_spectrum_list = []
        k_mode_number_list = []
        label_list = []
        k_p_edges = None
        k_v_edges = None
        for i in range(len(self.params['sections'])):
            for j in range(i+1, len(self.params['sections'])):
                sec1 = self.params['sections'][i]
                sec2 = self.params['sections'][j]

                ps_root_sec = ps_root%(sec1 + sec2 + 'x' + sec2 + sec1)
                power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
                    pss.load_power_spectrum(ps_root_sec)
                power_spectrum_list.append(power_spectrum)
                k_mode_number_list.append(k_mode_number_list)
                label_list.append(sec1 + sec2 + 'x' + sec2 + sec1)

        plot_2d_power_spectrum(power_spectrum_list, k_mode_number_list, 
                               k_p_edges, k_v_edges, filename='auto_ps_%dmode'%mode, 
                               label_list=libel_list, output=self.params['output_root'])

        # load the transfer function
        if os.path.exists(rf_root) and os.path.exists(tr_root):
            transfer_function = pss.load_transfer_function(rf_root, tr_root)[0]

            plot_2d_power_spectrum(transfer_function, None, k_p_edges, k_v_edges,
                                   filename = 'auto_tr_%dmode'%mode,
                                   label_list='', output=self.params['output_root'])
        else:
            print "Note: data for transfer function estimation not exists."
            transfer_function = None

        # load the noise weight
        weight = pss.load_weight_sec(inputroot + 'auto_ns_%dmode_2dpow'%mode, None)
        plot_2d_power_spectrum(weight, None, k_p_edges, k_v_edges,
                               filename = 'auto_wt_%dmode'%mode,
                               label_list='', output=self.params['output_root'])

        # compensate the power
        if transfer_function != None:
            pre += 'comp_'
            power_spectrum_list = [ ps * transfer_function for ps in power_spectrum_list]
            weight /= transfer_function ** 2

            plot_2d_power_spectrum(power_spectrum_list, k_mode_number_list, k_p_edges, 
                                   k_v_edges, filename='comp_auto_ps_%dmode'%mode, 
                                   label_list=libel_list, 
                                   output=self.params['output_root'])
            plot_2d_power_spectrum(weight, None, k_p_edges, k_v_edges,
                                   filename = 'comp_auto_wt_%dmode'%mode,
                                   label_list='', output=self.params['output_root'])

        # get the 1d power
        if self.params['truncate'] != []:
            pre += 'trun_'
            truncate_range = self.params['truncate']
        else:
            truncate_range = None
        ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_list, k_mode_number_list, weight, truncate_range)

        ps_1d_positive, ps_1d_positive_err, k_1d_centre_positive,\
        ps_1d_negative, ps_1d_negative_err, k_1d_centre_negative\
            = seperate_positive_and_negative_power(ps_1d_mean, ps_1d_std, k_centre)

        filename = '%sauto_ps_%dmode_1dpow_positive.txt'%(pre, mode)
        ps_1d = np.array([k_1d_centre_positive, ps_1d_positive, ps_1d_positive_err]).T
        ps_1d.savetxt('%s/%s'%(output, filename), ps_1d, fmt='%10.8f', delimiter=' ')
        filename = '%sauto_ps_%dmode_1dpow_negative.txt'%(pre, mode)
        ps_1d = np.array([k_1d_centre_negative, ps_1d_negative, ps_1d_negative_err]).T
        ps_1d.savetxt('%s/%s'%(output, filename), ps_1d, fmt='%10.8f', delimiter=' ')
