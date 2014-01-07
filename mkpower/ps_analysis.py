#! /usr/bin/env python 

import numpy as np
import ps_summary as pss
import os
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from kiyopy import parse_ini
from core import algebra
from mpi4py import MPI


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

    positive = np.concatenate([k_1d_centre_positive[None,:],
                               power_1d_positive[None, :], 
                               power_1d_positive_err_2], axis=0)

    k_1d_centre_negative = k_centre[power<0]
    power_1d_negative = -power[power<0]
    power_1d_negative_err = power_err[power<0]
    power_1d_negative_err_2 = np.zeros(shape=(2, power_1d_negative_err.shape[0]))
    power_1d_negative_err_2[0] = power_1d_negative_err
    power_1d_negative_err_2[1] = power_1d_negative_err
    power_1d_negative_err_2[0].put(
        np.where(power_1d_negative_err_2[0] >= power_1d_negative), 
        power_1d_negative[power_1d_negative_err_2[0] >= power_1d_negative] - lower)

    negative = np.concatenate([k_1d_centre_negative[None,:],
                               power_1d_negative[None, :], 
                               power_1d_negative_err_2], axis=0)

    return positive, negative

    #return power_1d_positive, power_1d_positive_err_2, k_1d_centre_positive,\
    #       power_1d_negative, power_1d_negative_err_2, k_1d_centre_negative

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
        #ax[i].set_xlabel('k vertical [h/Mpc]')
        ax[i].set_xlabel('k perpendicular [h/Mpc]')
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
        label_list = [label_list,]
    if k_mode_number_list != None:
        if not isinstance(k_mode_number_list, list):
            k_mode_number_list  = [k_mode_number_list, ]
    else:
        k_mode_number = k_mode_number_list

    for i in range(len(power_spectrum_list)):
        power_spectrum = power_spectrum_list[i]
        if k_mode_number != None:
            k_mode_number  = k_mode_number_list[i]
        label = label_list[i]

        power_spectrum_positive = np.ma.array(copy.deepcopy(power_spectrum))
        power_spectrum_positive[power_spectrum<=0] = np.ma.masked
        power_spectrum_positive = np.ma.log10(power_spectrum_positive)
        plot_list.append(power_spectrum_positive.T)
        x_list.append(k_v_edges)
        y_list.append(k_p_edges)

        title_list.append('%s\n%s positive'%(filename, label))

        if np.any(power_spectrum<0):
            power_spectrum_negative = np.ma.array(-copy.deepcopy(power_spectrum))
            power_spectrum_negative[power_spectrum>=0] = np.ma.masked
            power_spectrum_negative = np.ma.log10(power_spectrum_negative)
            plot_list.append(power_spectrum_negative.T)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)

            title_list.append('%s\n%s negative'%(filename, label))

        if k_mode_number != None:
            k_mode_number = np.ma.array(copy.deepcopy(k_mode_number))
            k_mode_number[k_mode_number==0] = np.ma.masked
            k_mode_number = np.ma.log10(k_mode_number.T)
            plot_list.append(k_mode_number)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)
            title_list.append('%s\n%s k mode'%(filename, label))

    cmax = np.ma.max(plot_list)#*0.5
    cmin = np.ma.min(plot_list)#*0.5
    
    fig, ax = image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                           n_row=int(round(float(len(power_spectrum_list))/2.)),
                           xlim=[min(k_v_edges), max(k_v_edges)],
                           ylim=[min(k_p_edges), max(k_p_edges)],
                           clim=[cmin, cmax])
    print "plot " + '%s/%s_2dpow.png'%(output, filename)
    plt.savefig('%s/%s_2dpow.png'%(output, filename), format='png')
    fig.clf()
    plt.close()

crosps_init = {
        "inputroot" : "./",
        "outputroot" : "./",
        "sections" : ["A", "B", "C", "D"],
        "truncate" : [],
        "mode" : 20,
        }
crospsprefix = 'crosps_'


class GBTxWiggleZPowerSpectrum_Analysis(object):
    r""" Get the power spectrum
    """

    def __init__(self, parameter_file_or_dict=None, feedback=1):

        self.params = parse_ini.parse(parameter_file_or_dict, 
                                      crosps_init, 
                                      prefix=crospsprefix,
                                      feedback=feedback)
        self.job_list = []

    def mpiexecute(self, processes=1):
        '''
        The MPI here is useless, bucause the caclulation is quick enough. 
        Just for making this module working with MPI.
        '''
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        comm.barrier()

        if rank == 0:
            self.execute(processes=1)

        comm.barrier()

    def execute(self, processes):

        pre = ''
        inputroot = self.params['inputroot']
        outputroot = self.params['outputroot']
        if not os.path.exists(outputroot):
            os.makedirs(outputroot)
        mode = self.params['mode']

        rf_root = inputroot + 'cros_rf_%dmode_2dpow'%mode
        tr_root = inputroot + 'cros_tr_%dmode_2dpow'%mode

        ps_root = inputroot + 'cros_ps_%dmode_2dpow'%mode
        ne_root = inputroot + 'cros_ne_%dmode_2dpow'%mode
        #ns_root = inputroot + 'cros_ns_%dmode_2dpow'%mode

        # load the power spectrum
        k_p_edges = None
        k_v_edges = None
        power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(ps_root)

        plot_2d_power_spectrum(power_spectrum, None, 
                               k_p_edges, k_v_edges, filename='cros_ps_%dmode'%mode, 
                               label_list='', output=outputroot)

        # load the short noise
        short_noise_root = inputroot + 'cros_sn_%dmode_2dpow'%mode
        short_noise, shrot_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(short_noise_root)
        plot_2d_power_spectrum(short_noise, None, 
                               k_p_edges, k_v_edges, filename='shortnoise_%dmode'%mode, 
                               label_list='', output=outputroot)

        power_spectrum -= short_noise

        # load the random power
        power_spectrum_list = []
        k_mode_number_list = []
        power_2d_raw_root = inputroot + 'cros_sn_%dmode_2draw'%mode
        power_2d_raw = np.load(power_2d_raw_root)
        power_2d_raw_kmn_root = inputroot + 'cros_sn_%dmode_2draw_kmn'%mode
        power_2d_raw_kmn = np.load(power_2d_raw_kmn_root)

        # load the transfer function
        if os.path.exists(rf_root) and os.path.exists(tr_root):
            transfer_function = pss.load_transfer_function(rf_root, tr_root, False)[0]

            plot_2d_power_spectrum(transfer_function, None, k_p_edges, k_v_edges,
                                   filename = 'cros_tr_%dmode'%mode,
                                   label_list='', output=outputroot)
        else:
            print "Note: data for transfer function estimation not exists."
            transfer_function = None

        # compensate the power
        if transfer_function != None:
            pre += 'comp_'
            power_spectrum *= transfer_function

            power_2d_raw *= transfer_function[None,...]

            plot_2d_power_spectrum(power_spectrum, None, k_p_edges, k_v_edges, 
                                   filename='comp_auto_ps_%dmode'%mode, 
                                   label_list=label_list, output=outputroot)

        # get the 1d power
        if self.params['truncate'] != []:
            pre += 'trun_'
            truncate_range = self.params['truncate']
        else:
            truncate_range = None
        ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum, k_mode_number, None, truncate_range)

        # get the 1d error
        if self.job_list == []:
            self.job_list = range(power_2d_raw.shape[0])
        for i in self.job_list:
            power_each = algebra.make_vect(power_2d_raw[i],
                    axis_names=power_spectrum.info['axes'])
            power_each.info = power_spectrum.info
            kmode_each = algebra.make_vect(power_2d_raw_kmn[i],
                    axis_names=power_spectrum.info['axes'])
            kmode_each.info = power_spectrum.info
            power_spectrum_list.append(power_each)
            k_mode_number_list.append(kmode_each)

        sn_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_list, k_mode_number_list, None, truncate_range)

        positive, negative = seperate_positive_and_negative_power(ps_1d_mean, 
                                                                  ps_1d_std, 
                                                                  k_centre)

        filename = '%scros_ps_%dmode_1dpow_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%10.8e')

        filename = '%scros_ps_%dmode_1dpow_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%10.8e')



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

    def __init__(self, parameter_file_or_dict=None, feedback=1):

        self.params = parse_ini.parse(parameter_file_or_dict, 
                                      autops_init, 
                                      prefix=autopsprefix,
                                      feedback=feedback)

    def mpiexecute(self, nprocesses=1):
        '''
        The MPI here is useless, bucause the caclulation is quick enough. 
        Just for making this module working with MPI.
        '''
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        comm.barrier()

        if rank == 0:
            self.execute(nprocesses=1)

        comm.barrier()


    def execute(self, processes):

        pre = ''
        inputroot = self.params['inputroot']
        outputroot = self.params['outputroot']
        if not os.path.exists(outputroot):
            os.makedirs(outputroot)
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
                print ps_root_sec
                power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
                    pss.load_power_spectrum(ps_root_sec)
                power_spectrum_list.append(power_spectrum)
                k_mode_number_list.append(k_mode_number)
                label_list.append(sec1 + sec2 + 'x' + sec2 + sec1)

        plot_2d_power_spectrum(power_spectrum_list, None, 
                               k_p_edges, k_v_edges, filename='auto_ps_%dmode'%mode, 
                               label_list=label_list, output=outputroot)

        # load the transfer function
        if os.path.exists(rf_root) and os.path.exists(tr_root):
            transfer_function = pss.load_transfer_function(rf_root, tr_root)[0]

            plot_2d_power_spectrum(transfer_function, None, k_p_edges, k_v_edges,
                                   filename = 'auto_tr_%dmode'%mode,
                                   label_list='', output=outputroot)
        else:
            print "Note: data for transfer function estimation not exists."
            transfer_function = None

        # load the noise weight
        weight = pss.load_weight_sec(inputroot + 'auto_ns_%dmode_2dpow'%mode, None)
        plot_2d_power_spectrum(weight, None, k_p_edges, k_v_edges,
                               filename = 'auto_wt_%dmode'%mode,
                               label_list='', output=outputroot)

        # compensate the power
        if transfer_function != None:
            pre += 'comp_'
            power_spectrum_list = [ ps * transfer_function for ps in power_spectrum_list]
            transfer_function[transfer_function==0] = np.inf
            weight /= transfer_function ** 2
            transfer_function[np.isinf(transfer_function)] = 0.

            plot_2d_power_spectrum(power_spectrum_list, None, k_p_edges, k_v_edges, 
                                   filename='comp_auto_ps_%dmode'%mode, 
                                   label_list=label_list, output=outputroot)
            plot_2d_power_spectrum(weight, None, k_p_edges, k_v_edges,
                                   filename = 'comp_auto_wt_%dmode'%mode,
                                   label_list='', output=outputroot)

        # get the 1d power
        if self.params['truncate'] != []:
            pre += 'trun_'
            truncate_range = self.params['truncate']
        else:
            truncate_range = None
        ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_list, k_mode_number_list, weight, truncate_range)

        positive, negative = seperate_positive_and_negative_power(ps_1d_mean, 
                                                                  ps_1d_std, 
                                                                  k_centre)

        filename = '%sauto_ps_%dmode_1dpow_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%10.8e')

        filename = '%sauto_ps_%dmode_1dpow_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%10.8e')

if __name__=="__main__":
    autops_init = {
        "inputroot" : "/home/ycli/data/ps_result/GBT_15hr_41-80_avg_fdgp_new_ABCD_1pt4_cov_auto_ps_20mode/",
        "outputroot" : "/home/ycli/data/ps_result/GBT_15hr_41-80_avg_fdgp_new_ABCD_1pt4_cov_auto_ps_20mode/",
        "sections" : ["A", "B", "C", "D"],
        "truncate" : [],
        "mode" : 20,
        }
    #GBTAutoPowerSpectrum_Analysis(autops_init).execute(1)

    crosps_init = {
        "inputroot" : "/home/ycli/data/ps_result/GBT_15hr_41-80_avg_fdgp_new_ABCD_1pt4_cov_cros_ps_20mode/",
        "outputroot" : "/home/ycli/data/ps_result/GBT_15hr_41-80_avg_fdgp_new_ABCD_1pt4_cov_cros_ps_20mode/",
        "sections" : ["A", "B", "C", "D"],
        "truncate" : [],
        "mode" : 20,
        }
    GBTxWiggleZPowerSpectrum_Analysis(autops_init).execute(1)
