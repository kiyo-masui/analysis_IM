#! /usr/bin/env python 

import numpy as np
import ps_summary as pss
import os
import copy
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import LogLocator
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from kiyopy import parse_ini
from core import algebra
from mpi4py import MPI
from simulations import corr21cm
import h5py
import math
import pdb

def seperate_positive_and_negative_power(power, power_err, k_centre, lower=1.e-11):
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
                 show=False, xlim=None, ylim=None, clim=None, log_scale=False):
    if n_col == None:
        n_col = int(np.ceil(float(len(plot_list))/float(n_row)))

    fig = plt.figure(figsize=(4*n_col+2,4*n_row+1))
    plt.gcf().subplots_adjust(bottom=0.15)
    ax = ImageGrid(fig, 111,
                   nrows_ncols = (n_row, n_col),
                   direction = "row",
                   axes_pad = 0.5,
                   add_all = True,
                   label_mode = "L",
                   share_all = True,
                   cbar_location = "right",
                   cbar_mode = "each",
                   cbar_size = "4%",
                   cbar_pad = 0.05,
                   )
    
    for i in range(len(plot_list)):
        x = x_list[i]
        y = y_list[i]
        image = plot_list[i]
        if log_scale:
            im = ax[i].pcolormesh(x, y, image, norm=matplotlib.colors.LogNorm(), vmin=6.e-8, vmax=6.e-3)
        else:
            im = ax[i].pcolormesh(x, y, image)

        if clim != None:
            im.set_clim(clim[0], clim[1])
        #ax[i].set_xlabel('k vertical [h/Mpc]')
        #ax[i].set_xlabel('k perpendicular [h/Mpc]')
        #ax[i].set_ylabel('k parallel [h/Mpc]')
        ax[i].set_xlabel('$k_\perp\/[h/\mathrm{Mpc}]$')
        ax[i].set_ylabel('$k_\parallel\/[h/\mathrm{Mpc}]$')
        if xlim == None:
            xlim = [x.min(), x.max()]
        if ylim == None:
            ylim = [y.min(), y.max()]
        ax[i].set_xlim(xlim[0], xlim[1])
        ax[i].set_ylim(ylim[0], ylim[1])
        if title_list != None:
            ax[i].set_title(title_list[i])
        ax[i].loglog()
        if log_scale:
            cbar = fig.colorbar(im, ax=ax[i], cax=ax[i].cax, 
                    ticks = LogLocator(subs=range(10)))
        else:
            cbar = fig.colorbar(im, ax=ax[i], cax=ax[i].cax, 
                    ticks = LinearLocator())
        cbar.ax.tick_params(length=3, width=1.)
        #cbar.ax.tick_params(which='minor', length=3, width=1.)
        #ax[i].cax.colorbar(im)
        #cbar.ax.set_ylabel(title_list[i])
        #cbar.ax.minorticks_on()
        #cbar.ax.set_yticklabels(cbar.ax.yaxis.get_ticklabels(), rotation=90)
        #ax[i].cax.yaxis.set_minor_formatter(FormatStrFormatter("%2.1e"))
        #ax[i].cax.yaxis.set_major_formatter(FormatStrFormatter("%2.1e"))

        #minorticks = im.norm(np.arange(1, 10, 2))
        #cbar.ax.yaxis.set_ticks(minorticks, minor=True)

    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)
    if show:
        plt.show()
    return fig, ax

def plot_2d_dipole(ps_2d_kpp, ps_2d_kpn, k_p_edges, k_v_edges, filename, 
        label=None, output='./'):

    fig = plt.figure(figsize=(8,6))

    ax_p_kpp = fig.add_axes([0.07, 0.50, 0.38, 0.38])
    ax_p_kpn = fig.add_axes([0.07, 0.10, 0.38, 0.38])
    cax_p = fig.add_axes([0.42, 0.30, 0.015, 0.40])
    #cax_p.set_ylabel('Positive Power')

    ax_n_kpp = fig.add_axes([0.52, 0.50, 0.38, 0.38])
    ax_n_kpn = fig.add_axes([0.52, 0.10, 0.38, 0.38])
    cax_n = fig.add_axes([0.87, 0.30, 0.015, 0.40])
    #cax_n.set_ylabel('Negative Power')

    ax_p_kpp.set_aspect('equal')
    ax_n_kpp.set_aspect('equal')
    ax_p_kpn.set_aspect('equal')
    ax_n_kpn.set_aspect('equal')

    ax_p_kpp.loglog()
    ax_p_kpn.loglog()
    ax_n_kpp.loglog()
    ax_n_kpn.loglog()

    ax_p_kpp.set_ylabel('$\mathrm{positive}\, k_\parallel \mathrm{ } [h/\mathrm{Mpc}]$')
    ax_p_kpn.set_ylabel('$\mathrm{negative}\, k_\parallel \mathrm{ } [h/\mathrm{Mpc}]$')
    ax_p_kpn.set_xlabel('$k_\perp \mathrm{ } [h/\mathrm{Mpc}]$')
    ax_p_kpp.set_xticks([])
    ax_n_kpp.set_ylabel('$\mathrm{positive}\, k_\parallel \mathrm{ } [h/\mathrm{Mpc}]$')
    ax_n_kpn.set_ylabel('$\mathrm{negative}\, k_\parallel \mathrm{ } [h/\mathrm{Mpc}]$')
    ax_n_kpn.set_xlabel('$k_\perp \mathrm{ } [h/\mathrm{Mpc}]$')
    ax_n_kpp.set_xticks([])

    ax_p_kpp.set_xlim(xmin=k_v_edges.min(), xmax=k_v_edges.max())
    ax_p_kpn.set_xlim(xmin=k_v_edges.min(), xmax=k_v_edges.max())
    ax_n_kpp.set_xlim(xmin=k_v_edges.min(), xmax=k_v_edges.max())
    ax_n_kpn.set_xlim(xmin=k_v_edges.min(), xmax=k_v_edges.max())
    ax_p_kpp.set_ylim(ymin=k_p_edges.min(), ymax=k_p_edges.max())
    ax_p_kpn.set_ylim(ymin=k_p_edges.min(), ymax=k_p_edges.max())
    ax_n_kpp.set_ylim(ymin=k_p_edges.min(), ymax=k_p_edges.max())
    ax_n_kpn.set_ylim(ymin=k_p_edges.min(), ymax=k_p_edges.max())

    X, Y = np.meshgrid(k_p_edges, k_v_edges)

    #print ps_2d_kpp.max(), ps_2d_kpp.min(), 
    #print ps_2d_kpn.max(), ps_2d_kpn.min(), 

    ps_p_kpp = np.ma.array(copy.deepcopy(ps_2d_kpp))
    ps_p_kpp[ps_2d_kpp<=0] = np.ma.masked
    ps_p_kpn = np.ma.array(copy.deepcopy(ps_2d_kpn))
    ps_p_kpn[ps_2d_kpn<=0] = np.ma.masked
    ps_n_kpp = np.ma.array(-copy.deepcopy(ps_2d_kpp))
    ps_n_kpp[ps_2d_kpp>=0] = np.ma.masked
    ps_n_kpn = np.ma.array(-copy.deepcopy(ps_2d_kpn))
    ps_n_kpn[ps_2d_kpn>=0] = np.ma.masked
    minmin = np.zeros_like(ps_p_kpp) + 6.e-8
    maxmax = np.zeros_like(ps_p_kpp) + 7.e-3
    vmin = np.ma.min([ps_p_kpp, ps_p_kpn, ps_n_kpp, ps_n_kpn])
    vmax = np.ma.max([ps_p_kpp, ps_p_kpn, ps_n_kpp, ps_n_kpn])
    im_p = ax_p_kpp.pcolormesh(X, Y, ps_p_kpp.T, cmap='jet', vmax=vmax, vmin=vmin,
            norm=matplotlib.colors.LogNorm())
    ax_p_kpn.pcolormesh(X, Y, ps_p_kpn.T, cmap='jet', vmax=vmax, vmin=vmin,
            norm=matplotlib.colors.LogNorm())
    im_n = ax_n_kpp.pcolormesh(X, Y, ps_n_kpp.T, cmap='jet', vmax=vmax, vmin=vmin,
            norm=matplotlib.colors.LogNorm())
    ax_n_kpn.pcolormesh(X, Y, ps_n_kpn.T, cmap='jet', vmax=vmax, vmin=vmin,
            norm=matplotlib.colors.LogNorm())

    cbar_p = fig.colorbar(im_p, cax=cax_p)
    cbar_n = fig.colorbar(im_n, cax=cax_n)

    ax_p_kpp.set_title('Positive Power')
    ax_n_kpp.set_title('Negative Power')

    ax_p_kpn.invert_yaxis()
    ax_n_kpn.invert_yaxis()

    #cax_p.minorticks_on()
    #cax_n.minorticks_on()

    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)

    print "plot " + '%s/%s_2dpow_dipole.png'%(output, filename)
    plt.savefig('%s/%s_2dpow_dipole.png'%(output, filename), format='png')
    fig.clf()
    plt.close()

def plot_2d_power_spectrum(power_spectrum_list, k_mode_number_list, k_p_edges, 
        k_v_edges, filename, label_list=None, output='./', logscale=True, 
        cmax = None, cmin = None):

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
        if k_mode_number_list != None:
            k_mode_number  = k_mode_number_list[i]
        label = label_list[i]

        power_spectrum_positive = np.ma.array(copy.deepcopy(power_spectrum))
        power_spectrum_positive[power_spectrum<=0] = np.ma.masked
        #if logscale:
        #    power_spectrum_positive = np.ma.log10(power_spectrum_positive)
        plot_list.append(power_spectrum_positive.T)
        x_list.append(k_v_edges)
        y_list.append(k_p_edges)

        #title_list.append('%s\n%s positive'%(filename, label))
        title_list.append('%s %s positive'%(filename, label))

        if cmax == None: # or np.ma.max(power_spectrum_positive) > cmax:
            cmax = np.ma.max(power_spectrum_positive)
        if cmin == None: # or np.ma.min(power_spectrum_positive) < cmin:
            cmin = np.ma.min(power_spectrum_positive)

        if np.any(power_spectrum<0):
            power_spectrum_negative = np.ma.array(-copy.deepcopy(power_spectrum))
            power_spectrum_negative[power_spectrum>=0] = np.ma.masked
            #if logscale:
            #    power_spectrum_negative = np.ma.log10(power_spectrum_negative)
            plot_list.append(power_spectrum_negative.T)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)

            #title_list.append('%s\n%s negative'%(filename, label))
            title_list.append('%s %s negative'%(filename, label))

            if cmax == None: # or np.ma.max(power_spectrum_negative) > cmax:
                cmax = np.ma.max(power_spectrum_negative)
            if cmin == None: # or np.ma.min(power_spectrum_negative) < cmin:
                cmin = np.ma.min(power_spectrum_negative)

        k_mode_number=None
        if k_mode_number != None:
            k_mode_number = np.ma.array(copy.deepcopy(k_mode_number))
            k_mode_number[k_mode_number==0] = np.ma.masked
            #if logscale:
            #    k_mode_number = np.ma.log10(k_mode_number.T)
            plot_list.append(k_mode_number.T)
            x_list.append(k_v_edges)
            y_list.append(k_p_edges)
            #title_list.append('%s\n%s k mode'%(filename, label))
            title_list.append('%s %s k mode'%(filename, label))

            # if plot k number, do not set color scale
            cmax = None

    #if cmax == None:
    #    cmax = np.ma.max(plot_list)#*0.5
    #if cmin == None:
    #    cmin = np.ma.min(plot_list)#*0.5
    if cmax == None:
        clim = None
    else:
        clim = [cmin, cmax]
    title_list = None
    fig, ax = image_box_2d(x_list, y_list, plot_list, title_list=title_list, 
                           n_row=int(round(float(len(power_spectrum_list))/2.)),
                           xlim=[min(k_v_edges), max(k_v_edges)],
                           ylim=[min(k_p_edges), max(k_p_edges)],
                           clim=clim, log_scale=logscale)
                           #clim=[cmin, cmax])
    print "plot " + '%s/%s_2dpow.png'%(output, filename)
    plt.savefig('%s/%s_2dpow.png'%(output, filename), format='png')
    fig.clf()
    plt.close()

optps_init = {
        "inputroot" : "./",
        "outputroot" : "./",
        "truncate" : [],
        }
optpsprefix = 'optps_'


class OpticalPowerSpectrum_Analysis(object):
    r""" Get the power spectrum
    """

    def __init__(self, parameter_file_or_dict=None, feedback=1):

        self.params = parse_ini.parse(parameter_file_or_dict, 
                                      crosps_init, 
                                      prefix=optpsprefix,
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

        ps_root = inputroot + '2df_ps_2dpow'

        # load the power spectrum
        k_p_edges = None
        k_v_edges = None
        power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(ps_root)

        plot_2d_power_spectrum(power_spectrum, None, 
                               k_p_edges, k_v_edges, filename='2df_ps',
                               label_list='', output=outputroot)

        # load the short noise
        short_noise_root = inputroot + '2df_sn_2dpow'
        short_noise, shrot_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(short_noise_root)
        plot_2d_power_spectrum(short_noise, None, 
                               k_p_edges, k_v_edges, filename='shortnoise',
                               label_list='', output=outputroot)

        power_spectrum -= short_noise

        # load the random power
        power_spectrum_list = []
        k_mode_number_list = []
        power_2d_raw_root = inputroot + '2df_sn_2draw'
        power_2d_raw = np.load(power_2d_raw_root)
        power_2d_raw_kmn_root = inputroot + '2df_sn_2draw_kmn'
        power_2d_raw_kmn = np.load(power_2d_raw_kmn_root)


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

        filename = '%s2df_ps_1dpow_positive.txt'%(pre)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

        filename = '%s2df_ps_1dpow_negative.txt'%(pre)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

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

        print 'Attempting to open ', inputroot + 'ps_result.hd5', ' from directory ', os.getcwd()
        ps_result = h5py.File(inputroot + 'ps_result.hd5', 'r')

        rf_root = 'cros_rf_%dmode_2dpow'%mode
        rb_root = 'cros_rb_%dmode_2dpow'%mode
        tr_root = 'cros_tr_%dmode_2dpow'%mode

        ps_root = 'cros_ps_%dmode_2dpow'%mode
        gal_root = '2df_ps_%dmode_2draw'%mode
        ne_root = 'cros_ne_%dmode_2dpow'%mode
        #ns_root = 'cros_ns_%dmode_2dpow'%mode

        k_p_edges = None
        k_v_edges = None
        # plot the simulation result
        si_root = 'cros_si_%dmode_2draw'%mode
        if rf_root in ps_result and tr_root in ps_result:

            power_spectrum_rf, k_mode_number_rf, k_p_edges, k_v_edges =\
                pss.load_power_spectrum(rf_root, ps_result)
            plot_2d_power_spectrum(power_spectrum_rf, k_mode_number_rf,
                                   k_p_edges, k_v_edges,
                                   filename = 'cros_rf_%dmode'%mode,
                                   label_list='', output=outputroot,)

            power_spectrum_rb, k_mode_number_rb, k_p_edges, k_v_edges =\
                pss.load_power_spectrum(rb_root, ps_result)
            plot_2d_power_spectrum(power_spectrum_rb, k_mode_number_rb,
                                   k_p_edges, k_v_edges,
                                   filename = 'cros_rb_%dmode'%mode,
                                   label_list='', output=outputroot,)

            power_spectrum_tr, k_mode_number_tr, k_p_edges, k_v_edges =\
                pss.load_power_spectrum(tr_root, ps_result)
            plot_2d_power_spectrum(power_spectrum_tr, k_mode_number_tr,
                                   k_p_edges, k_v_edges,
                                   filename = 'cros_tr_%dmode'%mode,
                                   label_list='', output=outputroot,)

            transfer_function_beam = \
                    pss.load_transfer_function(rf_root, rb_root, ps_result, False)[0]
            transfer_function_beam[ transfer_function_beam == 0 ] = np.inf
            plot_2d_power_spectrum(1./transfer_function_beam, None, k_p_edges, k_v_edges,
                                   filename = 'cros_transferf_beam_%dmode'%mode,
                                   label_list='', output=outputroot,
                                   logscale=False, cmax=1., cmin=0.3)

            transfer_function_mode = \
                    pss.load_transfer_function(rb_root, tr_root, ps_result, False)[0]
            transfer_function_mode[ transfer_function_mode == 0 ] = np.inf
            plot_2d_power_spectrum(1./transfer_function_mode, None, k_p_edges, k_v_edges,
                                   filename = 'cros_transferf_mode_%dmode'%mode,
                                   label_list='', output=outputroot,
                                   logscale=False, cmax=1., cmin=0.3)

            transfer_function = \
                    pss.load_transfer_function(rf_root, tr_root, ps_result, False)[0]
            transfer_function[ transfer_function == 0 ] = np.inf
            tr = algebra.make_vect(1./transfer_function,
                    axis_names=power_spectrum_tr.info['axes'])
            #algebra.save(outputroot + 'cros_transferf_%dmode'%mode, tr)
            plot_2d_power_spectrum(1./transfer_function, None, k_p_edges, k_v_edges,
                                   filename = 'cros_transferf_%dmode'%mode,
                                   label_list='', output=outputroot,
                                   logscale=False, cmax=1., cmin=0.3)
            transfer_function[ np.isinf(transfer_function)] = 0.
        else:
            print "Note: data for transfer function estimation not exists."
            transfer_function = None

        # the following is for the "density field slope" analysis. It computes
        # the ratio of the cross power to the 2dF auto power.
        if gal_root in ps_result:

        #    power_spectrum_sim, k_mode_number_sim, k_p_edges, k_v_edges =\
        #        pss.load_power_spectrum('cros_si_%dmode_2dpow'%mode, ps_result)

            power_spectrum_2df, k_mode_number_2df, k_p_edges, k_v_edges =\
                pss.load_power_spectrum('2df_ps_%dmode_2dpow'%mode, ps_result)
        #    power_spectrum_2df[power_spectrum_2df == 0] = np.inf


            plot_2d_power_spectrum(power_spectrum_2df, None,
                                   k_p_edges, k_v_edges,
                                   filename = '2df_ps_%dmode'%mode,
                                   label_list='', output=outputroot,
                                   cmin = 1.e-6, cmax=1.e-3)
            if tr_root in ps_result:
                power_spectrum_sim *= transfer_function
                plot_2d_power_spectrum(power_spectrum_sim / power_spectrum_2df, None,
                                       k_p_edges, k_v_edges,
                                       filename = 'comp_temp_si_%dmode'%mode,
                                       label_list='', output=outputroot,
                                       cmin = 1.e-6, cmax=1.e-3)

            #power_spectrum_2df[power_spectrum_2df == np.inf] = 0
            #plot_2d_power_spectrum(power_spectrum_2df, None,
            #                       k_p_edges, k_v_edges,
            #                       filename = '2df_si_%dmode'%mode,
            #                       label_list='', output=outputroot)

            #power_2d_raw_sim = ps_result[si_root].value
            #power_2d_raw_kmn_sim = ps_result['cros_si_%dmode_2draw_kmn'%mode].value
            #power_2d_raw_sim_kpp = ps_result['cros_si_%dmode_2draw_kpp'%mode].value
            #power_2d_raw_sim_kpn = ps_result['cros_si_%dmode_2draw_kpn'%mode].value
            #power_2d_sim_kmn_kpp = ps_result['cros_si_%dmode_2draw_kmn_kpp'%mode].value
            #power_2d_sim_kmn_kpn = ps_result['cros_si_%dmode_2draw_kmn_kpn'%mode].value
            #power_2d_raw_2df = ps_result[gal_root].value
            #power_2d_raw_kmn_2df = ps_result['2df_si_%dmode_2draw_kmn'%mode].value
            #power_2d_raw_2df_kpp = ps_result['2df_si_%dmode_2draw_kpp'%mode].value
            #power_2d_raw_2df_kpn = ps_result['2df_si_%dmode_2draw_kpn'%mode].value
            #power_2d_2df_kmn_kpp = ps_result['2df_si_%dmode_2draw_kmn_kpp'%mode].value
            #power_2d_2df_kmn_kpn = ps_result['2df_si_%dmode_2draw_kmn_kpn'%mode].value

            #if transfer_function != None:
            #    power_2d_raw_sim     *= transfer_function[None,...]
            #    power_2d_raw_sim_kpp *= transfer_function[None,...]
            #    power_2d_raw_sim_kpn *= transfer_function[None,...]

            sim_list = []
            sim_kpp_list = []
            sim_kpn_list = []
            sim_kmn_list = []
            sim_kmn_kpp_list = []
            sim_kmn_kpn_list = []
            ratio_list = []
            ratio_kmn_list = []
            for i in range(power_2d_raw_sim.shape[0]):
                #prepping for monopole, quadrupole
                sim_each = algebra.make_vect(power_2d_raw_sim[i],
                        axis_names=power_spectrum_sim.info['axes'])
                gal_each = algebra.make_vect(power_2d_raw_2df[i],
                        axis_names=power_spectrum_2df.info['axes'])
                gal_kmn_each = algebra.make_vect(power_2d_raw_kmn_2df[i],
                        axis_names=power_spectrum_2df.info['axes'])
                gal_kmn_each.info = power_spectrum_2df.info
                gal_each[gal_each == 0] = np.inf
                sim_each /= gal_each
                sim_each.info = power_spectrum_sim.info
                kmn_each = algebra.make_vect(power_2d_raw_kmn_sim[i],
                        axis_names=power_spectrum_sim.info['axes'])
                kmn_each.info = power_spectrum_sim.info
                sim_list.append(sim_each)
                sim_kmn_list.append(kmn_each)
                ratio_list.append(sim_each)
                ratio_kmn_list.append(gal_kmn_each)

            std = np.std(np.array(sim_list), axis=0)
            var = std**2
            var[var == 0] = np.inf
            var = 1 / var
            gal_var = var
            power_ratio = np.mean(np.array(sim_list), axis=0)

            filename = '%stemp_si_%dmode_value_std.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), np.vstack(\
                    (power_ratio.flatten()[power_ratio.flatten() \
                    != 0],std.flatten()[power_ratio.flatten() != 0]\
                    )),fmt='%15.12e')
            ratio_1d_mean, ratio_1d_std, ratio_k_centre = pss.convert_2dps_to_1dps_each(
                    ratio_list, ratio_kmn_list, var, None)
            ratio_positive, ratio_negative = seperate_positive_and_negative_power(
                    ratio_1d_mean, ratio_1d_std, ratio_k_centre)
            gal_var = var

            filename = '%stemp_si_%dmode_1dpow_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), ratio_positive.T, fmt='%15.12e')
            filename = '%stemp_si_%dmode_1dpow_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), ratio_negative.T, fmt='%15.12e')

        # load the short noise
        short_noise_root = 'cros_sn_%dmode_2dpow'%mode
        short_noise, shrot_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(short_noise_root, ps_result)
        plot_2d_power_spectrum(short_noise, None, 
                               k_p_edges, k_v_edges, filename='shortnoise_%dmode'%mode, 
                               label_list='', output=outputroot)
        sn_1d_mean, sn_1d_std, sn_k_centre = pss.convert_2dps_to_1dps_each(
                    short_noise, shrot_noise_kmn, None, None)
        sn_positive, sn_negative = seperate_positive_and_negative_power(
                    sn_1d_mean, sn_1d_std, sn_k_centre)

        filename = '%scros_sn_%dmode_1dpow_mean.txt'%(pre,mode)
        np.savetxt('%s/%s'%(outputroot, filename), sn_1d_mean, fmt='%15.12e')
        filename = '%scros_sn_%dmode_1dpow_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), sn_positive.T, fmt='%15.12e')
        filename = '%scros_sn_%dmode_1dpow_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), sn_negative.T, fmt='%15.12e')
        #if os.path.exists(si_root):
        if si_root in ps_result:

            power_spectrum_sim, k_mode_number_sim, k_p_edges, k_v_edges =\
                pss.load_power_spectrum('cros_si_%dmode_2dpow'%mode, ps_result)

            #power_spectrum_sim -= short_noise

            plot_2d_power_spectrum(power_spectrum_sim, k_mode_number_sim,
                                   k_p_edges, k_v_edges,
                                   filename = 'cros_si_%dmode'%mode,
                                   label_list='', output=outputroot,)

            if transfer_function != None:
                power_spectrum_sim *= transfer_function
                #ps_2d_kpp      *= transfer_function
                #ps_2d_kpn	   *= transfer_function

            plot_2d_power_spectrum(power_spectrum_sim, k_mode_number_sim, 
                                   k_p_edges, k_v_edges,
                                   filename = 'comp_cros_si_%dmode'%mode,
                                   label_list='', output=outputroot,)

            power_2d_raw_sim = ps_result[si_root].value
            power_2d_raw_kmn_sim = ps_result['cros_si_%dmode_2draw_kmn'%mode].value
            power_2d_raw_sim_kpp = ps_result['cros_si_%dmode_2draw_kpp'%mode].value
            power_2d_raw_sim_kpn = ps_result['cros_si_%dmode_2draw_kpn'%mode].value
            power_2d_sim_kmn_kpp = ps_result['cros_si_%dmode_2draw_kmn_kpp'%mode].value
            power_2d_sim_kmn_kpn = ps_result['cros_si_%dmode_2draw_kmn_kpn'%mode].value

            #power_2d_raw_sim -= short_noise

            if transfer_function != None:
                print 'Applying transfer to sims.'
                power_2d_raw_sim     *= transfer_function[None,...]
                #power_2d_raw_sim_kpp *= transfer_function[None,...]
                #power_2d_raw_sim_kpn *= transfer_function[None,...]

            sim_list = []
            sim_kpp_list = []
            sim_kpn_list = []
            sim_kmn_list = []
            sim_kmn_kpp_list = []
            sim_kmn_kpn_list = []
            for i in range(power_2d_raw_sim.shape[0]):
                #prepping for monopole, quadrupole
                sim_each = algebra.make_vect(power_2d_raw_sim[i], 
                        axis_names=power_spectrum_sim.info['axes'])
                sim_each.info = power_spectrum_sim.info
                kmn_each = algebra.make_vect(power_2d_raw_kmn_sim[i], 
                        axis_names=power_spectrum_sim.info['axes'])
                kmn_each.info = power_spectrum_sim.info
                sim_list.append(sim_each)
                sim_kmn_list.append(kmn_each)

                #following two blocks prep for dipole
                sim_each = algebra.make_vect(power_2d_raw_sim_kpp[i],
                        axis_names=power_spectrum_sim.info['axes'])
                sim_each.info = power_spectrum_sim.info
                kmn_each = algebra.make_vect(power_2d_sim_kmn_kpp[i],
                        axis_names=power_spectrum_sim.info['axes'])
                kmn_each.info = power_spectrum_sim.info
                sim_kpp_list.append(sim_each)
                sim_kmn_kpp_list.append(kmn_each)

                sim_each = algebra.make_vect(power_2d_raw_sim_kpn[i],
                        axis_names=power_spectrum_sim.info['axes'])
                sim_each.info = power_spectrum_sim.info
                kmn_each = algebra.make_vect(power_2d_sim_kmn_kpn[i],
                        axis_names=power_spectrum_sim.info['axes'])
                kmn_each.info = power_spectrum_sim.info
                sim_kpn_list.append(sim_each)
                sim_kmn_kpn_list.append(kmn_each)
                
            std = np.std(np.array(sim_list), axis=0)
            var = std**2
            var[var == 0] = np.inf
            var = var**(-1)

            #monopole
            si_1d_mean, si_1d_std, si_k_centre = pss.convert_2dps_to_1dps_each(
                    sim_list, sim_kmn_list, var, None)
            print si_k_centre, si_1d_mean, si_1d_std
            si_positive, si_negative = seperate_positive_and_negative_power(
                    si_1d_mean, si_1d_std, si_k_centre)
            si_1d_mean = np.concatenate((si_k_centre[None,:], si_1d_mean[None, :],\
                        si_1d_std[None,:]), axis=0)
            filename = '%scros_si_%dmode_1dpow.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_1d_mean.T, fmt='%15.12e')
            filename = '%scros_si_%dmode_1dk.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_k_centre, fmt='%15.12e')
            filename = '%scros_si_%dmode_1dpow_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_positive.T, fmt='%15.12e')
            filename = '%scros_si_%dmode_1dpow_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_negative.T, fmt='%15.12e')
            
            plot_2d_power_spectrum(std, None,
                                   k_p_edges, k_v_edges, 
                                   filename = 'cros_std_si_%dmode'%mode,
                                   label_list='', output=outputroot,)
            np.save(outputroot + 'weights', std)
            np.save(outputroot + 'k_vals.npy', k_v_edges)
            #single simulation
            #si_single, si_std, si_k_centre = pss.convert_2dps_to_1dps_each(
            #         sim_list[0], sim_kmn_list[0], var, None)
            #si_single, si_std, si_k_centre = pss.convert_2dps_to_1dps_each(
            #         sim_list[0], sim_kmn_list[0], var, None)
            #si_pos_single, si_neg_single = seperate_positive_and_negative_power(
            #        si_single, si_std, si_k_centre)
            #filename = '%scros_si_%dmode_1dpow_single_positive.txt'%(pre, mode)
            #np.savetxt('%s/%s'%(outputroot, filename), si_pos_single.T, fmt='%15.12e')
            #filename = '%scros_si_%dmode_1dpow_single_negative.txt'%(pre, mode)
            #np.savetxt('%s/%s'%(outputroot, filename), si_neg_single.T, fmt='%15.12e')

            #1D covariance
            si_1d_cov, si_k_centre, si_k_edges = pss.covar_1d(
                    sim_list, sim_kmn_list, var, None)
            si_1d_cov = np.nan_to_num(si_1d_cov)
            plot_2d_power_spectrum(si_1d_cov, None,
                                   si_k_edges, si_k_edges,
                                   filename = 'cros_cov_si_%dmode'%mode,
                                   label_list='', output=outputroot,
                                   logscale=False, cmax=1.0, cmin=-1.0)            

            #dipole
            si_kpp_mean, si_kpp_std, si_k_centre = pss.convert_2dps_to_1dps_each(
                    sim_kpp_list, sim_kmn_kpp_list, var, None, order=1)
            si_kpn_mean, si_kpn_std, si_k_centre = pss.convert_2dps_to_1dps_each(
                    sim_kpn_list, sim_kmn_kpn_list, var, None, order=1)
            si_dp_std = np.sqrt((si_kpp_std**2 + si_kpn_std**2))
            si_dp_mean = 0.5 * (si_kpp_mean + si_kpn_mean)
            si_dp_positive, si_dp_negative = seperate_positive_and_negative_power(
                    si_dp_mean, si_dp_std, si_k_centre)
            filename = '%scros_si_%dmode_dipole_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_dp_positive.T, fmt='%15.12e')
            filename = '%scros_si_%dmode_dipole_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_dp_negative.T, fmt='%15.12e')

            #quadrupole
            si_qp_mean, si_qp_std, si_k_centre = pss.convert_2dps_to_1dps_each(
                    sim_list, sim_kmn_list, var, None, order=2)
            si_qp_positive, si_qp_negative = seperate_positive_and_negative_power(
                    si_qp_mean, si_qp_std, si_k_centre)
            filename = '%scros_si_%dmode_quadrupole_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_qp_positive.T, fmt='%15.12e')
            filename = '%scros_si_%dmode_quadrupole_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_qp_negative.T, fmt='%15.12e')
            #print 'si_qp_positive = ', si_qp_positive, si_qp_positive.shape
        
        # load the transfer function
        #if os.path.exists(rf_root) and os.path.exists(tr_root):

        #if os.path.exists(inputroot + 'cros_nl_%dmode_2dpow'%mode):
        noise_weight = None
        if 'cros_nl_%dmode_2dpow'%mode in ps_result:

            noise_root = 'cros_nl_%dmode_2dpow'%mode
            noise_power, k_mode_number, k_p_edges, k_v_edges  = \
                    pss.load_power_spectrum(noise_root, ps_result)
            plot_2d_power_spectrum(noise_power, None, 
                    k_p_edges, k_v_edges, filename='noiselevel_%dmode'%mode, 
                    label_list='', output=outputroot)
            noise_1d_mean, noise_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    noise_power, k_mode_number, None, None)

            filename = '%scros_nl_%dmode_1dpow.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), np.concatenate(
                [k_centre[None,:], noise_1d_mean[None, :]], axis=0).T, fmt='%15.12e')

            if transfer_function != None:
                noise_power *= transfer_function**2
                plot_2d_power_spectrum(noise_power, None, 
                        k_p_edges, k_v_edges, filename='comp_noiselevel_%dmode'%mode, 
                        label_list='', output=outputroot)
                noise_1d_mean, noise_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                        noise_power, k_mode_number, None, None)

                filename = 'comp_%scros_nl_%dmode_1dpow.txt'%(pre, mode)
                np.savetxt('%s/%s'%(outputroot, filename), np.concatenate(
                    [k_centre[None,:], noise_1d_mean[None, :]], axis=0).T, fmt='%15.12e')

            noise_power[noise_power==0] = np.inf
            noise_weight = 1./noise_power

        #if not os.path.exists(ps_root):
        if not ps_root in ps_result:
            return
        # load the short noise
        short_noise_root = 'cros_sn_%dmode_2dpow'%mode
        short_noise, shrot_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(short_noise_root, ps_result)
        plot_2d_power_spectrum(short_noise, None, 
                               k_p_edges, k_v_edges, filename='shortnoise_%dmode'%mode, 
                               label_list='', output=outputroot)

        short_noise_root = '2df_sn_%dmode_2dpow'%mode
        #shot_noise_2df, sn_2df_kmn, k_p_edges, k_v_edges =\
        #    pss.load_power_spectrum(short_noise_root, ps_result)

        # load the power spectrum
        power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(ps_root, ps_result)
        power_spectrum -= short_noise
        plot_2d_power_spectrum(power_spectrum, None, 
                               k_p_edges, k_v_edges, filename='cros_ps_%dmode'%mode, 
                               label_list='', output=outputroot)

        #power_spectrum_2df, k_mode_number_2df, k_p_edges, k_v_edges =\
        #    pss.load_power_spectrum('2df_ps_%dmode_2dpow'%mode, ps_result)
        #power_spectrum_2df -= shot_noise_2df
        #plot_2d_power_spectrum(power_spectrum_2df, None,
        #                       k_p_edges, k_v_edges, filename='2df_ps_%dmode'%mode,
        #                       label_list='', output=outputroot)

        #power_spectrum_2df[power_spectrum_2df == 0] = np.inf
        #plot_2d_power_spectrum(power_spectrum / power_spectrum_2df, None,
        #                       k_p_edges, k_v_edges,
        #                       filename = 'temp_ps_%dmode'%mode,
        #                       label_list='', output=outputroot,
        #                       cmin = 1.e-6, cmax=1.e-3)

        power_spectrum_tr[power_spectrum_tr == 0] = np.inf
        plot_2d_power_spectrum(power_spectrum / power_spectrum_tr, None,
                               k_p_edges, k_v_edges, filename='ps_over_tr',
                               label_list='', output=outputroot)
        power_spectrum_rf[power_spectrum_rf == 0] = np.inf
        #power_spectrum_2df[power_spectrum_2df==np.inf] = 0
        simobj = corr21cm.Corr21cm()
        T_b = simobj.T_b(0.08) * 1e-3 * 0.85
        #plot_2d_power_spectrum(power_spectrum_2df / (power_spectrum_rf / T_b),
        #                       None, k_p_edges, k_v_edges, filename='2df_over_dm',
        #                       label_list='', output=outputroot)
        #bias, bias_std, bias_k_centre = pss.convert_2dps_to_1dps_each(
        #        np.sqrt(power_spectrum_2df / (power_spectrum_rf / T_b)), k_mode_number_2df, \
        #        gal_var / (power_spectrum_rf / T_b), None)
        #bias_pos, bias_neg = seperate_positive_and_negative_power(
        #        bias, bias_std, bias_k_centre)

        #filename='cros_bias_%dmode_1dpow_positive.txt'%(mode)
        #np.savetxt('%s/%s'%(outputroot, filename), bias_pos.T, fmt='%15.12e')
        #filename='cros_bias_%dmode_1dpow_negative.txt'%mode
        #np.savetxt('%s/%s'%(outputroot, filename), bias_neg.T, fmt='%15.12e')

        # plot 2d dipole
        short_noise_root = 'cros_sn_%dmode_2dpow_kpp'%mode
        short_noise_kpp, shrot_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(short_noise_root, ps_result)
        short_noise_root = 'cros_sn_%dmode_2dpow_kpn'%mode
        short_noise_kpn, shrot_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum(short_noise_root, ps_result)
        plot_2d_dipole(short_noise_kpp, short_noise_kpn, k_p_edges, k_v_edges, 
                filename='shortnoise_%dmode'%mode, 
                label='', output=outputroot)

        ps_2d_kpp, k_mode_number, k_p_edges, k_v_edges =\
                pss.load_power_spectrum(ps_root + '_kpp', ps_result)
        ps_2d_kpp -= short_noise_kpp
        ps_2d_kpn, k_mode_number, k_p_edges, k_v_edges =\
                pss.load_power_spectrum(ps_root + '_kpn', ps_result)
        ps_2d_kpn -= short_noise_kpn
        plot_2d_dipole(ps_2d_kpp, ps_2d_kpn, k_p_edges, k_v_edges, 
                filename='cros_ps_%dmode'%mode, 
                label='', output=outputroot)

        # load the random power
        power_2d_raw_root = 'cros_sn_%dmode_2draw'%mode
        power_2d_raw = ps_result[power_2d_raw_root]
        power_2d_raw_kmn_root = 'cros_sn_%dmode_2draw_kmn'%mode
        power_2d_raw_kmn = ps_result[power_2d_raw_kmn_root]
        # for dipole
        power_2d_raw_root = 'cros_sn_%dmode_2draw_kpp'%mode
        power_2d_raw_kpp = ps_result[power_2d_raw_root]
        power_2d_raw_kmn_root = 'cros_sn_%dmode_2draw_kmn_kpp'%mode
        power_2d_raw_kmn_kpp = ps_result[power_2d_raw_kmn_root]
        power_2d_raw_root = 'cros_sn_%dmode_2draw_kpn'%mode
        power_2d_raw_kpn = ps_result[power_2d_raw_root]
        power_2d_raw_kmn_root = 'cros_sn_%dmode_2draw_kmn_kpn'%mode
        power_2d_raw_kmn_kpn = ps_result[power_2d_raw_kmn_root]


        # truncate the 2df power
        if self.params['truncate'] != []:
            pre += 'trun_'
            truncate_range = self.params['truncate']
        else:
            truncate_range = None

        # get the 1d power without compensation
        ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum, k_mode_number, None, truncate_range)
        # get the 1d quadrupole without comensation
        qp_1d_mean, qp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum, k_mode_number, None, truncate_range, order=2)
        # get dipole
        dp_kpp_mean, dp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                ps_2d_kpp, k_mode_number, None, truncate_range, order=1)
        dp_kpn_mean, dp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                -ps_2d_kpn, k_mode_number, None, truncate_range, order=1)
        dp_mean = 0.5 * (dp_kpp_mean + dp_kpn_mean)

        # get the 1d error without compensation
        power_spectrum_list = []
        k_mode_number_list = []
        power_spectrum_kpp_list = []
        k_mode_number_kpp_list = []
        power_spectrum_kpn_list = []
        k_mode_number_kpn_list = []
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

            power_each = algebra.make_vect(power_2d_raw_kpp[i],
                    axis_names=power_spectrum.info['axes'])
            power_each.info = power_spectrum.info
            kmode_each = algebra.make_vect(power_2d_raw_kmn_kpp[i],
                    axis_names=power_spectrum.info['axes'])
            kmode_each.info = power_spectrum.info
            power_spectrum_kpp_list.append(power_each)
            k_mode_number_kpp_list.append(kmode_each)

            power_each = algebra.make_vect(power_2d_raw_kpn[i],
                    axis_names=power_spectrum.info['axes'])
            power_each.info = power_spectrum.info
            kmode_each = algebra.make_vect(power_2d_raw_kmn_kpn[i],
                    axis_names=power_spectrum.info['axes'])
            kmode_each.info = power_spectrum.info
            power_spectrum_kpn_list.append(-power_each)
            k_mode_number_kpn_list.append(kmode_each)

        sn_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_list, k_mode_number_list, None, truncate_range)


        filename = '%scros_ps_%dmode_1dshn.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), np.concatenate( 
                [k_centre[None,:], sn_1d_mean[None, :], 
                ps_1d_std[None, :], ps_1d_std[None, :]], 
                axis=0).T, fmt='%15.12e')

        sn_positive, sn_negative = seperate_positive_and_negative_power(sn_1d_mean, 
                                                                        ps_1d_std, 
                                                                        k_centre)
        filename = '%scros_ps_%dmode_1dshn_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), sn_positive.T, fmt='%15.12e')

        filename = '%scros_ps_%dmode_1dshn_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), sn_negative.T, fmt='%15.12e')

        positive, negative = seperate_positive_and_negative_power(ps_1d_mean, 
                                                                  ps_1d_std, 
                                                                  k_centre)

        filename = '%scros_ps_%dmode_1dpow_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

        filename = '%scros_ps_%dmode_1dpow_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

        # for quadrupole
        sn_1d_mean, qp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_list, k_mode_number_list, 
                None, truncate_range, order=2)
        positive, negative = seperate_positive_and_negative_power(qp_1d_mean, 
                                                                  qp_1d_std, 
                                                                  k_centre)
        filename = '%scros_ps_%dmode_quadrupole_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

        filename = '%scros_ps_%dmode_quadrupole_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')



        # for dipole
        sn_1d_mean, dp_kpp_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_kpp_list, k_mode_number_kpp_list, 
                None, truncate_range, order=1)
        sn_1d_mean, dp_kpn_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_kpn_list, k_mode_number_kpn_list, 
                None, truncate_range, order=1)
        dp_std = np.sqrt(0.5 * (dp_kpp_std**2 + dp_kpn_std**2))
        positive, negative = seperate_positive_and_negative_power(dp_mean, 
                                                                  dp_std, 
                                                                  k_centre)
        filename = '%scros_ps_%dmode_dipole_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

        filename = '%scros_ps_%dmode_dipole_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

        # compensate the power
        if transfer_function != None:
            pre += 'comp_'
            power_spectrum *= transfer_function
            ps_2d_kpp      *= transfer_function
            ps_2d_kpn      *= transfer_function

            power_2d_raw     *= transfer_function[None,...]
            power_2d_raw_kpp *= transfer_function[None,...]
            power_2d_raw_kpn *= transfer_function[None,...]

            plot_2d_power_spectrum(power_spectrum, None, k_p_edges, k_v_edges, 
                                   filename='comp_cros_ps_%dmode'%mode, 
                                   label_list='', output=outputroot)
            #power_spectrum_2df[power_spectrum_2df == 0] = np.inf
            #plot_2d_power_spectrum(power_spectrum / power_spectrum_2df, None,
            #                       k_p_edges, k_v_edges,
            #                       filename = 'comp_temp_ps_%dmode'%mode,
            #                       label_list='', output=outputroot,
            #                       cmin = 1.e-6, cmax=1.e-3)
            #np.save('ratio_array.npy', power_spectrum / power_spectrum_2df)

            # plot S/N
            if noise_weight != None:
                plot_2d_power_spectrum(power_spectrum**2 * noise_weight, None, 
                                       k_p_edges, k_v_edges, 
                                       filename='comp_cros_S2N_%dmode'%mode, 
                                       label_list='', output=outputroot,
                                       logscale=False, cmax=1., cmin=0.)

            plot_2d_dipole(ps_2d_kpp, ps_2d_kpn, k_p_edges, k_v_edges, 
                    filename='comp_cros_ps_%dmode'%mode, 
                    label='', output=outputroot)

            # get the 1d power
            ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum, k_mode_number, var, truncate_range)
            #temp_1d_mean, temp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
            #        power_spectrum / power_spectrum_2df, k_mode_number, gal_var, 
            #        truncate_range)

            # get the 1d quadrupole without comensation
            qp_1d_mean, qp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum, k_mode_number, None, truncate_range, order=2)
            # get dipole
            dp_kpp_mean, dp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    ps_2d_kpp, k_mode_number, None, truncate_range, order=1)
            dp_kpn_mean, dp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    -ps_2d_kpn, k_mode_number, None, truncate_range, order=1)
            dp_mean = 0.5 * (dp_kpp_mean + dp_kpn_mean)


            # get the 1d error
            power_spectrum_list = []
            k_mode_number_list = []
            power_spectrum_kpp_list = []
            k_mode_number_kpp_list = []
            power_spectrum_kpn_list = []
            k_mode_number_kpn_list = []
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

                power_each = algebra.make_vect(power_2d_raw_kpp[i],
                        axis_names=power_spectrum.info['axes'])
                power_each.info = power_spectrum.info
                kmode_each = algebra.make_vect(power_2d_raw_kmn_kpp[i],
                        axis_names=power_spectrum.info['axes'])
                kmode_each.info = power_spectrum.info
                power_spectrum_kpp_list.append(power_each)
                k_mode_number_kpp_list.append(kmode_each)

                power_each = algebra.make_vect(power_2d_raw_kpn[i],
                        axis_names=power_spectrum.info['axes'])
                power_each.info = power_spectrum.info
                kmode_each = algebra.make_vect(power_2d_raw_kmn_kpn[i],
                        axis_names=power_spectrum.info['axes'])
                kmode_each.info = power_spectrum.info
                power_spectrum_kpn_list.append(-power_each)
                k_mode_number_kpn_list.append(kmode_each)


            sn_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_list, k_mode_number_list, None, truncate_range)

            positive, negative = seperate_positive_and_negative_power(ps_1d_mean, 
                                                                      ps_1d_std, 
                                                                      k_centre)
            #temp_positive, temp_negative = seperate_positive_and_negative_power(temp_1d_mean,
            #                                                          temp_1d_std,
            #                                                          k_centre)

            filename = '%scros_ps_%dmode_1dpow_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

            filename = '%scros_ps_%dmode_1dpow_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

            #filename = '%stemp_ps_%dmode_1dpow_positive.txt'%(pre, mode)
            #np.savetxt('%s/%s'%(outputroot, filename), temp_positive.T, fmt='%15.12e')

            #filename = '%stemp_ps_%dmode_1dpow_negative.txt'%(pre, mode)
            #np.savetxt('%s/%s'%(outputroot, filename), temp_negative.T, fmt='%15.12e')

            # for quadrupole
            sn_1d_mean, qp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_list, k_mode_number_list, 
                    None, truncate_range, order=2)

            positive, negative = seperate_positive_and_negative_power(qp_1d_mean, 
                                                                      qp_1d_std, 
                                                                      k_centre)

            filename = '%scros_ps_%dmode_quadrupole_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

            filename = '%scros_ps_%dmode_quadrupole_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

            # for dipole
            sn_1d_mean, dp_kpp_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_kpp_list, k_mode_number_kpp_list, 
                    None, truncate_range, order=1)
            sn_1d_mean, dp_kpn_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_kpn_list, k_mode_number_kpn_list, 
                    None, truncate_range, order=1)
            dp_std = np.sqrt(0.5 * (dp_kpp_std**2 + dp_kpn_std**2))
            positive, negative = seperate_positive_and_negative_power(dp_mean, 
                                                                      dp_std, 
                                                                      k_centre)
            filename = '%scros_ps_%dmode_dipole_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

            filename = '%scros_ps_%dmode_dipole_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')


        # noise inv weight
        if var != None:
            pre += 'varweight_'
          
            # get the 1d power
            ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum, k_mode_number, var, truncate_range)
            #print '2d power = ', power_spectrum, power_spectrum.shape, '1d power = ', ps_1d_mean, ps_1d_mean.shape
            # get the 1d quadrupole 
            qp_1d_mean, qp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum, k_mode_number, var, truncate_range, order=2)
            # get dipole
            dp_kpp_mean, dp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    ps_2d_kpp, k_mode_number, var, truncate_range, order=1)
            dp_kpn_mean, dp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    -ps_2d_kpn, k_mode_number, var, truncate_range, order=1)
            dp_mean = 0.5 * (dp_kpp_mean + dp_kpn_mean)


            # get the 1d error
            power_spectrum_list = []
            k_mode_number_list = []
            power_spectrum_kpp_list = []
            k_mode_number_kpp_list = []
            power_spectrum_kpn_list = []
            k_mode_number_kpn_list = []
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

                power_each = algebra.make_vect(power_2d_raw_kpp[i],
                        axis_names=power_spectrum.info['axes'])
                power_each.info = power_spectrum.info
                kmode_each = algebra.make_vect(power_2d_raw_kmn_kpp[i],
                        axis_names=power_spectrum.info['axes'])
                kmode_each.info = power_spectrum.info
                power_spectrum_kpp_list.append(power_each)
                k_mode_number_kpp_list.append(kmode_each)

                power_each = algebra.make_vect(power_2d_raw_kpn[i],
                        axis_names=power_spectrum.info['axes'])
                power_each.info = power_spectrum.info
                kmode_each = algebra.make_vect(power_2d_raw_kmn_kpn[i],
                        axis_names=power_spectrum.info['axes'])
                kmode_each.info = power_spectrum.info
                power_spectrum_kpn_list.append(-power_each)
                k_mode_number_kpn_list.append(kmode_each)


            sn_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_list, k_mode_number_list, var, truncate_range)

            positive, negative = seperate_positive_and_negative_power(ps_1d_mean, 
                                                                      ps_1d_std, 
                                                                      k_centre)
            #print 'real positive = ', positive, positive.shape
            filename = '%scros_ps_%dmode_1dpow_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')
            filename = '%scros_ps_%dmode_1dpow_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

            # for quadrupole
            sn_1d_mean, qp_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_list, k_mode_number_list, 
                    var, truncate_range, order=2)

            positive, negative = seperate_positive_and_negative_power(qp_1d_mean, 
                                                                      qp_1d_std, 
                                                                      k_centre)

            filename = '%scros_ps_%dmode_quadrupole_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

            filename = '%scros_ps_%dmode_quadrupole_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')
            #print 'real quad = ', positive, positive.shape

            # for dipole
            sn_1d_mean, dp_kpp_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_kpp_list, k_mode_number_kpp_list, 
                    var, truncate_range, order=1)
            sn_1d_mean, dp_kpn_std, k_centre = pss.convert_2dps_to_1dps_each(
                    power_spectrum_kpn_list, k_mode_number_kpn_list, 
                    var, truncate_range, order=1)
            dp_std = np.sqrt(0.5 * (dp_kpp_std**2 + dp_kpn_std**2))
            positive, negative = seperate_positive_and_negative_power(dp_mean, 
                                                                      dp_std, 
                                                                      k_centre)

            filename = '%scros_ps_%dmode_dipole_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

            filename = '%scros_ps_%dmode_dipole_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')


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

        print 'Attempting to open ', inputroot + 'ps_result.hd5', ' from directory ', os.getcwd()
        ps_result = h5py.File(inputroot + 'ps_result.hd5', 'r')

        rf_root = 'auto_rf_%dmode_2dpow'%mode
        tr_root = 'auto_tr_%dmode_2dpow'%mode

        ps_root = 'auto_ps_%dmode_2dpow'%mode
        ne_root = 'auto_ne_%dmode_2dpow'%mode
        ns_root = 'auto_ns_%dmode_2dpow'%mode
        
        k_p_edges = None
        k_v_edges = None

        # load the transfer function
        if rf_root in ps_result and tr_root in ps_result:
            
            power_rf, k_mode_number, k_p_edges, k_v_edges =\
                    pss.load_power_spectrum(rf_root, ps_result)

            transfer_function = \
                    pss.load_transfer_function(rf_root, tr_root, ps_result, True)[0]

            transfer_function[transfer_function==0] = np.inf
           
            plot_2d_power_spectrum(1./transfer_function, None, k_p_edges, k_v_edges,
                                   filename = 'auto_transferf_%dmode'%mode,
                                   label_list='', output=outputroot,
                                   logscale=False, cmax=1., cmin=0.3)

            transfer_function[ np.isinf(transfer_function)] = 0.
        else:
            print "Note: data for transfer function estimation not exists."
            transfer_function = None


        # plot simulation
        si_root = 'auto_si_%dmode_2draw'%mode
        if si_root in ps_result:

            power_spectrum_sim, k_mode_number_sim, k_p_edges, k_v_edges =\
                pss.load_power_spectrum('auto_si_%dmode_2dpow'%mode, ps_result)

            if transfer_function !=None:
                power_spectrum_sim *= transfer_function

            plot_2d_power_spectrum(power_spectrum_sim, k_mode_number_sim, 
                                   k_p_edges, k_v_edges,
                                   filename = 'auto_si_%dmode'%mode,
                                   label_list='', output=outputroot,)

            power_2d_raw_sim = ps_result[si_root].value
            power_2d_raw_kmn_sim = ps_result['auto_si_%dmode_2draw_kmn'%mode].value

            if transfer_function != None:
                power_2d_raw_sim *= transfer_function[None,...]

            sim_list = []
            sim_kmn_list = []
            for i in range(power_2d_raw_sim.shape[0]):
                sim_each = algebra.make_vect(power_2d_raw_sim[i], 
                        axis_names=power_spectrum_sim.info['axes'])
                sim_each.info = power_spectrum_sim.info
                kmn_each = algebra.make_vect(power_2d_raw_kmn_sim[i], 
                        axis_names=power_spectrum_sim.info['axes'])
                kmn_each.info = power_spectrum_sim.info
                sim_list.append(sim_each)
                sim_kmn_list.append(kmn_each)

            std = np.std(np.array(sim_list), axis=0)
            var = std**2
            var[var == 0] = np.inf
            var = var**(-1)
           
            si_1d_mean, si_1d_std, si_k_centre = pss.convert_2dps_to_1dps_each(
                    sim_list, sim_kmn_list, var, None)
            si_positive, si_negative = seperate_positive_and_negative_power(
                    si_1d_mean, si_1d_std, si_k_centre)
            filename = '%sauto_si_%dmode_1dpow_positive.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_positive.T, fmt='%15.12e')
            filename = '%sauto_si_%dmode_1dpow_negative.txt'%(pre, mode)
            np.savetxt('%s/%s'%(outputroot, filename), si_negative.T, fmt='%15.12e')

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
                ps_root_sec = ps_root + '_%s'%(sec1 + sec2 + 'x' + sec2 + sec1)
                power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
                    pss.load_power_spectrum(ps_root_sec, ps_result)
                power_spectrum_list.append(power_spectrum)
                k_mode_number_list.append(k_mode_number)
                label_list.append(sec1 + sec2)

        #ps_root_sec = ps_root + '_%s'%(sec1 + sec2)
        #power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
        #        pss.load_power_spectrum(ps_root, ps_result)
        #power_spectrum_list.append(power_spectrum)
        #k_mode_number_list.append(k_mode_number)

        plot_2d_power_spectrum(power_spectrum_list, None, 
                               k_p_edges, k_v_edges, filename='auto_ps_%dmode'%mode,                         
                               label_list=label_list, output=outputroot)

        # load the noise weight
        #weight = pss.load_weight_sec(inputroot + 'auto_ns_%dmode_2dpow'%mode, None)
        #plot_2d_power_spectrum(weight, None, k_p_edges, k_v_edges,
        #                       filename = 'auto_wt_%dmode'%mode,
        #                       label='', output=outputroot)

        # truncate the 2d power
        if self.params['truncate'] != []:
            pre += 'trun_'
            truncate_range = self.params['truncate']
        else:
            truncate_range = None

        # compensate the power
        if transfer_function != None:
            pre += 'comp_'
            power_spectrum_list = [ ps * transfer_function for ps in power_spectrum_list]
            transfer_function[transfer_function==0] = np.inf
            #weight /= transfer_function ** 2
            transfer_function[np.isinf(transfer_function)] = 0.

            plot_2d_power_spectrum(power_spectrum_list, None, k_p_edges, k_v_edges, 
                                   filename='comp_auto_ps_%dmode'%mode, 
                                   label_list=label_list, output=outputroot)
            plot_2d_power_spectrum(var, None, k_p_edges, k_v_edges,
                                   filename = 'comp_auto_wt_%dmode'%mode,
                                   label_list='', output=outputroot)

        # get the 1d power
        ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
                power_spectrum_list, k_mode_number_list, var, truncate_range)

        positive, negative = seperate_positive_and_negative_power(ps_1d_mean, 
                                                                  ps_1d_std, 
                                                                  k_centre)

        filename = '%sauto_ps_%dmode_1dpow_positive.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')

        filename = '%sauto_ps_%dmode_1dpow_negative.txt'%(pre, mode)
        np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

if __name__=="__main__":
    #plot_2d_dipole()
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
    #GBTxWiggleZPowerSpectrum_Analysis(autops_init).execute(1)
