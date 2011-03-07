"""
Proceedure to estimate a GBT noise model from data.  This module is now
deprecated and will not work with new support code.
"""
import os

import scipy as sp
import scipy.fftpack as fft
import numpy.ma as ma
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT

# Define a dictionary with keys the names of parameters to be read from
# file and values the defaults.
params_init = {
               # Input and output.
               "input_root" : "./",
               "file_middles" : ("GBTdata",),
               "input_end" : ".raw.acs.fits",
               "output_root" : "noisemodel",
               "output_end" : '',
               # Select data to process.
               "scans" : (),
               "IFs" : (),
               # Algorithm
               "calculate_power_spectrum" : False,
               "calculate_covariance" : True,
               "subtract_freq_average" : False,
               "lags" : tuple(sp.arange(0.01, 61, 5.)),
               "normalize_to_average" : False,
               "norm_pol_weights" : ( 1., 0., 0., 1.),
               "norm_cal_weights" : ( 1., 1.),
               "normalize_dnudt" : True,
               "segment_length" : 0,
               # Polarizations and Cal States
               "pol_weights" : ( 1., 0., 0., 1.),
               "cal_weights" : (0., 1.),
               }
prefix = 'np_'


class NoisePower(object) :
    """Calculates time power spectrum and correlation function of data.
    """
    
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                          prefix=prefix)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
        
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        cal_weights = params['cal_weights']
        pol_weights = params['pol_weights']
        scan_len = params['segment_length']
        
        first_iteration = True
        # Loop over files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            # Read in the data, and loop over data blocks.
            Reader = core.fitsGBT.Reader(input_fname)
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)
            # Loop over scans.
            for Data in Blocks :
                if (Data.dims[1] != len(pol_weights) or
                    Data.dims[2] != len(cal_weights)) :
                    raise ValueError("pol_wieght or cal_weight parameter "
                                     "dimensions don' match data dimensions.")
                if scan_len == 0 :
                    scan_len = Data.dims[0]
                if scan_len > Data.dims[0] :
                    print "Warning: asked for segment length longer than scan."
                # Get the desired combination fo polarizations and cal states.
                data = ma.zeros((scan_len, Data.dims[-1]))
                for ii, pol_w in enumerate(pol_weights) :
                    for jj, cal_w in enumerate(cal_weights) :
                        data[:scan_len] += (Data.data[:scan_len,ii,jj,:]
                                            *pol_w*cal_w)
                # Calculate the time axis.
                Data.calc_time()
                time = Data.time[:scan_len]
                n_time = scan_len
                dt = abs(sp.mean(sp.diff(time)))
                tmean = ma.mean(data, 0)
                data -= tmean
                if params["normalize_dnudt"] :
                    dnu = abs(Data.field["CDELT1"])
                    data *= sp.sqrt(dnu*dt)
                if params['normalize_to_average'] :
                    # Get normalization for each frequency.
                    tmp_mean = ma.mean(Data.data, 0)
                    norm = sp.zeros((Data.dims[-1],))
                    for ii, pol_w in enumerate(params["norm_pol_weights"]) :
                        for jj, cal_w in enumerate(params["norm_cal_weights"]) :
                            norm += tmp_mean[ii,jj,:] *pol_w*cal_w
                    data /= norm
                if params['subtract_freq_average'] :
                    data -= ma.mean(data, 1)[:,sp.newaxis]
                if params["calculate_covariance"] :
                    if first_iteration :
                        lags = params['lags']
                        nlags = len(lags)
                        covariance = sp.zeros(nlags, Data.dims[-1])
                    # Loop over times to accumulate covariances
                    for ii in range(ntime) :
                        # Take products of data not yet considered and find
                        # thier lag.
                        squares = data[ii,:]*data[ii:,:]
                        # Loop over lags.
                        for jj in range(nlags) :
                            if jj == 0 :
                                (inds,) = sp.where(sp.logical_and(
                                                   dt <= lags[jj], dt >= 0))
                            else :
                                (inds,) = sp.where(sp.logical_and(
                                            dt <= lags[jj], dt > lags[jj-1]))
                            if len(inds) > 0 :
                                tempvar = ma.sum(squares[inds,:], 0)
                                covariance[jj,:] += tempvar.filled(0)
                            counts[jj,:] += squares[inds,:].count(0)
                        # End ii loop over times
                    # End calculate_covariance if
                if params["calculate_power_spectrum"] :
                    # For now just assume all the blocks are the same
                    # length (which will be the case in practice).  This
                    # will get more sophisticated eventually.
                    if first_iteration :
                        # Allowcate memory and calculate the FT time axis.
                        power_spectrum = sp.zeros((n_time//2 + 1, 
                                                   Data.dims[-1]))
                        power_counts = sp.zeros(Data.dims[-1], dtype=int)
                        ps_freqs = sp.arange(n_time//2 + 1, dtype=float)
                        ps_freqs /= (n_time//2 + 1)*dt*2
                    # Throw frequencies with masked data.
                    good_freqs = ma.count_masked(data, 0)==0
                    # Calculate power spectrum.
                    this_power = abs(fft.fft(data, axis=0)
                                     [range(n_time//2+1)])
                    this_power = this_power**2/n_time
                    power_spectrum[:, good_freqs] += this_power[:, good_freqs]
                    power_counts[good_freqs] += 1
                    # End power spectrum if
                first_iteration = False
                # End loop over files scans.
            # End loop over files
        # Store outputs in the class.
        if params["calculate_covariance"] :
            self.covariance = covariance/counts
            self.covariance_counts = counts
            self.lags = lags
        if params["calculate_power_spectrum"] :
            self.power_spectrum = power_spectrum/power_counts
            self.power_counts = power_counts
            self.ps_freqs = ps_freqs


def make_plots(data, ini) :
    """
    Make plots for noise model.

    This function takes the data from noise_model.noise_model and plots it.  It
    can then either show the plots, or save them to file.

    Arguments:
        data : A dictionary of data that is the output of noise_model() (or
               compatable).
        ini : Either a dictionary or file name containing the parameters
              for this program.  This is parsed using the kiyopy.parser.
              All parameters begin with the prefix plot_ and so common
              practice should be to include the parameters for this
              program in the parameter file for noise_model()
    """
    
    # Get input parameters.
    params_init = {
                      "plot_save_file" : False,
                      "plot_norm_to_first_lag" : False,
                      "plot_fbins" : [],
                      "plot_show" : False,
                      "plot_output_root" : "",
                      "plot_output_dir" : "./"
                      }
    parameters = parse_ini.parse(ini, params_init)
    # Unpack data.
    block_LST = data["block_LST"]
    block_medians = data["block_medians"]
    block_RMS = data["block_RMS"]
    if data.has_key("covariance") :
        covariance = data["covariance"]
        lags = data["lags"]
        nlags = len(lags)
        # Get center bins for plotting
        plot_lags = [lags[0]/2]
        for ii in range(1,nlags) :
            plot_lags.append((lags[ii]+lags[ii-1])/2)
        plot_lags = sp.array(plot_lags)
        if parameters["plot_norm_to_first_lag"] :
            covariance /= covariance[0,:]
    if data.has_key("power_spectrum") :
        power_spectrum = data["power_spectrum"]
        ps_freqs = data["ps_freqs"]
    # Isolate only the frequency bins we need for plotting. 
    fbins_plot = parameters["plot_fbins"]
    nfbins_plot = len(fbins_plot)
    if nfbins_plot == 0 :
        nfbins_plot = len(block_medians[0,:])
        fbins_plot = range(nfbins_plot)
    plot_block_RMS = block_RMS[:,fbins_plot]
    plot_block_medians = block_medians[:,fbins_plot]
    markers = ['o','s','^','h','p','v']
    colours = ['b','g','r','c','m','y']
    nmarkers = len(markers)
    ncolours = len(colours)

    plt.figure(1)
    for ii in range(nfbins_plot) :
        plt.plot(block_LST, plot_block_medians[:,ii],
                 linestyle='None',
                 marker=markers[ii%nmarkers],
                 markerfacecolor=colours[ii%ncolours])
    plt.xlabel("File Start LST")
    plt.ylabel("File Time Median")

    plt.figure(2)
    for ii in range(nfbins_plot) :
        plt.plot(block_LST, plot_block_RMS[:,ii],
                 linestyle='None',
                 marker=markers[ii%nmarkers],
                 markerfacecolor=colours[ii%ncolours])
    plt.xlabel("File Start LST")
    plt.ylabel("File RMS")
    
    if data.has_key("covariance") :
        plt.figure(3)
        for ii in range(nfbins_plot) :
            inds, = sp.where(covariance[:,fbins_plot[ii]] > 0)
            a = sp.array(plot_lags[inds])
            b = sp.array(covariance[inds,fbins_plot[ii]])
            plt.loglog(plot_lags[inds], covariance[inds,fbins_plot[ii]],
                       linestyle='None',
                       marker=markers[ii%nmarkers],
                       markerfacecolor=colours[ii%ncolours])
            inds, = sp.where(covariance[:,fbins_plot[ii]] < 0)
            plt.loglog(plot_lags[inds], -covariance[inds,fbins_plot[ii]],
                       linestyle='None',
                       marker=markers[ii%nmarkers],
                       markeredgecolor=colours[ii%ncolours],
                       markerfacecolor='w')
        plt.xlabel("Lag (sidereal seconds)")
        plt.ylabel("covariance")

    if data.has_key("power_spectrum") :
        plt.figure(4)
        for ii in range(nfbins_plot) :
            plt.loglog(ps_freqs, power_spectrum[:,fbins_plot[ii]],
                       linestyle='None',
                       marker=markers[ii%nmarkers],
                       markerfacecolor=colours[ii%ncolours])
        plt.xlabel("frequency (sidereal Hz)")
        plt.ylabel("power")
        

    if parameters['plot_save_file'] :
        os.system('mkdir -p ' + parameters["plot_output_dir"])
        out_root = (parameters["plot_output_dir"] + '/' +
                   parameters["plot_output_root"])
        plt.figure(1)
        plt.savefig(out_root+"medians.eps")
        plt.figure(2)
        plt.savefig(out_root+"RMS.eps")
        if data.has_key("covariance") :
            plt.figure(3)
            plt.savefig(out_root+"covariance.eps")
        if data.has_key("power_spectrum") :
            plt.figure(4)
            plt.savefig(out_root+"power_spectrum.eps")

    if parameters['plot_show'] :
        plt.show()


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    NoisePower(str(sys.argv[1])).execute()

