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
import fitsGBT

def noise_model(ini_fname=None) :
    """
    Main function for solving for a noise model from GBT Spectrometer data.
    Depricated and won't work right now.

    Should be passed an input parameter file name.
    """
    
    # Define a dictionary with keys the names of parameters to be read from
    # file and values the defaults.
    params_init = {
                 # Input and output.
                 "input_root" : "./",
                 "input_middles" : ["GBTdata"],
                 "input_ends" : ".raw.acs.fits",
                 "output_dir" : "./",
                 "output_root" : "noisemodel",
                 "make_figures" : False,
                 # Select data to process.
                 "files_scans_per_block" : 1,
                 "scans" : [],
                 "IFs" : [],
                 "first_freq" : 128,
                 "nfreqs":1792,
                 "freqs_per_fbin": 64,
                 # Algorithm
                 "calculate_power_spectrum" : False,
                 "calculate_covariance" : True,
                 "subtract_freq_average" : False,
                 "lags" : list(sp.arange(0.01, 61, 5.)),
                 # Polarizations and Cal States
                 "pol_weights" : [ 1., 0., 0., 1.],
                 "cal_weights" : [0., 1.],
                 "norm_pol_weights" : [ 0., 0., 0., 0.],
                 "norm_cal_weights" : [0., 0.]
                 }
    # Read the parameter file, store in dictionary named parameters.
    parameters, undeclared = parse_ini.parse(ini_fname, params_init,
                                             return_undeclared = True)
    ncal = 2
    npol = 4
    norm_pol_weights = parameters['norm_pol_weights']
    norm_cal_weights = parameters['norm_cal_weights']
    if norm_cal_weights.count(0.) == 2 :
        norm_cal_weights = parameters['cal_weights']
    if norm_pol_weights.count(0.) == 4 :
        norm_pol_weights = parameters['pol_weights']

    # Set up frequencies
    nfreqs = parameters["nfreqs"]
    freqs_per_fbin = parameters["freqs_per_fbin"]
    nfreqs = nfreqs - nfreqs%freqs_per_fbin
    nfbins = nfreqs//freqs_per_fbin
    first_freq = parameters["first_freq"]
    
    out_root = parameters["output_dir"]+'/'+parameters['output_root']
    os.system('mkdir -p ' + parameters["output_dir"])

    # Devide files and scans in data into blocks.  A block will be processes as
    # a continuouse stream of data with the time median subtracted once per
    # block.
    fname_middles = parameters["input_middles"]
    # Figure out how the files and scans are to be grouped
    files_scans_per_block = parameters["files_scans_per_block"]
    scans = parameters["scans"]
    nfiles = len(fname_middles)
    nscans = len(scans)
    if nscans == 0 and files_scans_per_block < 0 :
        raise Exception('If dividing block size by scans'
                        '(files_scans_per_block) then you must explicitly '
                        'give a scan list (scans != []). In parameter file:' +
                        params_fname)
    if nscans == 0 :
        nscans = 1;  # This is a lie, makes the file indexing work.
        all_scans_flag == True # Used for scan indexing.
    else :
        all_scans_flag = False
    if files_scans_per_block == 0 :
        scans_block = nfiles*nscans
    elif files_scans_per_block > 0 :
        scans_block = nscans*files_scans_per_block
    elif files_scans_per_block < 0 :
        scans_block = -files_scans_per_block
    nblocks = nfiles*nscans//scans_block # Excess is truncated, not processed.
    
    # Arrays to store results data 
    block_medians = sp.empty((nblocks,nfbins), dtype=float)
    block_RMS = sp.empty((nblocks,nfbins), dtype=float)
    block_starts_LST = sp.empty(nblocks, dtype=float)
    if parameters["calculate_covariance"] :
        # Set up lags.
        lags = parameters["lags"]
        nlags = len(lags)
        counts = sp.zeros((nlags,nfbins), dtype=int)
        covariance = sp.zeros((nlags,nfbins), dtype=float)
    elif not parameters["calculate_power_spectrum"] :
        raise Exception("Not calculating covariance or power spectrum."
                        " Nothing to do.")
    if (parameters["calculate_power_spectrum"] and not 
        (files_scans_per_block == -1 ) ):
        raise Exception("Power Spectrum calculation only implimented with "
                        "single scan blocks. Set files_scans_per_block = -1.")
    
    first_block_iteration = True
    for kk in range(nblocks) :
        first_file_iteration = True # Reset flag for each block.
        if ((kk+1)*scans_block)%nscans :
            files_this_block = fname_middles[(kk*scans_block)//nscans:
                                             ((kk+1)*scans_block)//nscans+1]
        else :
            files_this_block = fname_middles[(kk*scans_block)//nscans:
                                             ((kk+1)*scans_block)//nscans]
        nfiles_this_block = len(files_this_block)
        scans_sofar = 0 # scans read sofar this block.

        for jj in range(nfiles_this_block) :
            file_middle = files_this_block[jj]
            if all_scans_flag :
                scans_this_file = []
            else :
                scans_this_file = scans[(kk*scans_block + scans_sofar)%nscans :
                                        (kk*scans_block + scans_sofar)%nscans +
                                        scans_block - scans_sofar]
                scans_sofar += len(scans_this_file)
            # Read Fits data from file.  fitsGBT.Processor also does some
            # preliminary processing.
            fname = (parameters['input_root'] + file_middle + 
                        parameters['input_ends'])
            Proc = fitsGBT.Processor(fname)
            data_dict = Proc.read_process(IFs=parameters["IFs"],
                                          scans=scans_this_file)
            # eventually should protect from transiting LST = 0
            LST_this_file = data_dict["LST"]
            print "File LST start: ", LST_this_file[0]
            ntimes = len(LST_this_file)
            # T is an numpy masked array type.  It has internal flags for each
            # data point.
            T_all_freqs = ma.zeros((ntimes,nfreqs))
            T_all_freqs_norm = ma.zeros((ntimes,nfreqs))
            for polindex in range(npol) :
                for calindex in range(ncal) :
                    T_all_freqs +=  ( parameters['cal_weights'][calindex] * 
                        parameters['pol_weights'][polindex] *
                        data_dict["P"][:,polindex,calindex,
                        first_freq:first_freq + nfreqs] )
                    T_all_freqs_norm +=  ( norm_cal_weights[calindex] * 
                        norm_pol_weights[polindex] *
                        data_dict["P"][:,polindex,calindex,
                        first_freq:first_freq + nfreqs] )
           
            # Rebin data frequencies into courser fbins.
            T_this_file = ma.empty((ntimes,nfbins), dtype=float, order='F') 
            T_this_file_norm = ma.empty((ntimes,nfbins), dtype=float,
                                        order='F') 
            # 'F' <- premature optimization.
            for ii in range(ntimes) :
                for jj in range(nfbins) :
                    fbin_inds = sp.arange(jj*freqs_per_fbin, 
                                          (jj+1)*freqs_per_fbin)
                    T_this_file[ii,jj] = ma.median(T_all_freqs[ii,fbin_inds])
                    T_this_file_norm[ii,jj] = ma.median(T_all_freqs_norm
                                                        [ii,fbin_inds])
           
            # Put files together into one big array.
            if first_file_iteration :
                T = ma.array(T_this_file)
                T_norm = ma.array(T_this_file_norm)
                LST = sp.array(LST_this_file)
                block_starts_LST[kk] = LST_this_file[0]
            else :
                T = ma.concatenate((T,T_this_file),0)
                T_norm = ma.concatenate((T_norm,T_this_file_norm),0)
                LST = ma.concatenate((LST,LST_this_file))
            first_file_iteration = False
            # End file loop.

        print 'Taking median, dumping covarience.'
        # subtract off time median
        time_median = ma.median(T, 0)
        time_median_norm = ma.median(T_norm, 0)
        block_medians[kk,:] = time_median
        T = (T - time_median)/time_median_norm
        if parameters['subtract_freq_average'] :
            # T = T/block_RMS[kk,:] now I scale by the median.
            tempmed = ma.median(T, 1)
            T = T - tempmed[:,sp.newaxis]
        block_RMS[kk,:] = sp.sqrt(ma.mean(T**2,0))
        ntimes_block = len(LST)
        
        if parameters["calculate_covariance"] :
            # Loop over times to accumulate covariences
            for ii in range(ntimes_block) :
                # Take products of data not yet considered and find thier lag
                variance = T[ii,:]*T[ii:,:] # This is also a masked array
                dST = abs(LST[ii:]-LST[ii])
                # Loop over lags
                for jj in range(nlags) :
                    if jj == 0 :
                        (inds,) = sp.where(sp.logical_and(dST <= lags[jj], 
                                                          dST >= 0))
                    else :
                        (inds,) = sp.where(sp.logical_and(dST <= lags[jj],
                                                          dST > lags[jj-1]))
                    if len(inds) > 0 :
                        tempvar = ma.sum(variance[inds,:], 0)
                        covariance[jj,:] += tempvar
                    counts[jj,:] += variance[inds,:].count(0)
                # End ii loop over times
            # End calculate_covarience if
            if parameters["calculate_power_spectrum"] :
                # For now just assume all the blocks are the same length (which
                # will be the case in practice).  This will get more
                # sophisticated eventually.
                if first_block_iteration :
                    power_spectrum = sp.zeros((ntimes_block//2 + 1,nfbins))
                    dt = sp.median(sp.diff(LST))
                    ps_freqs = sp.arange(ntimes_block//2 + 1, dtype=float)
                    ps_freqs /= (ntimes_block//2 + 1)*dt*2
                power_spectrum += abs(fft.fft(T, axis=0)[range(ntimes_block//2
                                                         +1)])**2/ntimes_block
            # End calculate_power_spectrum if
            first_block_iteration = False
        # End kk loop over blocks.
    # Gather output data and write outputs.
    out_dict =  {
                 "block_LST" : block_starts_LST,
                 "block_medians" : block_medians,
                 "block_RMS" :block_RMS
                }
    
    # Normalize the covarience and the power spectrum, add them to the outputs.
    if parameters["calculate_covariance"] :
        inds = sp.where(counts<2)
        counts -= 1
        counts[inds] = 1
        covariance /= counts
        out_dict.update( {
                          "lags" : lags,
                          "covariance" : covariance,
                          "counts" : counts
                          } )
    if parameters["calculate_power_spectrum"] :
        power_spectrum /= nblocks
        out_dict.update( {
                          "power_spectrum" : power_spectrum,
                          "ps_freqs" : ps_freqs
                          } )
    
    # Making and saving plots.
    if parameters["make_figures"] :
        undeclared = make_plots(out_dict, undeclared)
    return out_dict



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
    noisemodel(string(sys.argv[1]))

