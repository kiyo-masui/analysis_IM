"""Proceedure to estimate a GBT noise model from data."""

import os

import scipy as sp
import scipy.fftpack as fft
import numpy.ma as ma
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
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
               # Whether to stack the results for all the files togther or to
               # hold everything in a list.
               "stack_files" : True
               }
prefix = 'np_'


class NoisePower(object) :
    """Calculates time power spectrum and correlation function of data.
    """
    
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                      prefix=prefix, feedback=feedback)
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

def make_masked_time_stream(Blocks) :
    """Converts Data Blocks into a single uniformly sampled time stream.
    
    Also produces the mask giving whether elements are valid entries or came
    from a zero pad.  This produes the required inputs for calculating a
    windowed power spectrum.

    Parameters
    ----------
    Blocks : tuple of DataBlock objects.

    Returns
    -------
    time_stream : array
        All the data in `Blocks` but concatenated along the time axis and
        padded with zeros such that the time axis is uniformly sampled and
        uninterupted.
    mask : array same shape as `time_stream`
        1.0 if data in the correspoonding `time_stream` element is filled 
        and 0 if the data was missing.  This is like a window where 
        time_stream = mask*real_data.
    """

    # Shape of all axes except the time axis.
    back_shape = Blocks[0].dims[1:]
    # Get the time sample spacing.
    Blocks[0].calc_time()
    dt = abs(sp.mean(sp.diff(Blocks[0].time)))
    # Find the beginning and the end of the time axis by looping through
    # blocks.
    min_time = float('inf')
    max_time = 0.0
    for Data in Blocks :
        Data.calc_time()
        min_time = min(min_time, min(Data.time))
        max_time = max(min_time, max(Data.time))
        # Ensure that the time sampling is uniform.
        if not (sp.allclose(abs(sp.diff(Data.time)), dt, rtol=0.1)
                and sp.allclose(abs(sp.mean(sp.diff(Data.time))), dt)) :
            msg = ("Time sampling not uniformly spaced or Data Blocks don't "
                   "agree on sampling.")
            raise ce.DataError(msg)
        # Ensure the shapes are right.
        if Data.dims[1:] != back_shape :
            msg = ("All data blocks must have the same shape except the time "
                   "axis.")
            raise ce.DataError(msg)
    # Calculate the time axis.
    time = sp.arange(min_time, max_time + dt, dt)
    ntime = len(time)
    # Allowcate memory for the outputs.
    time_stream = sp.zeros((ntime,) + back_shape, dtype=float)
    mask = sp.zeros((ntime,) + back_shape, dtype=sp.float32)
    # Loop over all times and fill in the arrays.
    for Data in Blocks :
        # Apply an offset to the time in case the start of the Data Block
        # doesn't line up with the time array perfectly.
        offset = time[sp.argmin(abs(time - Data.time[0]))] - Data.time[0]
        for ii in range(Data.dims[0]) :
            ind = sp.argmin(abs(time - (Data.time[ii] + offset)))
            if sp.any(mask[ind, ...]) :
                msg = "Overlapping times in Data Blocks."
                raise ce.DataError(msg)
            time_stream[ind, ...] = Data.data[ii, ...].filled(0.0)
            mask[ind, ...] = sp.logical_not(Data.data.mask[ii, ...])
    return time_stream, mask

def windowed_power(data1, window1, data2=None, window2=None) :
    """Calculates a windowed cross power spectrum.

    Calculates the cross power spectrum of uniformly and continuously
    sampled time stream data.  There may be missing data samples in where both
    the data and the window should be 0, otherwise the window should be unity.
    
    All input arrays must be 1D and the same length.

    Parameters
    ----------
    data1 : array
        First data time stream.
    window1 : array
        Window for first time stream.
    data2 : array
        Second data time stream. If None, set equal to `data1`.
    window2 : array
        Window for second time stream. If None, set equal to `window1`.  
        Ignored if data2 is None.

    Returns
    -------
    power_spectrum : 1D array
        Cross power of the input time streams accounting for the windows.
    
    Notes
    -----
    This function calculates the best estimate for the cross power of
    real_data1 and real_data2 given data1 = window1*real_data1 and 
    data2 = window2*real_data2.  The most common case will be when the windows
    are simply masks: arrays of 0s and 1s, but the function should work fien
    for more general windows.
    """

    if data2 is None :
        data2 = data1
        window2 = window1

    data_power = fft.fft(data1)*fft.fft(data2).conj()
    window_power = fft.fft(window1)*fft.fft(window2).conj()

    true_corr = fft.ifft(data_power)/fft.ifft(window_power)
    true_power = fft.fft(true_corr)

    return true_power



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    NoisePower(str(sys.argv[1])).execute()

