"""Proceedure to estimate a GBT noise model from data."""

import os
import sys

import scipy as sp
import scipy.fftpack as fft
import scipy.signal as sig
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
    # Polarizations and Cal States
    "pol_weights" : ( 1.,),
    "cal_weights" : (1.,),
    # Whether to put into units of the thermal expectation.
    "norm_to_thermal" : False,
    # How many time bins to use.  By default, match to the first file.
    "n_time_bins" : 0
    }

class NoisePower(object) :
    """Calculates the Noise power spectrum and other noise statistics."""
  
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
        n_time = params["n_time_bins"]
        n_files = len(params["file_middles"])
        first_iteration = True
        # Loop over files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            # Read in the data, and loop over data blocks.
            Reader = core.fitsGBT.Reader(input_fname)
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)
            for Data in Blocks:
                data = Data.data
                data_selected = ma.zeros((Data.dims[0], 1, 1, Data.dims[3]),
                                          dtype=float)
                data_selected.mask = ma.getmaskarray(data_selected)
                for ii in range(len(pol_weights)) :
                    for jj in range(len(cal_weights)) :
                        data_selected[:,0,0,:] += (data[:,ii,jj,:]
                                                   * pol_weights[ii]
                                                   * cal_weights[jj])
                Data.set_data(data_selected)
            full_data, mask, this_dt = make_masked_time_stream(Blocks, n_time)
            # TODO: Figure out a better way to deal with this (only drop
            # affected frequencies).
            if sp.any(sp.allclose(mask[:,0,0,:], 0.0, 0)):
                n_files -= 1
                continue
            if first_iteration and not n_time :
                n_time = full_data.shape[0]
                dt = this_dt
            elif abs((this_dt - dt)/dt) > 0.001 :
                msg = "Files have different time samplings."
                raise ce.DataError(msg)

            # XXX to speed things up for testing.
            full_data = full_data[...,0::10]
            mask = mask[...,0::10]
            #full_data = full_data[...,[80]]
            #mask = mask[...,[80]]
            # XXX

            full_mean = sp.mean(full_data, 0)/sp.mean(mask, 0)
            full_data -= full_mean
            full_mean.shape = full_mean.shape[-1:]
            full_data *= mask
            full_data1 = full_data[:,0,0,:,None]
            full_data2 = full_data[:,0,0,None,:]
            mask1 = mask[:,0,0,:,None]
            mask2 = mask[:,0,0,None,:]

            n_freq = full_data.shape[-1]
            thermal_expectation = (full_mean / sp.sqrt(dt) 
                                / sp.sqrt(abs(Blocks[0].field['CDELT1'])))
            # Loop to do only a bit of the calculation at a time.  Reduces
            # memory use by a factor of a few.
            power_mat = sp.empty((n_time//2, n_freq, n_freq),
                                      dtype=float)
            for ii in xrange(n_freq):
                # `out` argument must be a view, not a copy.
                windowed_power(full_data1[:,[ii],:], mask1[:,[ii],:],
                               full_data2, mask2, axis=0,
                               out=(power_mat[:,ii,:])[:,None,:])
            frequency = ps_freq_axis(dt, full_data.shape[0])
            if params['norm_to_thermal'] :
                power_mat /= (thermal_expectation[:,None] 
                              * thermal_expectation[None,:])
            # Combine across files.
            if first_iteration :
                self.power_mat = power_mat
                self.thermal_expectation = thermal_expectation
                self.frequency = frequency
            else :
                self.power_mat += power_mat
                self.thermal_expectation += thermal_expectation
            first_iteration = False
        if n_files > 0 :
            self.power_mat /= n_files
            self.thermal_expectation / n_files





params_init_old = {
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


class NoisePowerOld(object) :
    """Calculates time power spectrum and correlation function of data.
    """
    
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init_old, 
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

def ps_freq_axis(dt, n):
    """Calculate the frequency axis for a power spectrum.

    Parameters
    ----------
    dt : float
        time step for time stream that went into power spectrum calculation.
    n : int
        Number of time stream data points that went into the power spectrum
        calculation.

    returns
    -------
    frequencies : 1D array of floats.
        The frequency axis.  Has length `n`//2.
    """

    frequency = sp.arange(n//2, dtype=float)
    df = 1.0/dt/n
    frequency *= df
    return frequency

def make_masked_time_stream(Blocks, ntime=None, window=None) :
    """Converts Data Blocks into a single uniformly sampled time stream.
    
    Also produces the mask giving whether elements are valid entries or came
    from a zero pad.  This produes the required inputs for calculating a
    windowed power spectrum.

    Parameters
    ----------
    Blocks : tuple of DataBlock objects.
    ntime : int
        Total number of time bins in output arrays.  If shorter than required
        extra data is truncated.  If longer, extra data is masked.  Default is
        to use exactly the number that fits all the data.
    window : string or tuple
        Type of window to apply to each DataBlock.  Valid options are the valid
        arguments to scipy.signal.get_window().  By default, don't window.

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
    dt : float
        The time step of the returned time stream.
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
                and sp.allclose(abs(sp.mean(sp.diff(Data.time))), dt,
                                rtol=0.001)) :
            msg = ("Time sampling not uniformly spaced or Data Blocks don't "
                   "agree on sampling.")
            raise ce.DataError(msg)
        # Ensure the shapes are right.
        if Data.dims[1:] != back_shape :
            msg = ("All data blocks must have the same shape except the time "
                   "axis.")
            raise ce.DataError(msg)
    # Calculate the time axis.
    if not ntime :
        time = sp.arange(min_time, max_time + dt, dt)
        ntime = len(time)
    else :
        time = sp.arange(ntime)*dt + min_time
    # Allowcate memory for the outputs.
    time_stream = sp.zeros((ntime,) + back_shape, dtype=float)
    mask = sp.zeros((ntime,) + back_shape, dtype=sp.float32)
    # Very important to subtract the mean out of the signal, otherwise the
    # window coupling to the mean (0) mode will dominate everything.
    total_sum = 0
    total_counts = 0
    for Data in Blocks:
        total_sum += sp.sum(Data.data.filled(0), 0)
        total_counts += ma.count(Data.data, 0)
    total_mean = total_sum / total_counts
    # Loop over all times and fill in the arrays.
    for Data in Blocks :
        # Subtract the mean calculated above.
        Data.data -= total_mean
        # Apply an offset to the time in case the start of the Data Block
        # doesn't line up with the time array perfectly.
        offset = time[sp.argmin(abs(time - Data.time[0]))] - Data.time[0]
        # Generate window function.
        if window:
            window_function = sig.get_window(window, Data.dims[0])
        for ii in range(Data.dims[0]) :
            ind = sp.argmin(abs(time - (Data.time[ii] + offset)))
            if abs(time[ind] - (Data.time[ii])) < 0.5*dt :
                if sp.any(mask[ind, ...]) :
                    msg = "Overlapping times in Data Blocks."
                    raise ce.DataError(msg)
                if window:
                    window_value = window_function[ii]
                else :
                    window_value = 1.0
                time_stream[ind, ...] = (window_value 
                                         * Data.data[ii, ...].filled(0.0))
                mask[ind, ...] = window_value * sp.logical_not(ma.getmaskarray(
                                     Data.data)[ii, ...])
    return time_stream, mask, dt

def windowed_power(data1, window1, data2=None, window2=None, axis=-1, out=None) :
    """Calculates a windowed cross power spectrum.

    Calculates the cross power spectrum of uniformly and continuously
    sampled time stream data.  There may be missing data samples in where both
    the data and the window should be 0, otherwise the window should be unity.
    
    All input arrays must be the same shape or broadcastable to the same
    shape.

    Parameters
    ----------
    data1 : array
        First data time stream.
    window1 : array
        Window for first time stream.
    data2 : array, default `data1`
        Second data time stream. If None, set equal to `data1`.
    window2 : array, default `window1`
        Window for second time stream. If None, set equal to `window1`.  
        Ignored if data2 is None.
    axis : int
        For multidimentional inputs, the axis along which to calculate the
        power spectrum.
    out : array
        Optional preallowcated output array for storing result.  Must be
        compatible with return value.

    Returns
    -------
    power_spectrum : 1D array
        Cross power of the input time streams accounting for the windows.  It
        will be the same shape as the input arrays accept for the `axis` axis
        will be shortened to n//2.
    
    Notes
    -----
    This function calculates the best estimate for the cross power of
    real_data1 and real_data2 given data1 = window1*real_data1 and 
    data2 = window2*real_data2.  The most common case will be when the windows
    are simply masks: arrays of 0s and 1s, but the function should work fine
    for more general windows.
    """

    if data2 is None :
        data2 = data1
    if window2 is None :
        window2 = window1
    
    # Delete variables ASAP as they are memory intensive.
    data_power = fft.fft(data1, axis=axis)*fft.fft(data2, axis=axis).conj()
    window_power = (fft.fft(window1, axis=axis)
                    * fft.fft(window2, axis=axis).conj())
    data_corr = fft.ifft(data_power, axis=axis, overwrite_x=True)
    window_corr = fft.ifft(window_power, axis=axis, overwrite_x=True)
    true_corr = data_corr / window_corr
    n = data1.shape[axis]
    # XXX
    #plt.figure()
    #plt.semilogy(window_power[:,0,0], '.')
    #plt.semilogy(-window_power[:,0,0], '.r')
    #plt.figure()
    #plt.semilogy(window_corr[:, 0, 0])
    #plt.semilogy(data_corr[:, 0, 0])
    #plt.semilogy(true_corr[:, 0, 0])
    #plt.semilogy(true_corr[:, 0, 1])
    #plt.show()
    #print sp.where(window_corr == 0)
    # XXX
    del data_power, window_power, data_corr, window_corr
    true_power = fft.fft(true_corr, axis=axis, overwrite_x=True)
    del true_corr
    # Truncate the power spectrum to only the unaliased and positive
    # frequencies.
    s = slice(n//2)
    indices = [slice(sys.maxsize)] * data1.ndim
    indices[axis] = s
    indices = tuple(indices)
    if not out is None :
        out[...] = true_power[indices]
    else :
        # Copy to release extra memory.
        out = sp.copy(true_power[indices])
    return out

def full_power_mat(Blocks, n_time=None, window=None, deconvolve=True) :
    """Calculate the full power spectrum of a data set with channel
    correlations.
    
    Only one cal state and pol state assumed.
    """

    full_data, mask, dt = make_masked_time_stream(Blocks, n_time, window)
    n_time = full_data.shape[0]
    # Broadcast to shape such that all pairs of channels are calculated.
    full_data1 = full_data[:,0,0,:,None]
    full_data2 = full_data[:,0,0,None,:]
    mask1 = mask[:,0,0,:,None]
    mask2 = mask[:,0,0,None,:]
    # XXX Accutally do this efficiently.
    if not deconvolve:
        mask1[:] = 1
        mask2[:] = 1
    # Frequencies for the power spectrum (Hz)
    frequency = ps_freq_axis(dt, n_time)
    n_chan = full_data.shape[-1]
    # Loop to do only a bit of the calculation at a time.  Reduces
    # memory use by a factor of a few.
    power_mat = sp.empty((n_time//2, n_chan, n_chan),
                              dtype=float)
    for ii in xrange(n_chan):
        # `out` argument must be a view, not a copy.
        windowed_power(full_data1[:,[ii],:], mask1[:,[ii],:], full_data2,
                       mask2, axis=0, out=(power_mat[:,ii,:])[:,None,:])

    return power_mat, frequency

def diag_power_mat(Blocks, n_time=None, window=None, deconvolve=True) :
    """Calculate the full power spectrum of a data set without channel
    correlations.
    
    Only one cal state and pol state assumed.
    """

    full_data, mask, dt = make_masked_time_stream(Blocks, n_time, window)
    n_time = full_data.shape[0]
    # Broadcast to shape such that all pairs of channels are calculated.
    full_data1 = full_data[:,0,0,:]
    mask1 = mask[:,0,0,:]
    # XXX Accutally do this efficiently.
    if not deconvolve:
        mask1[:] = 1
        mask2[:] = 1
    # Frequencies for the power spectrum (Hz)
    frequency = ps_freq_axis(dt, n_time)
    n_chan = full_data.shape[-1]
    # Loop to do only a bit of the calculation at a time.  Reduces
    # memory use by a factor of a few.
    power_mat = sp.empty((n_time//2, n_chan),
                              dtype=float)
    windowed_power(full_data1, mask1, axis=0, out=power_mat)

    return power_mat, frequency




# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    NoisePower(str(sys.argv[1])).execute()

