"""Proceedure to estimate a GBT noise model from data."""

import os
import sys
import math

import scipy as sp
import scipy.fftpack as fft
import scipy.signal as sig
import numpy.random as rand
import numpy.ma as ma
#import matplotlib.pyplot as plt

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
    "n_time_bins" : 0,
    "window" : "hanning",
    "deconvolve" : True,
    "subtract_slope" : False
    }

prefix = 'np_'

class NoisePower(object) :
    """Calculates the Noise power spectrum and other noise statistics."""
  
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                      prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
        
        # Some set up.
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        cal_weights = params['cal_weights']
        pol_weights = params['pol_weights']
        n_time = params["n_time_bins"]
        n_files = len(params["file_middles"])
        window = params['window']
        deconvolve = params['deconvolve']
        subtract_slope = params['subtract_slope']
        # Loop over files to process.
        first_iteration = True
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            # Read in the data.
            Reader = core.fitsGBT.Reader(input_fname)
            Blocks = Reader.read(params['scans'], params['IFs'],
                                 force_tuple=True)
            # Loop over the Blocks to select the channel polarizations and cal
            # state that we want to process.
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
            # Calculate the raw power spectrum.
            power_mat, window_function, this_dt, full_mean = \
                full_power_mat(Blocks, n_time=n_time, window=window, 
                deconvolve=deconvolve, subtract_slope=subtract_slope)
            # Get rid of the extra cal and polarization axes.
            power_mat = power_mat[:,0,0,:,:]
            window_function = window_function[:,0,0,:,:]
            full_mean = full_mean[0,0,:]
            # TODO: Figure out a better way to deal with this (only drop
            # affected frequencies).
            #if sp.any(sp.allclose(mask[:,0,0,:], 0.0, 0)):
            #    n_files -= 1
            #    continue
            if first_iteration :
                n_time = power_mat.shape[0]
                dt = this_dt
            elif not sp.allclose(dt, this_dt, rtol=0.001):
                msg = "Files have different time samplings."
                raise ce.DataError(msg)
            # Format the power spectrum.
            power_mat = prune_power(power_mat, 0)
            power_mat = make_power_physical_units(power_mat, this_dt)
            # TODO In the future the thermal expectation could include
            # polarization factors (be 'I' aware) and cal factors.
            thermal_expectation = (full_mean / sp.sqrt(this_dt)  /
                                   sp.sqrt(abs(Blocks[0].field['CDELT1'])) /
                                   sp.sqrt(1.0 / 2.0 / this_dt))
            if params['norm_to_thermal'] :
                power_mat /= (thermal_expectation[:,None] 
                              * thermal_expectation[None,:])
            # Combine across files.
            if first_iteration :
                self.power_mat = power_mat
                self.thermal_expectation = thermal_expectation
            else :
                self.power_mat += power_mat
                self.thermal_expectation += thermal_expectation
            first_iteration = False
        self.frequency = ps_freq_axis(dt, n_time)
        if n_files > 0 :
            self.power_mat /= n_files
            self.thermal_expectation / n_files
        #plt.loglog(self.frequency, self.power_mat[:,20,20])
        #plt.show()
        # Next steps: Take mean over frequency and then factor the chan-chan
        # matrix.  Pick out the top n eigenvalues, then rotate the unaveraged
        # power mat to that basis.  Look at the spectrum of thoses modes and
        # the spectrums of leftover parts.  Save thermal parameters, the modes
        # and the parameters for the spectra of those modes.

def ps_freq_axis(dt, n):
    """Calculate the frequency axis for a power spectrum.

    Parameters
    ----------
    dt : float
        time step for time stream that went into power spectrum calculation.
    n : int
        Number of time stream data points that went into the power spectrum
        calculation.

    Returns
    -------
    frequencies : 1D array of floats.
        The frequency axis.  Has length `n`//2.
    """

    frequency = sp.arange(n//2, dtype=float)
    df = 1.0/dt/n
    frequency *= df
    return frequency

def make_masked_time_stream(Blocks, ntime=None, window=None, 
                            return_means=False, subtract_slope=False) :
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
        to use exactly the number that fits all the data.  Set to a negitive
        factor to zero pad to a power of 2 and by at least at least the factor.
    window : string or tuple
        Type of window to apply to each DataBlock.  Valid options are the valid
        arguments to scipy.signal.get_window().  By default, don't window.
    return_means : bool
        Whether to return an array of the channed means.
    subtract_slope : bool
        Whether to subtract a linear function of time from each channel.

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
    means : array (optional)
        The mean from each channel.
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
    mean_time = 0.0
    n_data_times = 0
    for Data in Blocks :
        Data.calc_time()
        min_time = min(min_time, min(Data.time))
        max_time = max(min_time, max(Data.time))
        mean_time += sp.sum(Data.time)
        n_data_times += len(Data.time)
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
    mean_time /= n_data_times
    # Calculate the time axis.
    if not ntime :
        ntime = (max_time - min_time) // dt + 1
    elif ntime < 0:
        # 0 pad by a factor of at least -ntime, but at most 10% more than this.
        time_min = -ntime * (max_time - min_time) / dt
        n_block = 1
        while n_block < time_min/20.0:
            n_block *= 2
        ntime = (time_min//n_block  + 1) * n_block

    time = sp.arange(ntime)*dt + min_time
    # Allowcate memory for the outputs.
    time_stream = sp.zeros((ntime,) + back_shape, dtype=float)
    mask = sp.zeros((ntime,) + back_shape, dtype=sp.float32)
    # Very important to subtract the mean out of the signal, otherwise the
    # window coupling to the mean (0) mode will dominate everything.
    total_sum = 0.0
    total_counts = 0
    total_slope = 0.0
    time_norm = 0.0
    for Data in Blocks:
        total_sum += sp.sum(Data.data.filled(0), 0)
        total_counts += ma.count(Data.data, 0)
        total_slope += sp.sum(Data.data.filled(0) 
                              * (Data.time[:,None,None,None] - mean_time), 0)
        time_norm += sp.sum(sp.logical_not(ma.getmaskarray(Data.data))
                            * (Data.time[:,None,None,None] - mean_time)**2, 0)
    total_counts[total_counts == 0] = 1
    time_norm[time_norm == 0.0] = 1
    total_mean = total_sum / total_counts
    total_slope /= time_norm
    # Loop over all times and fill in the arrays.
    for Data in Blocks :
        # Subtract the mean calculated above.
        Data.data -= total_mean
        # If desired, subtract of the linear function of time.
        if subtract_slope:
            Data.data -= (total_slope 
                          * (Data.time[:,None,None,None] - mean_time))
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
    if return_means:
        return time_stream, mask, dt, total_mean
    else :
        return time_stream, mask, dt

def windowed_power(data1, window1, data2=None, window2=None, axis=-1) :
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

    Returns
    -------
    power_spectrum : 1D array
        Cross power of the input time streams accounting for the windows.  It
        will be the same shape as the input arrays accept for the `axis` axis
        will be shortened to n//2.
        If input data has units Kelvin, power spectrum has units Kelvin**2
        (not Kelvin**2 / bin).
    
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
    
    data_power = calculate_power(data1, data2, axis)
    window_power = calculate_power(window1, window2, axis)
    true_power = deconvolve_power(data_power, window_power, axis=axis)
    return true_power

def calculate_power(data1, data2=None, axis=-1):
    """Calculate the cross or auto power spectum.
    
    If input data has units Kelvin, power spectrum has units Kelvin**2 (not
    Kelvin**2 / bin).
    """
    
    if data2 is None :
        data2 = data1
    power = fft.fft(data1, axis=axis)*fft.fft(data2, axis=axis).conj()
    power = power / data1.shape[axis]
    return power

def prune_power(power_spectrum, axis=-1):
    """Discard negitive frequencies from power spectrum.

    Parameters
    ----------
    power_spectrum : array
        A power spectrum to be pruned of negitive frequencies.
    axis : int
        The axis representing the frequency direction.

    Returns
    -------
    pruned_power_spectrum : array
        The modified power spectrum.  Similar shape as `power_spectrum` except
        for the `axis` axis, which is reduced to n//2.
    """

    n = power_spectrum.shape[axis]
    # Build up the slice for this axis.
    s = slice(n//2)
    indices = [slice(sys.maxsize)] * power_spectrum.ndim
    indices[axis] = s
    indices = tuple(indices)
    # Truncate the power spectrum to only the unaliased and positive
    # frequencies.
    pruned_power_spectrum = power_spectrum[indices]
    return pruned_power_spectrum

def make_power_physical_units(power, dt):
    """Converts power spectrum to physical units given.

    Given the time stream time step, devide by the bandwidth to put the power
    spectrum in the most physical units of amp**2/Hz.
    """
    
    return power*dt*2

def deconvolve_power(power_spectrum, window_power, axis=-1):
    """Deconvolve a power spectrum with a window function."""
    
    data_corr = fft.ifft(power_spectrum, axis=axis)
    window_corr = fft.ifft(window_power, axis=axis)
    true_corr = data_corr / window_corr
    true_power = fft.fft(true_corr, axis=axis, overwrite_x=True)
    return true_power

def convolve_power(power_spectrum, window_power, axis=-1):
    """Convolve a power spectrum with a window function."""
    
    data_corr = fft.ifft(power_spectrum, axis=axis)
    window_corr = fft.ifft(window_power, axis=axis)
    true_corr = data_corr * window_corr
    true_power = fft.fft(true_corr, axis=axis, overwrite_x=True)
    return true_power

def calculate_overf_correlation(amp, index, f0, dt, n_lags):
    """Calculates the correation function for a 1/f power spectrum.
    """
    
    # Cast inputs as floats as I do a bunch of division.
    dt = float(dt)
    f0 = float(f0)
    index = float(index)
    # Number of points used in calculation needs to be at least 10 times bigger
    # than final number of point returned.  Need the frequency resolution for
    # numerical accuracy.
    n = 20*n_lags
    n_return = n_lags
    # Generate the powerspectrum.
    power = overf_power_spectrum(amp, index, f0, dt, n)
    # FFT it to the correlation function.
    corr = fft.ifft(power)
    corr = corr[:n_return].real
    # Bump the whole thing up to the last lag is at 0 correlation (as opposed
    # to being negitive).
    corr -= corr[-1]
    # To normalize, need to multiply by the bandwidth.
    corr *= 1.0/2/dt
    return corr

def overf_power_spectrum(amp, index, f0, dt, n):
    """Calculates the theoredical f**`index` power spectrum.
    """
    
    # Get the frequencies represented in the FFT.
    df = 1.0/dt/n
    freq = sp.arange(n, dtype=float)
    freq[n//2+1:] -= freq[-1] + 1
    freq = abs(freq)*df
    # 0th mode is meaningless.  Set to unity to avoid errors.
    freq[0] = 1
    # Make the power spectrum.
    power = amp*(freq/f0)**index
    # Set the mean mode explicitly to 0.
    power[0] = 0
    return power

def generate_overf_noise(amp, index, f0, dt, n):
    """Generate noise time series with a f**`index` power spectrum."""

    white_noise = rand.normal(size=n)
    power_spectrum = overf_power_spectrum(amp, index, f0, dt, n)
    # Power spectrum is in physical units of T**2/Hz.  Put in discrete units by
    # multiplying by the bandwidth.
    power_spectrum *= 1.0/2.0/dt
    noise = fft.ifft(fft.fft(white_noise)*sp.sqrt(power_spectrum)).real
    return noise

def full_power_mat(Blocks, n_time=None, window=None, deconvolve=True,
                   subtract_slope=False) :
    """Calculate the full power spectrum of a data set with channel
    correlations.
    
    Only one cal state and pol state assumed.
    """
    
    full_data, mask, dt, channel_means = make_masked_time_stream(Blocks, n_time, 
        window=window, return_means=True, subtract_slope=subtract_slope)
    n_time = full_data.shape[0]
    n_pol = full_data.shape[1]
    n_cal = full_data.shape[2]
    # Broadcast to shape such that all pairs of channels are calculated.
    full_data1 = full_data[:,:,:,:,None]
    full_data2 = full_data[:,:,:,None,:]
    mask1 = mask[:,:,:,:,None]
    mask2 = mask[:,:,:,None,:]
    n_chan = full_data.shape[-1]
    # Loop to do only a bit of the calculation at a time.  Reduces
    # memory use by a factor of a few.
    power_mat = sp.empty((n_time, n_pol, n_cal, n_chan, n_chan),
                              dtype=complex)
    window_function = sp.empty((n_time, n_pol, n_cal, n_chan, n_chan),
                              dtype=complex)
    # Break the calculation up over a loop to conserve memory.
    for ii in xrange(n_chan):
        w = calculate_power(mask1[:,:,:,[ii],:], mask2, axis=0)
        # Calculate the window normalizations.
        window_norms = sp.mean(w, 0)
        # Do nothing to fully masked channels to keep things finite.
        window_norms[window_norms < 1e-7] = 1
        if deconvolve:
            p = windowed_power(full_data1[:,:,:,[ii],:],
                               mask1[:,:,:,[ii],:], full_data2, mask2, axis=0)
            power_mat[:,:,:,[ii],:] = p
        else:
            p = calculate_power(full_data1[:,:,:,[ii],:], full_data2, axis=0)
            # Normalize from windowing.
            power_mat[:,:,:,[ii],:] = p / window_norms
        w /= window_norms
        window_function[:,:,:,[ii],:] = w
    return power_mat, window_function, dt, channel_means


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    NoisePower(str(sys.argv[1])).execute()
