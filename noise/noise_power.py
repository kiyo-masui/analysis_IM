"""Proceedure to estimate a GBT noise model from data."""

import os
import sys
import math

import scipy as sp
import scipy.fftpack as fft
import scipy.signal as sig
import scipy.linalg as linalg
import scipy.special
import numpy.random as rand
import numpy.ma as ma
# XXX for testing, but needs to be commented for production.
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
import core.fitsGBT
from utils import misc


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
    
    def get_data(self, file_middle):
        params = self.params
        cal_weights = params['cal_weights']
        pol_weights = params['pol_weights']
        n_time = self.n_time
        window = params['window']
        subtract_slope = params['subtract_slope']
        input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
        # Read in the data.
        Reader = core.fitsGBT.Reader(input_fname)
        Blocks = Reader.read(params['scans'], params['IFs'],
                             force_tuple=True)
        # On the first pass, set the channel width.
        if not hasattr(self, "chan_width"):
            self.chan_width = Blocks[0].field['CDELT1']
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
        # Convert the data to the proper format and return it.
        return make_masked_time_stream(Blocks, n_time, window=window, 
                                       return_means=True, 
                                       subtract_slope=subtract_slope)
    
    def execute(self):
        self.power_mat, self.thermal_expectation = self.full_calculation()
        n_chan = self.power_mat.shape[1]
        n_freq = self.power_mat.shape[0]
        # Calculate the the mean channel correlations at low frequencies.
        low_f_mat = sp.mean(self.power_mat[1:4 * n_chan + 1,:,:], 0).real
        # Factorize it into preinciple components.
        e, v = linalg.eigh(low_f_mat)
        self.low_f_mode_values = e
        # Make sure the eigenvalues are sorted.
        if sp.any(sp.diff(e) < 0):
            raise RuntimeError("Eigenvalues not sorted.")
        self.low_f_modes = v
        # Now subtract out the noisiest channel modes and see what is left.
        n_modes_subtract = 10
        mode_subtracted_power_mat = sp.copy(self.power_mat.real)
        mode_subtracted_auto_power = sp.empty((n_modes_subtract, n_freq))
        for ii in range(n_modes_subtract):
            mode = v[:,-ii]
            amp = sp.sum(mode[:,None] * mode_subtracted_power_mat, 1)
            amp = sp.sum(amp * mode, 1)
            to_subtract = amp[:,None,None] * mode[:,None] * mode
            mode_subtracted_power_mat -= to_subtract
            auto_power = mode_subtracted_power_mat.view()
            auto_power.shape = (n_freq, n_chan**2)
            auto_power = auto_power[:,::n_chan + 1]
            mode_subtracted_auto_power[ii,:] = sp.mean(auto_power, -1)
        self.subtracted_auto_power = mode_subtracted_auto_power

    def full_calculation(self, chan_modes_subtract=None, n_poly_subtract=0):
        
        # Some set up.
        params = self.params
        deconvolve = params['deconvolve']
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        self.n_time = params["n_time_bins"]
        n_files = len(params["file_middles"])
        # Loop over files to process.
        first_iteration = True
        for file_middle in params['file_middles'] :
            # Get the data.
            full_data, mask, this_dt, full_mean = self.get_data(file_middle)
            if first_iteration :
                self.n_time = full_data.shape[0]
                n_chan = full_data.shape[-1]
                dt = this_dt
            elif not sp.allclose(dt, this_dt, rtol=0.001):
                msg = "Files have different time samplings."
                raise ce.DataError(msg)
            # Subtract out any channel modes passed in.
            if not (chan_modes_subtract is None):
                for v in chan_modes_subtract:
                    full_data -= v * sp.sum(v * full_data, -1)[:, None]
            # Subtract out polynomials from each channel if desired.
            if first_iteration:
                # Generate basis polynomials.
                basis_poly = sp.empty((n_poly_subtract, self.n_time))
                time_scaled = ((sp.arange(self.n_time, dtype=float) * 2
                                - self.n_time + 1.0) / self.n_time)
                for ii in range(n_poly_subtract):
                    #tmp_poly = scipy.special.eval_chebyu(ii, time_scaled)
                    tmp_poly = sp.cos(sp.pi*ii*time_scaled)
                    tmp_poly *= 1.0/sp.sqrt(sp.sum(tmp_poly**2))
                    basis_poly[ii,:] = tmp_poly
                # Allocate memory to hold the amplitude spectrum.
                poly_spectrum = sp.zeros((n_poly_subtract, n_chan), 
                                         dtype=float)
            # Fit for the polynomials in each channel.
            for ii in range(n_chan):
                weighted_poly = basis_poly * mask[:,0,0,ii]
                poly_corr = sp.sum(full_data[:,0,0,ii]
                                   * basis_poly[:,:], -1)
                poly_covar = sp.sum(weighted_poly[:,None,:] 
                                    * basis_poly[None,:,:], -1)
                if n_poly_subtract:
                    poly_amps = linalg.solve(poly_covar, poly_corr, 
                                             sym_pos=True, overwrite_a=True,
                                             overwrite_b=True)
                    poly_spectrum[:,ii] += poly_amps**2
            # Calculate the raw power spectrum.
            power_mat, window_function = calculate_full_power_mat(full_data, 
                                                mask, deconvolve=deconvolve)

            # Get rid of the extra cal and polarization axes.
            power_mat = power_mat[:,0,0,:,:]
            window_function = window_function[:,0,0,:,:]
            full_mean = full_mean[0,0,:]
            # TODO: Figure out a better way to deal with this (only drop
            # affected frequencies).
            #if sp.any(sp.allclose(mask[:,0,0,:], 0.0, 0)):
            #    n_files -= 1
            #    continue
            # Format the power spectrum.
            power_mat = prune_power(power_mat, 0)
            power_mat = make_power_physical_units(power_mat, this_dt)
            # TODO In the future the thermal expectation could include
            # polarization factors (be 'I' aware) and cal factors.
            thermal_expectation = (full_mean / sp.sqrt(this_dt)  /
                                   sp.sqrt(abs(self.chan_width)) /
                                   sp.sqrt(1.0 / 2.0 / this_dt))
            if params['norm_to_thermal'] :
                power_mat /= (thermal_expectation[:,None] 
                              * thermal_expectation[None,:])
            # Combine across files.
            if first_iteration :
                total_power_mat = power_mat
                total_thermal_expectation = thermal_expectation
            else :
                total_power_mat += power_mat
                total_thermal_expectation += thermal_expectation
            first_iteration = False
        if not hasattr(self, 'frequency'):
            self.frequency = ps_freq_axis(dt, self.n_time)
        if n_files > 0 :
            total_power_mat /= n_files
            total_thermal_expectation / n_files
            poly_spectrum /= n_files
        if n_poly_subtract:
            return total_power_mat, total_thermal_expectation, poly_spectrum
        else:
            return total_power_mat, total_thermal_expectation
        

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
    # Also get the time axis and the mask
    # for calculating basis polynomials.
    unmask = sp.zeros((0,) + back_shape, dtype=bool)
    time = sp.zeros((0,), dtype=float)
    start_ind = []
    min_time = float('inf')
    max_time = 0.0
    #mean_time = 0.0
    #n_data_times = 0
    for Data in Blocks :
        Data.calc_time()
        start_ind.append(len(time))
        time = sp.concatenate((time, Data.time))
        this_unmask = sp.logical_not(ma.getmaskarray(Data.data))
        unmask = sp.concatenate((unmask, this_unmask), 0)
        # Often the start or the end of a scan is completly masked.  Make sure
        # we don't start till the first unmasked time and end at the last
        # unmasked time.
        time_unmask = sp.alltrue(ma.getmaskarray(Data.data), -1)
        time_unmask = sp.alltrue(time_unmask, -1)
        time_unmask = sp.alltrue(time_unmask, -1)
        if sp.alltrue(time_unmask):
            continue
        time_unmask = sp.logical_not(time_unmask)
        min_time = min(min_time, min(Data.time[time_unmask]))
        max_time = max(min_time, max(Data.time[time_unmask]))
        #mean_time += sp.sum(Data.time[time_unmask])
        #n_data_times += len(Data.time[time_unmask])
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
    # Now calculate basis polynomials for the mean mode and the slope mode.
    polys = misc.ortho_poly(time[:,None,None,None], 2, unmask, 0)
    #mean_time /= n_data_times
    #if n_data_times == 0:
    #    n_data_times = 1
    # Very important to subtract the mean out of the signal, otherwise the
    # window coupling to the mean (0) mode will dominate everything. Can also
    # optionally take out a slope.
    # Old algorithm.
    #total_sum = 0.0
    #total_counts = 0
    #total_slope = 0.0
    #time_norm = 0.0
    #for Data in Blocks:
    #    total_sum += sp.sum(Data.data.filled(0), 0)
    #    total_counts += ma.count(Data.data, 0)
    #    total_slope += sp.sum(Data.data.filled(0) 
    #                          * (Data.time[:,None,None,None] - mean_time), 0)
    #    time_norm += sp.sum(sp.logical_not(ma.getmaskarray(Data.data))
    #                        * (Data.time[:,None,None,None] - mean_time)**2, 0)
    #total_counts[total_counts == 0] = 1
    #time_norm[time_norm == 0.0] = 1
    #total_mean = total_sum / total_counts
    #total_slope /= time_norm
    # New algorithm.
    mean_amp = 0
    slope_amp = 0
    for ii, Data in enumerate(Blocks):
        si = start_ind[ii]
        this_nt = Data.dims[0]
        data = Data.data.filled(0)
        mean_amp += sp.sum(data * unmask[si:si + this_nt,...]
                           * polys[0,si:si + this_nt,...], 0)
        slope_amp += sp.sum(data * unmask[si:si + this_nt,...]
                            * polys[1,si:si + this_nt,...], 0)
    polys[0,...] *= mean_amp
    polys[1,...] *= slope_amp
    # Calculate the time axis.
    if min_time > max_time:
        min_time = 0
        max_time = 6 * dt
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
    # Loop over all times and fill in the arrays.
    for ii, Data in enumerate(Blocks):
        this_nt = Data.dims[0]
        si = start_ind[ii]
        # Subtract the mean calculated above.
        this_data = Data.data.copy()
        this_data -= polys[0,si:si + this_nt,...]
        # If desired, subtract of the linear function of time.
        if subtract_slope:
            #this_data -= (total_slope 
            #              * (Data.time[:,None,None,None] - mean_time))
            this_data -= polys[1,si:si + this_nt,...]
        # Find the first and last unmasked times.
        time_unmask = sp.alltrue(ma.getmaskarray(this_data), -1)
        time_unmask = sp.alltrue(time_unmask, -1)
        time_unmask = sp.alltrue(time_unmask, -1)
        if sp.alltrue(time_unmask):
            continue
        time_unmask = sp.logical_not(time_unmask)
        unmasked_ind, = sp.where(time_unmask)
        first_ind = min(unmasked_ind)
        last_ind = max(unmasked_ind)
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
        # Apply an offset to the time in case the start of the Data Block
        # doesn't line up with the time array perfectly.
        offset = (time[sp.argmin(abs(time - Data.time[first_ind]))]
                  - Data.time[first_ind])
        # Generate window function.
        if window:
            window_function = sig.get_window(window, last_ind - first_ind + 1)
        for ii in range(first_ind, last_ind + 1) :
            ind = sp.argmin(abs(time - (Data.time[ii] + offset)))
            if abs(time[ind] - (Data.time[ii])) < 0.5*dt :
                if sp.any(mask[ind, ...]) :
                    msg = "Overlapping times in Data Blocks."
                    raise ce.DataError(msg)
                if window:
                    window_value = window_function[ii - first_ind]
                else :
                    window_value = 1.0
                time_stream[ind, ...] = (window_value 
                                         * this_data[ii, ...].filled(0.0))
                mask[ind, ...] = window_value * sp.logical_not(ma.getmaskarray(
                                     this_data)[ii, ...])
    if return_means:
        return time_stream, mask, dt, polys[0,0,...]
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

    Given the time stream time step, devide by the twice the bandwidth to 
    put the power spectrum in the most physical units of amp**2/Hz.
    """
    
    return power*dt

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
    """Calculates the correlation function for a 1/f power spectrum.
    """
    
    # Cast inputs as floats as I do a bunch of division.
    dt = float(dt)
    f0 = float(f0)
    index = float(index)
    # Number of points used in calculation needs to be at least 10 times bigger
    # than final number of point returned.  This is so we are not affected by
    # the periodicity of the correlation function.
    buff_factor = 64
    n = buff_factor * n_lags
    n_return = n_lags
    # Generate the power spectrum.
    # Need to add a low frequency cut off, since there is an IR divergence.
    # Choose to cut off at 1/2df (so we get a bit of slope mode).
    power = overf_power_spectrum(amp, index, f0, dt, n,
                                 cut_off=1./n_lags/dt/2.0)
    # FFT it to the correlation function.
    corr = fft.ifft(power)
    # Complex part should be zero.
    corr = corr.real
    # In previous versions of this function, we shifted the output function.
    # however this screws up positive definiteness of the correlation matrix
    # and is unnecessary if you have the IR cut off.
    #corr -= corr[2 * n_return]
    # Trim to return size.
    corr = corr[:n_return]
    # To normalize, need to multiply by twice the bandwidth.
    corr *= 1.0/dt
    return corr

def overf_power_spectrum(amp, index, f0, dt, n, cut_off=0):
    """Calculates the theoretical f**`index` power spectrum.
    """
    
    if cut_off < 0:
        raise ValueError("Low frequency cut off must not be negative.")
    # Sometimes the fitting routines do something weird that causes
    # an overflow from a ridiculous index.  Limit the index.
    index = max(index, -20)
    # Get the frequencies represented in the FFT.
    df = 1.0/dt/n
    freq = sp.arange(n, dtype=float)
    freq[n//2+1:] -= freq[-1] + 1
    freq = abs(freq)*df
    # 0th (mean) mode is meaningless has IR divergence.  Deal with it later (in
    # the cut off.
    freq[0] = 1
    # Make the power spectrum.
    power = (freq/f0)**index
    power *= amp
    # Restore frequency of mean mode.
    freq[0] = 0
    # Find the power just above the cut off frequency.
    p_cut = power[sp.amin(sp.where(freq > cut_off)[0])]
    # Flatten off the power spectrum.
    power[freq <= cut_off] = p_cut
    return power

def generate_overf_noise(amp, index, f0, dt, n):
    """Generate noise time series with a f**`index` power spectrum."""

    white_noise = rand.normal(size=n)
    power_spectrum = overf_power_spectrum(amp, index, f0, dt, n)
    # Power spectrum is in physical units of T**2/Hz.  Put in discrete units by
    # multiplying by twice the bandwidth.
    power_spectrum *= 1.0/dt
    noise = fft.ifft(fft.fft(white_noise)*sp.sqrt(power_spectrum)).real
    return noise

def calculate_full_power_mat(full_data, mask, deconvolve=True, normalize=True):
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
        # Calculate the window normalizations.  This is critical as several
        # functions that call this one assume this normalization.
        window_norms = sp.sum(w.real, 0) / float(n_time)
        # Do nothing to fully masked channels to keep things finite.
        window_norms[window_norms < 1e-10] = 1
        if deconvolve:
            p = windowed_power(full_data1[:,:,:,[ii],:],
                               mask1[:,:,:,[ii],:], full_data2, mask2, axis=0)
            power_mat[:,:,:,[ii],:] = p
        else:
            p = calculate_power(full_data1[:,:,:,[ii],:], full_data2, axis=0)
            # Normalize from windowing.
            if normalize:
                p /= window_norms
            power_mat[:,:,:,[ii],:] = p
        if normalize and not deconvolve:
            w /= window_norms
        window_function[:,:,:,[ii],:] = w
    return power_mat, window_function

def full_power_mat(Blocks, n_time=None, window=None, deconvolve=True,
                   subtract_slope=False, normalize=True, split_scans=False) :
    """Calculate the full power spectrum of a data set with channel
    correlations.
    
    Only one cal state and pol state assumed... Don't think this is true -KM.
    """
    
    if split_scans:
        if n_time < 0:
            nt = 0
            for Data in Blocks:
                nt = max(nt, Data.dims[0])
            time_min = -n_time * nt
            n_block = 1
            while n_block < time_min/20.0:
                n_block *= 2
            n_time = (time_min//n_block  + 1) * n_block
        back_dims = Blocks[0].dims[1:]
        n_chan = back_dims[-1]
        power_mat = sp.zeros((n_time,) + back_dims + (n_chan,), dtype=complex)
        window_function = sp.zeros((n_time,) + back_dims + (n_chan,),
                                   dtype=complex)
        channel_means = sp.zeros(back_dims, dtype=float)
        for ii, Data in enumerate(Blocks):
            this_data, this_mask, this_dt, this_means = \
                    make_masked_time_stream((Data,), n_time, window=window,
                                            return_means=True, 
                                            subtract_slope=subtract_slope)
            channel_means += this_means
            if ii == 0:
                dt = this_dt
            else:
                if not sp.allclose(dt, this_dt, rtol=0.001):
                    raise RuntimeError("Time sampling doesn't line up.")
            this_power, this_window = calculate_full_power_mat(this_data, 
                        this_mask, deconvolve=deconvolve, normalize=normalize)
            power_mat += this_power / len(Blocks)
            window_function += this_window / len(Blocks)
        channel_means /= len(Blocks)
    else:
        full_data, mask, dt, channel_means = make_masked_time_stream(Blocks, 
            n_time, window=window, return_means=True,
            subtract_slope=subtract_slope)
        power_mat, window_function = calculate_full_power_mat(full_data, mask, 
                        deconvolve=deconvolve, normalize=normalize)
    return power_mat, window_function, dt, channel_means


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    NoisePower(str(sys.argv[1])).execute()
