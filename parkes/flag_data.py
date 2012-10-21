#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import os
import copy

import numpy.ma as ma
import scipy as sp
import scipy.signal as sig

import fitsGBT
import base_single

import kiyopy.custom_exceptions as ce
from time_stream import hanning
from time_stream import cal_scale
from time_stream import rotate_pol

# XXX
#import matplotlib.pyplot as plt

class FlagData(base_single.BaseSingle) :
    '''Pipeline module that flags RFI and other forms of bad data.

    '''

    prefix = 'fd_'
    params_init = {
                   # In multiples of the standard deviation of the whole block
                   # once normalized to the time median.
                   'perform_hanning' : False,
                   'cal_scale' : False,
                   # Rotate to XX,XY,YX,YY is True.
                   'rotate' : False,
                   # Any frequency with variance > sigma_thres sigmas will be 
                   # flagged (recursively).
                   'sigma_thres' : 6.,
                   # A Data that has more than badness_thres frequencies flagged
                   # (as a fraction) will be considered bad.
                   'badness_thres' : 0.1,
                   # How many times to hide around a bad time.
                   'time_cut' : 40
                   }
    feedback_title = 'New flags each data Block: '
    
    def action(self, Data):
        '''Prepares Data and flags RFI.
        
        Parameters
        ----------
        Data : DataBlock
            Contains information in a usable format direct from GBT. 

        Returns
        -------
        Data : DataBlock
            The input `Data` with RFI flagged. Will also be cal scaled and
            rotated to XX,YY... if so chosen.

        '''
        params = self.params
        # Keep track of how many pre existing flags there are for feedback
        # purposes.
        already_flagged = ma.count_masked(Data.data)
        if params["rotate"]:
            if (tuple(Data.field['CRVAL4']) == (1, 2, 3, 4)):
                rotate_pol.rotate(Data, (-5,-7,-8,-6))
                Data.add_history('Rotated to XX,XY,YX,YY')
        # Few operations to be performed before flagging.
        if params["perform_hanning"] :
            hanning.hanning_smooth(Data)
            Data.add_history('Hanning smoothed.')
        if params["cal_scale"] :
            cal_scale.scale_by_cal(Data, True, False, False)
            Data.add_history('Converted to units of noise cal temperture.')
        # Flag the data.
        apply_cuts(Data, sigma_thres=params['sigma_thres'], 
                    badness_thres=params['badness_thres'],
                    time_cut=params['time_cut'])
        Data.add_history('Flagged Bad Data.', ('Sigma threshold: '
                    + str(self.params['sigma_thres']), 'Badness threshold: '
                    + str(self.params['badness_thres']), 'Time mask size: '
                    + str(self.params['time_cut'])))
        # Report the number of new flags.
        new_flags = ma.count_masked(Data.data) - already_flagged
        percent = float(new_flags) / Data.data.size * 100
        self.block_feedback = '%d (%f%%), ' % (new_flags, percent)

        # time tsys, make data in unit of Jy
        Data.data -= Data.data.mean(axis=0)
        Data.data *= Data.field['TSYS'][:,:,:,:,None]
        return Data

def apply_cuts(Data, sigma_thres=6, badness_thres=0.1, time_cut=40):
    '''Flags bad data from RFI and far outliers.

    See `flag_data()` for parameter explanations and more info.
    '''
    badness = flag_data(Data, sigma_thres, badness_thres, time_cut)
    # Can print or return badness here if you would like
    # to see if the Data had a problem in time or not.
    return

def flag_data(Data, sigma_thres, badness_thres, time_cut):
    '''Flag bad data from RFI and far outliers.

    Parameters
    ----------
    Data : DataBlock
        Contains information in a usable format direct from GBT. Bad
        frequencies will be flagged in all polarizations and cal states.
    sigma_thres : int or float
        Any frequency with variance > `sigma_thres` sigmas will be 
        flagged (recursively).
    badness_thres : float
        A `Data` that has more than `badness_thres` frequencies flagged
        (as a fraction) will be considered 'bad'. `0` means that everything
        will be considered bad while `1` means nothing will be. 
    time_cut : int
        How many time bins (as an absolute number) to flag if `Data` has been
        considered 'bad'. See `destroy_time_with_mean_arrays` for more 
        infomation on this.

    Returns
    -------
    badness : bool
        Returns `True` iff a `Data` has been considered 'bad'.
    
    Notes
    -----
    'badness' is when more than a certain fraction of freqs has been flagged
    from `Data`. This certain fraction comes from `badness_thres`. `Data` that
    is 'bad' has a lot of frequencies flagged and this can because a lot of 
    frequencies are actually bad or because there was a blip in time (maybe
    the machine choked for a second).
    If a `Data` gets considered 'bad' then the algorithm tries to find
    something wrong in time (and masks those bad times) and redoes the RFI
    flagging. If there is a significant decrease (5%) in the number of 
    frequencies flagged, then the problem was in time and it uses the mask
    from this second run with bad times flagged. If not, then the `Data` is
    bad either way and it uses the mask from the first run. Increasing the
    `time_cut` in this situation is not recommended since you lose a lot more
    data (there are 10 times as many freq. bins as time bins in `Data`). 
    '''
    # Flag data on a [deep]copy of Data. If too much destroyed,
    # check if localized in time. If that sucks too, then just hide freq.

    Data1 = copy.deepcopy(Data)
    itr = 0            # For recursion
    max_itr = 20       # For recursion
    bad_freqs = []
    amount_masked = -1 # For recursion
    while not (amount_masked == 0) and itr < max_itr:
        amount_masked = destroy_with_variance(Data1, sigma_thres, bad_freqs) 
        itr += 1
    bad_freqs.sort()
    # Remember the flagged data.
    mask = Data1.data.mask
    # Check for badness.
    percent_masked1 = (float(len(bad_freqs)) / Data1.dims[-1])
    badness = (percent_masked1 > badness_thres)
    # If too many frequencies flagged, it may be that the problem
    # happens in time, not in frequency.
    if badness:
        Data2 = copy.deepcopy(Data)
        # Mask the bad times.
        destroy_time_with_mean_arrays(Data2, flag_size=time_cut)
        # Then try to flag again with bad times masked.
        # Bad style for repeating as above, sorry.
        itr = 0
        bad_freqs = []
        amount_masked = -1
        while not (amount_masked == 0) and itr < max_itr:
            amount_masked = destroy_with_variance(Data2, bad_freq_list=bad_freqs) 
            itr += 1
        bad_freqs.sort()
        percent_masked2 = (float(len(bad_freqs)) / Data2.dims[-1])
        # If the data is 5% or more cleaner this way <=> it is not bad.
        badness = (percent_masked1 - percent_masked2) < 0.05
        # If this data does not have badness, that means there was
        # a problem in time and it was solved, so use this mask.
        # If the data is still bad, then the mask from Data1 will be used.
        if not badness:
            itr = 0
            while not (amount_masked == 0) and itr < max_itr:
                amount_masked = destroy_with_variance(Data2, sigma_thres,
                                                      bad_freqs) 
                itr += 1
            Data1 = Data2
    # We've flagged the RFI down to the foreground limit.  Filter out the
    # foregrounds and flag again to get below the foreground limit.
    # TODO, hard coded time_bins_smooth acctually depends on the scan speed and
    # the time sampling.
    filter_foregrounds(Data1, n_bands=40, time_bins_smooth=10)
    itr = 0 
    while not (amount_masked == 0) and itr < max_itr:
        amount_masked = destroy_with_variance(Data1, sigma_thres, bad_freqs) 
        itr += 1
    mask = Data1.data.mask
    # Finally copy the mask to origional data block.
    Data.data.mask = mask
    return badness

def destroy_with_variance(Data, sigma_thres=6, bad_freq_list=[]):
    '''Mask frequencies with high variance.

    Since the signal we are looking for is much weaker than what is in `Data`,
    any frequency that is 'too spiky' is not signal and is RFI instead. Using
    variance as a test really makes this 'spikyness' stand out.

    Parameters
    ----------
    Data : DataBlock
        Contains information in a usable format direct from GBT. Bad
        frequencies will be flagged in all polarizations and cal states.
    sigma_thres : int or float
        Any frequency with variance > `sigma_thres` sigmas will be 
        flagged (recursively).
    bad_freq_list : list of int
        A list of bad frequencies. Since this method is called over and over,
        this list keeps track of what has been flagged. Bad frequencies that
        are found will be appended to this list.

    Returns
    -------
    amount_masked : int
        The amount of frequencies masked.
    amount_masked_max : int
        The lagest amount of frequences masked in different beams

    Notes
    -----
    Data must have 5 dim: [time beam pol cal freq]
    Polarizations must be in XX,YY format.
    Calibration must be off
    Operate on each beam seperated

    '''
    shape = Data.data.shape
    amount_masked_max = 0
    # Loop over each beam
    for ii in range(shape[1]):
        # Get the normalized variance array for each polarization.
        a = ma.var(Data.data[:,ii,0,0,:],0)/(ma.mean(Data.data[:,ii,0,0,:],0)**2)#XX
        b = ma.var(Data.data[:,ii,1,0,:],0)/(ma.mean(Data.data[:,ii,1,0,:],0)**2)#YY
        # Get the mean and standard deviation [sigma].
        means = sp.array([ma.mean(a), ma.mean(b)]) 
        sig   = sp.array([ma.std(a), ma.std(b)])
        # Get the max accepted value [sigma_thres*sigma, sigma_thres=6 works really well].
        max_sig = sigma_thres*sig
        max_accepted = means + max_sig
        amount_masked = 0
        for freq in range(0, len(a)):
            if ((a[freq] > max_accepted[0]) or (b[freq] > max_accepted[1])):
                # mask
                amount_masked += 1
                bad_freq_list.append(freq)
                Data.data[:,:,:,:,freq].mask = True
        if amount_masked > amount_masked_max:
            amount_masked_max = amount_masked 
    return amount_masked_max

def destroy_with_variance_4pol2cal(Data, sigma_thres=6, bad_freq_list=[]):
    '''Mask frequencies with high variance.

    Since the signal we are looking for is much weaker than what is in `Data`,
    any frequency that is 'too spiky' is not signal and is RFI instead. Using
    variance as a test really makes this 'spikyness' stand out.

    Parameters
    ----------
    Data : DataBlock
        Contains information in a usable format direct from GBT. Bad
        frequencies will be flagged in all polarizations and cal states.
    sigma_thres : int or float
        Any frequency with variance > `sigma_thres` sigmas will be 
        flagged (recursively).
    bad_freq_list : list of int
        A list of bad frequencies. Since this method is called over and over,
        this list keeps track of what has been flagged. Bad frequencies that
        are found will be appended to this list.

    Returns
    -------
    amount_masked : int
        The amount of frequencies masked.

    Notes
    -----
    Polarizations must be in XX,XY,YX,YY format.

    '''
    XX_YY_0 = ma.mean(Data.data[:, 0, 0, :], 0) * ma.mean(Data.data[:, 3, 0, :], 0)
    XX_YY_1 = ma.mean(Data.data[:, 0, 1, :], 0) * ma.mean(Data.data[:, 3, 1, :], 0)
    # Get the normalized variance array for each polarization.
    a = ma.var(Data.data[:, 0, 0, :], 0) / (ma.mean(Data.data[:, 0, 0, :], 0)**2) # XX
    b = ma.var(Data.data[:, 1, 0, :], 0) / XX_YY_0                                # XY
    c = ma.var(Data.data[:, 2, 0, :], 0) / XX_YY_0                                # YX
    d = ma.var(Data.data[:, 3, 0, :], 0) / (ma.mean(Data.data[:, 3, 0, :], 0)**2) # YY
    # And for cal off.
    e = ma.var(Data.data[:, 0, 1, :], 0) / (ma.mean(Data.data[:, 0, 1, :], 0)**2) # XX
    f = ma.var(Data.data[:, 1, 1, :], 0) / XX_YY_1                                # XY
    g = ma.var(Data.data[:, 2, 1, :], 0) / XX_YY_1                                # YX
    h = ma.var(Data.data[:, 3, 1, :], 0) / (ma.mean(Data.data[:, 3, 1, :], 0)**2) # YY
    # Get the mean and standard deviation [sigma].
    means = sp.array([ma.mean(a), ma.mean(b), ma.mean(c), ma.mean(d),
                        ma.mean(e), ma.mean(f), ma.mean(g), ma.mean(h)]) 
    sig = sp.array([ma.std(a), ma.std(b), ma.std(c), ma.std(d),
                      ma.std(e), ma.std(f), ma.std(g), ma.std(h)])
    # Get the max accepted value [sigma_thres*sigma, sigma_thres=6 works really well].
    max_sig = sigma_thres*sig
    max_accepted = means + max_sig
    amount_masked = 0
    for freq in range(0, len(a)):
        if ((a[freq] > max_accepted[0]) or
            (b[freq] > max_accepted[1]) or
            (c[freq] > max_accepted[2]) or
            (d[freq] > max_accepted[3]) or
            (e[freq] > max_accepted[4]) or
            (f[freq] > max_accepted[5]) or
            (g[freq] > max_accepted[6]) or
            (h[freq] > max_accepted[7])):
            # mask
            amount_masked += 1
            bad_freq_list.append(freq)
            Data.data[:,:,:,freq].mask = True
    return amount_masked

def destroy_time_with_mean_arrays(Data, flag_size=40):
    '''Mask times with high means.
    
    If there is a problem in time, the mean over all frequencies
    will stand out greatly [>10 sigma has been seen]. Flag these bad
    times and +- `flag_size` times around it. Will only be called if `Data`
    has 'badness'.

    Parameters
    ----------
    Data : DataBlock
        Contains information in a usable format direct from GBT. Bad
        times will be flagged in all polarizations and cal states.
    time_cut : int
        How many frequency bins (as an absolute number) to flag in time.
    '''
    # Loop over beams
    for ii in range(Data.dims[1]):
        # Get the means over all frequencies. (for all pols. and cals.)
        a = ma.mean(Data.data[:, ii, 0, 0, :], -1)
        b = ma.mean(Data.data[:, ii, 1, 0, :], -1)
        # Get means and std for all arrays.
        means = sp.array([ma.mean(a), ma.mean(b)])
        sig = sp.array([ma.std(a), ma.std(b)])
        # Get max accepted values.
        max_accepted = means + 3*sig
        # Find bad times.
        bad_times = []
        for time in range(0,len(a)):
            if ((a[time] > max_accepted[0]) or (b[time] > max_accepted[1])):
                bad_times.append(time)
        # Mask bad times and those +- flag_size around.
        for time in bad_times:
            if time-flag_size < 0:
                Data.data[0:(time+flag_size),ii,:,:,:].mask = True
            else:
                Data.data[(time-flag_size):(time+flag_size),ii,:,:,:].mask = True
    return

def destroy_time_with_mean_arrays_4pol2cal(Data, flag_size=40):
    '''Mask times with high means.
    
    If there is a problem in time, the mean over all frequencies
    will stand out greatly [>10 sigma has been seen]. Flag these bad
    times and +- `flag_size` times around it. Will only be called if `Data`
    has 'badness'.

    Parameters
    ----------
    Data : DataBlock
        Contains information in a usable format direct from GBT. Bad
        times will be flagged in all polarizations and cal states.
    time_cut : int
        How many frequency bins (as an absolute number) to flag in time.
    '''
    # Get the means over all frequencies. (for all pols. and cals.)
    a = ma.mean(Data.data[:, 0, 0, :], -1)
    b = ma.mean(Data.data[:, 1, 0, :], -1)
    c = ma.mean(Data.data[:, 2, 0, :], -1)
    d = ma.mean(Data.data[:, 3, 0, :], -1)
    e = ma.mean(Data.data[:, 0, 1, :], -1)
    f = ma.mean(Data.data[:, 1, 1, :], -1)
    g = ma.mean(Data.data[:, 2, 1, :], -1)
    h = ma.mean(Data.data[:, 3, 1, :], -1)
    # Get means and std for all arrays.
    means = sp.array([ma.mean(a), ma.mean(b), ma.mean(c), ma.mean(d),
                        ma.mean(e), ma.mean(f), ma.mean(g), ma.mean(h)])
    sig = sp.array([ma.std(a), ma.std(b), ma.std(c), ma.std(d),
                      ma.std(e), ma.std(f), ma.std(g), ma.std(h)])
    # Get max accepted values.
    max_accepted = means + 3*sig
    # Find bad times.
    bad_times = []
    for time in range(0,len(a)):
        if ((a[time] > max_accepted[0]) or
            (b[time] > max_accepted[1]) or
            (c[time] > max_accepted[2]) or
            (d[time] > max_accepted[3]) or
            (e[time] > max_accepted[4]) or
            (f[time] > max_accepted[5]) or
            (g[time] > max_accepted[6]) or
            (h[time] > max_accepted[7])):
            bad_times.append(time)
    # Mask bad times and those +- flag_size around.
    for time in bad_times:
        Data.data[(time-flag_size):(time+flag_size),:,:,:].mask = True
    return

def filter_foregrounds(Data, n_bands=20, time_bins_smooth=10.):
    """Gets an estimate of the foregrounds and subtracts it out of the data.
    
    The Foreground removal is very rough, just used to push the foreground down
    a bunch so the RFI can be more easily found.
    Two things are done to estimate the foregrounds: averaging over a fairly
    wide band, and smoothing to just below the beam crossing time scale.

    Parameters
    ----------
    Data : DataBolock object
        Data from which to remove the foregrounds.
    n_bands : int
        Number of bands to split the data into.  Forgrounds are assumed to
        be the same throughout this band.
    time_bins : float
        Number of time bins to smooth over to find the foregrounds (full width
        half max of the filter kernal). Should be
        shorter than the beam crossing time (by about a factor of 2).
    """
    
    # Some basic numbers.
    for jj in range(Data.dims[1]):
        n_chan = Data.dims[-1]
        sub_band_width = float(n_chan)/n_bands
        # First up, initialize the smoothing kernal.
        width = time_bins_smooth/2.355
        # Two sigma edge cut off.
        nk = round(4*width) + 1
        smoothing_kernal = sig.gaussian(nk, width)
        smoothing_kernal /= sp.sum(smoothing_kernal)
        smoothing_kernal.shape = (nk, 1, 1)
        # Now loop through the sub-bands. Foregrounds are assumed to be identical
        # within a sub-band.
        for subband_ii in range(n_bands):
            # Figure out what data is in this subband.
            band_start = round(subband_ii * sub_band_width)
            band_end = round((subband_ii + 1) * sub_band_width)
            data = Data.data[:,jj,:,:,band_start:band_end]
            # Estimate the forgrounds.
            # Take the band mean.
            foregrounds = ma.mean(data, -1)
            # Now low pass filter.
            fore_weights = (sp.ones(foregrounds.shape, dtype=float)
                            - ma.getmaskarray(foregrounds))
            foregrounds -= ma.mean(foregrounds, 0)
            foregrounds = foregrounds.filled(0)
            foregrounds = sig.convolve(foregrounds, smoothing_kernal, mode='same')
            fore_weights = sig.convolve(fore_weights, smoothing_kernal,
                                        mode='same')
            foregrounds /= fore_weights
            # Subtract out the foregrounds.
            data[...] -= foregrounds[:,:,:,None]


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

