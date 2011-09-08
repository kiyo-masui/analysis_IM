#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import os
import copy

import numpy.ma as ma
import scipy as sp
import scipy.signal as sig

import core.fitsGBT
import kiyopy.custom_exceptions as ce
import base_single
import hanning
import cal_scale
from time_stream import rotate_pol

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
                   'sigma_thres' : 6,
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
        self.block_feedback = str(new_flags) + ', '
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
            mask = Data2.data.mask
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

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

