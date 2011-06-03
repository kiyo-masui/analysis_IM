#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import os
import copy

import numpy.ma as ma
import scipy as sp
import scipy.signal as sig
import matplotlib.pyplot as plt

import core.fitsGBT
import kiyopy.custom_exceptions as ce
import base_single
import hanning
import cal_scale
from time_stream import rotate_pol

class FlagData(base_single.BaseSingle) :
    """Pipeline module that flags rfi and other forms of bad data.

    For lots of information look at the doc-string for 
    flag_data_aravind.apply_cuts.
    """

    prefix = 'fd_'
    params_init = {
                   # In multiples of the standard deviation of the whole block
                   # once normalized to the time median.
                   'perform_hanning' : False,
                   'cal_scale' : False,
                   # Rotate to XX,XY,YX,YY is True.
                   'rotate' : False,
                   # A Data that has more than badness_thres frequencies flagged
                   # (as a %) will be considered bad.
                   'badness_thres' : 0.1,
                   # How many times to hide around a bad time.
                   'time_cut' : 40
                   }
    feedback_title = 'New flags each data Block: '
    
    def action(self, Data) :
        print 'i am here'
        params = self.params
        # Keep track of how many pre existing flags there are for feedback
        # purposes.
        already_flagged = ma.count_masked(Data.data)
        if params["rotate"]:
            if (tuple(Data.field['CRVAL4']) == (1, 2, 3, 4)):
                print 'switched'
                rotate_pol.rotate(Data, (-5,-7,-8,-6))
        # Few operations to be performed before flagging.
        if params["perform_hanning"] :
            hanning.hanning_smooth(Data)
            Data.add_history('Hanning smoothed.')
        if params["cal_scale"] :
            cal_scale.scale_by_cal(Data, True, False, False)
            Data.add_history('Converted to units of noise cal temperture.')
        # Flag the data.
        print Data.field['CRVAL4']
        apply_cuts(Data, badness_thres=params['badness_thres'],
                    time_cut=params['time_cut'])
        Data.add_history('Flagged Bad Data.', ('Rotated to XX,XY,YX,YY: ' + 
                     str(self.params['rotate']))) #, ('Sigma threshold: ' +
           #         str(self.params['sigma_thres']), 'Polarization threshold: '
           #         + str(self.params['pol_thres'])))
        # Report the number of new flags.
        new_flags = ma.count_masked(Data.data) - already_flagged
        self.block_feedback = str(new_flags) + ', '
        print Data.field['CRVAL4']
        return Data

def apply_cuts(Data, badness_thres=0.1, time_cut=40):
    """Flags bad data from RFI and far outliers.
    See flag_data for other parameters."""
    # Make Data XX,XY,YX,YY only if in Stokes' parameters.
#    switched = False
#    if (tuple(Data.field['CRVAL4']) == (1, 2, 3, 4)):
#        print 'switched'
#        rotate_pol.rotate(Data, (-5,-7,-8,-6))
#        switched = True
    # Flag data and store if it was bad or not [see determine badness].
    badness = flag_data(Data, badness_thres, time_cut)
    print badness
#    # Switch Data back to I,Q,U,V if they were so originally.
#    if switched:
#        rotate_pol.rotate(Data, (1,2,3,4))
    return


def flag_data(Data, badness_thres, time_cut):
    '''Flag bad data from RFI and far outliers. See params_init dictionary
    for badness_thres. See 'flag_size' in destroy_time_with_mean_arrays
    for time_cut.'''
    # Flag data on a [deep]copy of Data. If too much destroyed,
    # check if localized in time. If that sucks too, then just hide freq.
    Data1 = copy.deepcopy(Data)
    itr = 0            # For recursion
    max_itr = 20       # For recursion
    bad_freqs = []
    amount_masked = -1 # For recursion
    while not (amount_masked == 0) and itr < max_itr:                         
        amount_masked = destroy_with_variance(Data1, bad_freq_list=bad_freqs) 
        itr += 1
    bad_freqs.sort()
    print bad_freqs
    # Remember the flagged data.
    mask = Data1.data.mask
##    badness = determine_badness(bad_freqs)
    percent_masked1 = (float(len(bad_freqs)) / Data1.dims[-1])
    badness = percent_masked1 > badness_thres
    # If too many frequencies flagged, it may be that the problem
    # happens in time, not in frequency.
    if badness:
        Data2 = copy.deepcopy(Data)
        #flag_across_time(Data2)
        destroy_time_with_mean_arrays(Data2, flag_size=time_cut)
        # Bad style for repeating as above, sorry.
        itr = 0
        bad_freqs = []
        amount_masked = -1
        while not (amount_masked == 0) and itr < max_itr:
            amount_masked = destroy_with_variance(Data2, bad_freq_list=bad_freqs) 
            itr += 1
        bad_freqs.sort()
##        badness = determine_badness(bad_freqs)
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

def destroy_with_variance(Data, level=1, bad_freq_list=[]):
    '''Mask spikes in Data using variance. Polarizations must be in
    XX,XY,YX,YY format.
    level represents how sensitive the flagger is (smaller = more masking).
    The flagged frequencies are appended to bad_freq_list.'''
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
    # Get the max accepted value [6*sigma, with level thrown in for fun].
    max_sig = 6*sig*level
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
#           for pols in range(0, Data.dims[1]):
#                Data.data[:, pols, 0, freq].mask = True
#                Data.data[:, pols, 1, freq].mask = True
            Data.data[:,:,:,freq].mask = True
    print amount_masked
    return amount_masked

def determine_badness(bad_freqs):
    '''Return True if too many frequencies have been masked [hence, bad].
    If True, the data must be flagged in time, not frequency.
    "Too many frequencies" means that there exists a slice of bad_freqs
    that contains at least 200 elements and the difference between adjacent
    elements is less than 10. Something like this would happen if the 
    telescope goes wonky for a few seconds.
    bad_freqs is a SORTED list of integers representing flagged frequencies.'''
    count = 0
    high_score = 0
    if len(bad_freqs) > 200:
        for i in range(1, len(bad_freqs)):
            if ((bad_freqs[i] - bad_freqs[i-1]) < 10):
                count += 1
            else:
                high_score = count
                count = 0
    badness = (high_score >= 200)
    return badness

def flag_across_time(Data, section_size=15):
    '''Will only be called if a Data has "badness" [see determine_badness].
    Splits the Data into different pieces in time, section_size long, and
    find the badness of each. The sections in time that have badness will be
    flagged, so that [hopefully] the Data will no longer have badness overall.'''
    bad_times = []
    for splits in range(0, Data.dims[0], section_size):
        # Make an array that mimics a Data block.
        arr = ma.masked_array([Data.data[splits:(splits+section_size), 0, 0, :],
            Data.data[splits:(splits+section_size), 1, 0, :],
            Data.data[splits:(splits+section_size), 2, 0, :],
            Data.data[splits:(splits+section_size), 3, 0, :],
            Data.data[splits:(splits+section_size), 0, 1, :],
            Data.data[splits:(splits+section_size), 1, 1, :],
            Data.data[splits:(splits+section_size), 2, 1, :],
            Data.data[splits:(splits+section_size), 3, 1, :]], mask=False)
        # If there is masked data already in the split, it will ruin
        # accuracy so all data from there is ignored.
        if (True in Data.data[splits:(splits+section_size),:,:,:].mask):
            bad_times.append(splits) 
        itr = 0            # For recursion
        max_itr = 20       # For recursion
        bad_freqs = []
        amount_masked = -1 # For recursion
        while not (amount_masked == 0) and itr < max_itr:                         
            amount_masked = destroy_with_variance_arrays(arr, bad_freq_list=bad_freqs) 
            itr += 1
        bad_freqs.sort()
##        badness = determine_badness(bad_freqs)
#        badness = (float(len(bad_freqs)) / Data.dims[-1]) > 0.05
#        if badness:
        # If section in time is bad, remember it.
        if len(bad_freqs) > 50:
            bad_times.append(splits)
    # Flag Data in time for all bad time sections found.
    print bad_times
    for time in bad_times:
#        for pols in range(0, Data.dims[1]):
#            Data.data[time:(time+section_size), pols, 0, :].mask = True
#            Data.data[time:(time+section_size), pols, 1, :].mask = True
        Data.data[time:(time+section_size),:,:,:].mask = True
    return

def destroy_with_variance_arrays(arr, level=1, bad_freq_list=[]):
    '''Same as "destroy_with_variance" but for arrays arr, not Data.
    arr should contain 8 arrays for each XX,XY,YX,YY polarization
    [in that order] with cal=0, then the same again with cal=1.'''
    XX_YY_0 = ma.mean(arr[0], 0) * ma.mean(arr[3], 0)
    XX_YY_1 = ma.mean(arr[4], 0) * ma.mean(arr[7], 0)
    a = ma.var(arr[0], 0) / (ma.mean(arr[0], 0)**2) # XX
    b = ma.var(arr[1], 0) / XX_YY_0                 # XY
    c = ma.var(arr[2], 0) / XX_YY_0                 # YX
    d = ma.var(arr[3], 0) / (ma.mean(arr[3], 0)**2) # YY
    e = ma.var(arr[4], 0) / (ma.mean(arr[4], 0)**2) # XX
    f = ma.var(arr[5], 0) / XX_YY_1                 # XY
    g = ma.var(arr[6], 0) / XX_YY_1                 # YX
    h = ma.var(arr[7], 0) / (ma.mean(arr[7], 0)**2) # YY
    means = sp.array([ma.mean(a), ma.mean(b), ma.mean(c), ma.mean(d),
                        ma.mean(e), ma.mean(f), ma.mean(g), ma.mean(h)])
    sig = sp.array([ma.std(a), ma.std(b), ma.std(c), ma.std(d),
                      ma.std(e), ma.std(f), ma.std(g), ma.std(h)])
    max_sig = 6*sig*level
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
            arr.mask[:,:,freq] = True
    print amount_masked
    return amount_masked

def destroy_time_with_mean_arrays(Data, flag_size=40):
    '''If there is a problem in time, the mean over all frequencies
    will stand out greatly [>10 sigma has been seen]. Flag these bad
    times and +- flag_size times around it. Will only be called if a Data has
    "badness" [see determine_badness].'''
    # Get the means over all frequencies.
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
    print bad_times
    for time in bad_times:
        Data.data[(time-flag_size):(time+flag_size),:,:,:].mask = True
    return

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

