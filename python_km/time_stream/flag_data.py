#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import os

import numpy.ma as ma
import scipy as sp
import scipy.signal as sig

import core.fitsGBT
import kiyopy.custom_exceptions as ce
import base_single
import hanning

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
                   'sigma_thres' : 5,
                   # In multiples of the measured standard deviation.
                   'pol_thres' : 5
                   }
    feedback_title = 'New flags each data Block: '
    
    def action(self, Data) :
        already_flagged = ma.count_masked(Data.data)
        params = self.params
        if params["perform_hanning"] :
            hanning.hanning_smooth(Data)
            Data.add_history('Hanning smoothed.')
        apply_cuts(Data, sig_thres=params['sigma_thres'], 
                   pol_thres=params['pol_thres'])
        new_flags = ma.count_masked(Data.data) - already_flagged
        self.block_feedback = str(new_flags) + ', '

        Data.add_history('Flagged Bad Data.', ('Sigma threshold: ' +
                    str(self.params['sigma_thres']), 'Polarization threshold: '
                    + str(self.params['pol_thres'])))
        return Data

def apply_cuts(Data, sig_thres=5.0, pol_thres=5.0) :
    """Flags bad data from RFI and far outliers..
    
    Masks any data that is further than 'sig_thres' sigmas from the mean of the
    entire data block (all frequencies, all times).  Each cal, fequency and
    polarization is first normalized to the time mean.  Cut only checked for XX
    and YY polarizations but all polarizations are masked.  This cut is
    deactivated by setting 'sig_thres' < 0.

    Aravind rfi cut is also performed, cutting data with polarization lines, as
    well as lines in XX and YY.

    Arguments:
        Data -      A DataBlock object containing data to be cleaned.
        sig_thres - (float) Any XX or YY data that deviates by this more than 
                    this many sigmas is flagged.  Data first normalized to 
                    time median.
        pol_thres - (float) Any data cross polarized by more than this 
                    many sigmas is flagged.
        width -     (int) In the polaridezation cut, flag data within this many
                    frequency bins of offending data.
        flatten -   (Bool) Preflatten the polarization spectrum.
    """

    # Here we check the polarizations and cal indicies
    if tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6) :
        xx_ind = 0
        yy_ind = 3
        pol_inds = [1,2]
        norm_inds = [0, 3]
    elif tuple(Data.field['CRVAL4']) == (1, 2, 3, 4) :
        pol_inds = [2,3]
        xx_ind = 0
        yy_ind = 1
        norm_inds = [0,0]
    else :
        raise ce.DataError('Polarization types not as expected,'
                           ' function needs to be generalized.')
    if tuple(Data.field['CAL']) == ('T', 'F') :
        on_ind = 0
        off_ind = 1
    else :
        raise ce.DataError('Cal states not as expected.')
    if pol_thres > 0 :    
        # Always flag cal on and cal off together. Don't use ma.sum because we
        # want to transfer flags.
        data = Data.data[:,:,0,:] + Data.data[:,:,1,:]
        # Calculate the polarization cross coefficient.
        cross = ma.sum(data[:,pol_inds,:]**2, 1)
        #cross /= (data[:,norm_inds[0],:] *
        #          data[:,norm_inds[1],:]).filled(1.0e-10)
        # Run it throught the flagger.
        mask = filter_flagger(cross, Data.dims[-1]//200, pol_thres, axis=-1)
        for ii in range(Data.dims[1]) :
            Data.data[:, ii, 0, :][mask] = ma.masked
            Data.data[:, ii, 1, :][mask] = ma.masked
    if sig_thres > 0 :
        data = Data.data[:,:,0,:] + Data.data[:,:,1,:]
        # Figure out the order of the polynomial by the range of azimuths.
        beam_width = 0.3 # Approximate.
        order = round(abs(max(Data.field['CRVAL2']) - min(Data.field['CRVAL2']))
                 / beam_width)
        order = max([order, 1])
        # First flag the frequency average, looking for obviouse blips.
        tdata = ma.mean(data/ma.mean(data, 0), -1)
        for ii in range(Data.dims[1]) :
            mask = filter_flagger(tdata[:, [ii]], order, 
                                  sig_thres, axis=0)
            Data.data[mask, ...] = ma.masked
        # Now flag at each frequency.
        data = Data.data[:,:,0,:] + Data.data[:,:,1,:]
        for ii in range(Data.dims[1]) :
            this_mask = filter_flagger(data[:, ii, :], order,
                                  sig_thres, axis=0)
            if ii == 0:
                mask = this_mask
            else :
                mask[this_mask] = True
        for ii in range(Data.dims[1]) :
            Data.data[:, ii, 0, :][mask] = ma.masked
            Data.data[:, ii, 1, :][mask] = ma.masked

def filter_flagger(arr, order, threshold, axis=-1) :
    """arr should be 2 d array, flagging done along the `axis` dimension."""
    
    # Maximum allowible number of iterations.
    max_itr = 20
    if arr.ndim == 1 :
        arr = sp.reshape(arr, arr.shape + (1,))
        axis = 0
    elif arr.ndim == 2 :
        pass
    else :
        raise ValueError("Input array must be 1 or 2 dimensional.")
    # Bring axis to the first index.
    arr = sp.rollaxis(arr, axis)
    # Allowcate memory for the mask array.
    mask = sp.zeros(arr.shape, dtype=bool)
    # Convert to a normal array.
    if isinstance(arr, ma.MaskedArray) :
        mask[...] = arr.mask[...]
        arr = arr.filled(0.0)
    # We make progressivly tighter and tighter cuts, starting with the very
    # weak sqrt(n/2)*stdev.
    n = arr.shape[0]
    coarse_thres = max(sp.sqrt(n/2), threshold)
    thresholds = [coarse_thres, sp.sqrt(coarse_thres*threshold), threshold]
    thres_ind = 0
    # Initialize a few intermediate arrays.
    tmp_arr = arr.copy()
    smoothed = sp.empty(arr.shape, dtype=float)
    smoothed[...] = sp.mean(tmp_arr, 0)
    # 'x' values for the polyfit.
    x = sp.arange(n)
    # Iteritivly flag data.
    for ii in range(max_itr) :
        # Get rid of bad data.
        tmp_arr[mask] = smoothed[mask]
        # Smooth the data.
        p = sp.polyfit(x, tmp_arr, order)
        for jj in range(p.shape[1]) :
            smoothed[:, jj] = sp.polyval(p[:, jj], x)
        # Flag the data.
        new_masks = (abs(tmp_arr - smoothed)
                     > thresholds[thres_ind]*sp.std(tmp_arr - smoothed, 0))
        if sp.all(sp.sum(new_masks, 0) == 0) :
            thres_ind += 1
            if thres_ind == len(thresholds) :
                break
        else :
            mask[new_masks] = True
    if axis == 1 or axis == -1 :
        mask = sp.rollaxis(mask, 0, 2)
    return mask

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

