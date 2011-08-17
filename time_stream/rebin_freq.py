#!/usr/bin/python
"""This module down samples the frequency axis of the data.

This fails for if your numpy version is < 1.3.
"""

import math

import scipy as sp
import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce

class RebinFreq(base_single.BaseSingle) :
    """Pipeline module that down samples the frequnecy axis of the data."""

    prefix = 'rf_'
    params_init = {
                   'channel_width' : 1.0, # MHz
                   'mean_instead_median' : False,
                   # If 0, rebin by width using `channel_width`.  Otherwise
                   # ignore `chanel_width` and combine this many bins.
                   'n_bins_combined' : 0
                   }
    def action(self, Data, ) :
        if self.params['n_bins_combined'] :
            rebin(Data, self.params['n_bins_combined'], 
                  mean=self.params['mean_instead_median'], by_nbins=True)
            Data.add_history('Rebinned Frequency axis.', 
                             ('Number of bins averaged: '
                                     + str(self.params['n_bins_combined']), ))
        else :
            rebin(Data, self.params['channel_width'], 
                  mean=self.params['mean_instead_median'])
            Data.add_history('Rebinned Frequency axis.', ('Channel width: '
                                     + str(self.params['channel_width']), ))
        return Data

def rebin(Data, width, mean=False, by_nbins=False) :
    """The function that acctually does the rebinning on a Data Block."""
    
    # Input tells us whether to use mean or median.
    if mean :
        method = ma.mean
    else :
        method = ma.median

    if by_nbins :
        width = int(width)
        if width <= 1 :
            raise ValueError("Invalid number of bins to average")
        # Get new axis parameters.
        new_cdelt = width*Data.field['CDELT1']
        nbins = int(sp.ceil(float(Data.dims[-1])/width))
        new_centre = nbins//2 + 1
        Data.calc_freq()
        Data.field['CRVAL1'] = Data.freq[int((new_centre+0.5)*width)]
        # Case where evenly divisable (much more efficient).
        if Data.dims[-1] % width == 0:
            new_data = Data.data
            new_data.shape = Data.data.shape[:-1] + (nbins, width)
            new_data = method(new_data, -1)
        else :
            # Allowcate memory for Data array.
            new_data = ma.empty(Data.dims[:3] + (nbins,))
            # Loop over new bins and rebin.
            for ii in xrange(nbins) :
                new_data[:,:,:,ii] = method(Data.data[:,:,:,ii*width:(ii+1)*width],
                                        3)
        Data.set_data(new_data)
    else :
        # Convert to Hertz.
        width = width*1.0e6
        new_cdelt = width * sp.sign(Data.field['CDELT1'])
        # Figure out some basics.
        Data.calc_freq()
        freq = sp.array(Data.freq)
        # Extra bit on the bandwidth is because frequency labels are channel 
        # centre.
        bandwidth = abs(freq[-1] - freq[0]) + abs(Data.field['CDELT1'])
        nbins = int(bandwidth//width)
        new_centre = int((Data.field['CRPIX1']-1)
                         * abs(Data.field['CDELT1'])/width)
        new_dims = Data.dims[0:-1] + (nbins, )
        # Get old data and allowcate memory for new data.
        old_data = ma.array(Data.data, copy=True)
        Data.set_data(ma.zeros(new_dims))
        new_freq = Data.field['CRVAL1'] + new_cdelt*(sp.arange(nbins)
                                                     - new_centre)
        for ii in range(1,nbins-1) :
            inds = (sp.logical_and(
                        abs(freq - new_freq[ii]) <= abs(freq - new_freq[ii+1]),
                        abs(freq - new_freq[ii]) < abs(freq - new_freq[ii-1])))
            subdata = (old_data[:,:,:,inds])
            Data.data[:,:,:,ii] = method(subdata, 3)
        # Above loop breaks for end points... deal with them.
        inds, = sp.where(abs(freq - new_freq[0]) <= abs(freq - new_freq[1]))
        subdata = old_data[:,:,:,inds]
        Data.data[:,:,:,0] = method(subdata, 3)
        inds, = sp.where(abs(freq-new_freq[nbins-1])
                         < abs(freq-new_freq[nbins-2]))
        subdata = old_data[:,:,:,inds]
        Data.data[:,:,:,nbins-1] = method(subdata, 3)
        
        Data.freq = new_freq
    Data.field['CRPIX1'] = sp.array(new_centre + 1, dtype=int)
    Data.field['CDELT1'] = sp.array(new_cdelt, dtype=float)



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    RebinFreq(str(sys.argv[1])).execute()

