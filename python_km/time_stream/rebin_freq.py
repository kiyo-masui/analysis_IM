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
            raise NotImplementedError('Right now can only rebin by channel'
                                      ' width.')
        rebin(Data, self.params['channel_width'], 
              mean=self.params['mean_instead_median'])
        Data.add_history('Rebinned Frequency axis.', ('Channel width: '
                                 + str(self.params['channel_width']), ))
        return Data

def rebin(Data, width, mean=False) :
    """The function that acctually does the rebinning on a Data Block."""
    
    # Convert to Hertz (in python this is safe).
    width = width*1.0e6
    new_cdelt = width * round(Data.field['CDELT1']/abs(Data.field['CDELT1']))
    # Figure out some basics.
    Data.calc_freq()
    freq = sp.array(Data.freq)
    # Extra bit on the bandwidth is because frequency labels are channel centre.
    bandwidth = abs(freq[-1] - freq[0]) + abs(Data.field['CDELT1'])
    nbins = int(bandwidth//width)
    new_centre = int((Data.field['CRPIX1']-1)*abs(Data.field['CDELT1'])/width)
    new_dims = Data.dims[0:-1] + (nbins, )
    # Get old data and allowcate memory for new data.
    old_data = ma.array(Data.data, copy=True)
    Data.set_data(ma.zeros(new_dims))
    # Input tells us mean or median.
    if mean :
        method = ma.mean
    else :
        method = ma.median

    new_freq = Data.field['CRVAL1'] + new_cdelt*(sp.arange(nbins) - new_centre)
    for ii in range(1,nbins-1) :
        inds = (sp.logical_and(
                    abs(freq - new_freq[ii]) <= abs(freq - new_freq[ii+1]),
                    abs(freq - new_freq[ii]) < abs(freq - new_freq[ii-1]) ))
        subdata = (old_data[:,:,:,inds])
        Data.data[:,:,:,ii] = method(subdata, 3)
    # Above loop breaks for end points... deal with them.
    inds, = sp.where(abs(freq - new_freq[0]) <= abs(freq - new_freq[1]))
    subdata = old_data[:,:,:,inds]
    Data.data[:,:,:,0] = method(subdata, 3)
    inds, = sp.where(abs(freq-new_freq[nbins-1]) < abs(freq-new_freq[nbins-2]))
    subdata = old_data[:,:,:,inds]
    Data.data[:,:,:,nbins-1] = method(subdata, 3)
    
    Data.freq = new_freq
    Data.field['CDELT1'] = sp.array(new_cdelt, dtype=float)
    Data.field['CRPIX1'] = sp.array(new_centre + 1, dtype=int)



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    RebinFreq(str(sys.argv[1])).execute()

