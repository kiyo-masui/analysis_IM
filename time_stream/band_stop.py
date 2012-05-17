#!/usr/bin/python
"""This module masks frequency ranges."""

import math

import scipy as sp
import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce
import utils.misc as utils
import rebin_time
import rebin_freq


class BandStop(base_single.BaseSingle) :
    """Pipeline module that down samples the frequnecy axis of the data."""

    params_init = {
                   'band_stops' : []
                   }
    prefix = 'bs_'

    def action(self, Data):
        
        for limits in self.params['band_stops']:
            mask_frequency_range(Data, limits)
        
        #Data.add_history('Masked stop bands.',
        #                 ('Stop bands: ' + str(self.params['band_stops'])))
        # Less detailed history message until pyfits bug gets fixed.
        Data.add_history('Masked stop bands.')
        return Data

def mask_frequency_range(Data, limits):

    if len(limits) != 2:
        raise ValueError("Limits must be length 2 sequence.")
    lower = min(limits)
    upper = max(limits)

    Data.calc_freq()
    freq = Data.freq
    delta_f = abs(sp.mean(sp.diff(freq)))
    mask = sp.logical_and(freq + delta_f/2 <= upper, freq - delta_f/2 >= lower)
    Data.data[...,mask] = ma.masked

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    BandStop(str(sys.argv[1])).execute()

