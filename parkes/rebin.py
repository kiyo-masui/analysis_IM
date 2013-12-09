#!/usr/bin/python
"""This module down samples the time axis of the data.
"""

import math

import scipy as sp
import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce
#import utils.misc as utils
import misc as utils
import rebin_time
import rebin_freq


class Rebin(base_single.BaseSingle) :
    """Pipeline module that down samples the frequnecy axis of the data."""

    params_init = {
                   'n_bins_combined_time' : 2,
                   'n_bins_combined_freq' : 16
                   }
    prefix = 'rb_'

    def action(self, Data, ) :
        if self.params['n_bins_combined_freq'] > 1:
            rebin_freq.rebin(Data, self.params['n_bins_combined_freq'],
                             True, True)
            Data.add_history('Rebinned freq axis.', 
                             ('Number of freq bins averaged: '
                              + str(self.params['n_bins_combined_freq']), ))
        if self.params['n_bins_combined_time'] > 1:
            rebin_time.rebin(Data, self.params['n_bins_combined_time'])
            Data.add_history('Rebinned time axis.', 
                             ('Number of time bins averaged: '
                              + str(self.params['n_bins_combined_time']), ))
        return Data



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Rebin(str(sys.argv[1])).execute()

