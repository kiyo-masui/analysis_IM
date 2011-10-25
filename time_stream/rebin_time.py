#!/usr/bin/python
"""This module down samples the frequency axis of the data.

This fails for if your numpy version is < 1.3.
"""

import math

import scipy as sp
import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce

prefix = 'rt_'

class RebinTime(base_single.BaseSingle) :
    """Pipeline module that down samples the frequnecy axis of the data."""

    params_init = {
                   'n_bins_combined' : 2
                   }
    def action(self, Data, ) :
        rebin(Data, self.params['n_bins_combined'])
        Data.add_history('Rebinned time axis.', 
                         ('Number of time bins averaged: '
                          + str(self.params['n_bins_combined']), ))
        return Data

def rebin(Data, n_bins) :
    """The function that acctually does the rebinning on a Data Block."""
    pass
    



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    RebinTime(str(sys.argv[1])).execute()

