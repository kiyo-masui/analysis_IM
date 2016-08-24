#!/usr/bin/python
"""This module down samples the time axis of the data.
"""

import math

import scipy as sp
import numpy.ma as ma

from parkes import base_single
from parkes import rebin_time
from parkes import rebin_freq
import kiyopy.custom_exceptions as ce
import utils.misc as utils


class CheckFitsFile(base_single.BaseSingle) :
    """Pipeline module that down samples the frequnecy axis of the data."""

    params_init = {
                   'n_bins_combined_time' : 2,
                   'n_bins_combined_freq' : 16
                   }
    prefix = 'rb_'

    def action(self, Data, ) :
        freq = Data.calc_freq()
        print freq



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    rawdatapath = ('/mnt/data-pen3/ycli/map_result/rebinned/parkes_2008_09_12_west_P641.fits',)
    CheckFitsFile(rawdatapath).execute()

