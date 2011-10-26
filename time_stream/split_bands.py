#!/usr/bin/python
"""This module down samples the frequency axis of the data.

This fails for if your numpy version is < 1.3.
"""

import math

import scipy as sp
import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce
from core import utils

prefix = 'rt_'

params_init = {
        }

class RebinTime(base_single.BaseSingle) :
    """Pipeline module that down samples the frequnecy axis of the data."""

    def action(self, Data, ) :
        split(Data, )
        Data.add_history('Split into bands.', 
                         ('Number of time bins averaged: '
                          + str(self.params['n_bins_combined']), ))
        return Data
