#!/usr/bin/python
"""Module splits up the frequency axis into multiple blocks (bands, IFs,
whatver you want to call them).
"""

import math
import copy

import scipy as sp
import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce
from core import utils
from core import data_block


class SplitBands(base_single.BaseSingle) :
    """Pipeline module that splits data into multiple bands."""

    prefix = 'sb_'

    params_init = {
            'n_bands' : 2,
            'n_bins_band' : 32,
            'offset' : 0
            }

    def action(self, Data, ):
        params = self.params
        Blocks = split(Data, params['n_bands'], params['n_bins_band'],
                       params['offset'])
        Data.add_history('Split into bands.', 
                         ('(n_bands, n_bins_band, offset): '
                          + str((params['n_bands'], params['n_bins_band'],
                                 params['offset'])),))
        return Blocks


def split(Data, n_bands=2, bins_band=None, offset=0):
    
    n_chan = Data.dims[-1]
    if offset < 0:
        raise ValueError("offset parameter must be greater or equal to 0.")
    if bins_band is None:
        bins_band = (n_chan - offset) // n_bands
    # Check inputs.
    if n_bands * bins_band + offset > n_chan:
        msg = "Requested bands don't fit in input frequency axis."
        raise ValueError(msg)
    Data.calc_freq()
    freq = Data.freq
    out_Blocks = ()
    # Make the Blocks.
    for ii in range(n_bands):
        # Copy the data.
        new_data = Data.data[...,offset + ii * bins_band
                             :offset + (ii + 1) * bins_band]
        this_Block = data_block.DataBlock(new_data)
        # Copy the other fields
        this_Block.field = copy.deepcopy(Data.field)
        this_Block.field_axes = copy.deepcopy(Data.field_axes)
        this_Block.field_formats = copy.deepcopy(Data.field_formats)
        # Fix the frequency axis.
        this_Block.field['CRPIX1'] = sp.array(bins_band//2 + 1, dtype=int)
        this_Block.field['CRVAL1'] = sp.array(freq[offset + ii * bins_band +
                                              bins_band//2], dtype=float)
        # Check the new block and add it to the outputs.
        this_Block.verify()
        out_Blocks += (this_Block,)
    return out_Blocks

        
