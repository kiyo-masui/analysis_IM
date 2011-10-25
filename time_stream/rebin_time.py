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

def rebin(Data, n_bins_combined) :
    """The function that acctually does the rebinning on a Data Block."""
    
    nt = Data.data.shape[0]
    new_nt = nt // n_bins_combined
    new_shape = (new_nt,) + Data.data.shape[1:]
    unmask = sp.logical_not(ma.getmaskarray(Data.data))
    data = Data.data.filled(0)
    # Allowcate memeory for the rebinned data.
    new_data = ma.zeros(new_shape, dtype=data.dtype)
    counts = sp.zeros(new_shape, dtype=int)
    # Add up the bins to be combined.
    for ii in range(n_bins_combined):
        new_data += data[ii:new_nt * n_bins_combined:n_bins_combined,...]
        counts += unmask[ii:new_nt * n_bins_combined:n_bins_combined,...]
    new_data[counts == 0] = ma.masked
    counts[counts == 0] = 1
    new_data /= counts
    Data.set_data(new_data)
    # Now deal with all the other records that aren't the main data.
    for field_name in Data.field.iterkeys():
        # DATE-OBS is a string field so we have to write special code for it.
        if field_name == "DATE-OBS":
            time_field = Data.field[field_name]
            new_field = sp.empty(new_nt, dtype=Data.field[field_name].dtype)
            # Convert to float, average, then convert back to a string.
            time_float = utils.time2float(time_field)
            for ii in range(new_nt):
                tmp_time = sp.mean(time_float[n_bins_combined * ii
                                              : n_bins_combined * (ii + 1)])
                new_field[ii] = utils.float2time(tmp_time)
            Data.set_field(field_name, new_field, 
                       axis_names=Data.field_axes[field_name],
                       format=Data.field_formats[field_name])
            continue
        # Only change fields that have a 'time' axis.
        try:
            time_axis = list(Data.field_axes[field_name]).index('time')
        except ValueError:
            continue
        # For now, the time axis has to be the first axis.
        if time_axis != 0:
            msg = "Expected time to be the first axis for all fields."
            raise NotImplementedError(msg)
        field_data = Data.field[field_name]
        if not field_data.dtype.name == "float64":
            msg = "Field data type is not float. Handle explicitly."
            raise NotImplementedError(msg)
        new_field = sp.empty(field_data.shape[:time_axis] + (new_nt,) 
                             + field_data.shape[time_axis + 1:],
                             dtype=field_data.dtype)
        for ii in range(new_nt):
            tmp_data = sp.sum(field_data[n_bins_combined * ii
                                         :n_bins_combined * (ii + 1),...], 0)
            tmp_data /= n_bins_combined
            new_field[ii,...] = tmp_data
        Data.set_field(field_name, new_field, 
                       axis_names=Data.field_axes[field_name],
                       format=Data.field_formats[field_name])



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    RebinTime(str(sys.argv[1])).execute()

