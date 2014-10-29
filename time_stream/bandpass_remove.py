#! /usr/bin/python
"""This module remove the bandpass gain
"""

import os
import copy
import gc

import numpy as np
import numpy.ma as ma
import scipy as sp
import scipy.signal as sig

from core import fitsGBT
import base_single
from utils import batch_handler

import kiyopy.custom_exceptions as ce
from time_stream import hanning
from time_stream import cal_scale
from time_stream import rotate_pol

# XXX
#import matplotlib.pyplot as plt

class BandpassRemove(base_single.BaseSingle) :
    '''Pipeline module that flags RFI and other forms of bad data.

    '''

    prefix = 'bp_'
    params_init = {
                   # In multiples of the standard deviation of the whole block
                   # once normalized to the time median.
                   'perform_hanning' : False,
                   'time_cut' : 60,
                   }
    feedback_title = 'New flags each data Block: '

    def action(self, Data):
        '''Prepares Data and flags RFI.
        
        Parameters
        ----------
        Data : DataBlock
            Contains information in a usable format direct from GBT. 

        Returns
        -------
        Data : DataBlock
            The input `Data` with RFI flagged. Will also be cal scaled and
            rotated to XX,YY... if so chosen.

        '''
        params = self.params
        bandpass_rm(Data, time_n=params['time_cut'])
        Data.add_history('Flagged Bad Data.',('time cut:'+str(params['time_cut']),))
        if params["perform_hanning"] :
            hanning.hanning_smooth(Data)
            Data.add_history('Hanning smoothed.')
        return Data

@batch_handler.log_timing
def bandpass_rm(Data, time_n=60):
    '''bandpass remove

    Parameters
    ----------
    Data : DataBlock
        Contains information in a usable format direct from GBT. Bad
        frequencies will be flagged in all polarizations and cal states.
    time_cut : int
        How many time bins (as an absolute number) are used for estimate bandpass.

    '''
    # get the tsys
    tsys = Data.field['TSYS']

    # get the data
    spec = Data.data

    # flag by tsys
    flag_times = 20
    for i in range(flag_times):
        spec, tsys = flag_by_tsys(spec, tsys)
    spec = np.array(spec)
    tsys = np.array(tsys)

    spec[np.isnan(spec)] = ma.masked
    tsys[np.isnan(tsys)] = ma.masked

    shape = spec.shape

    spec_m = np.zeros(spec.shape)
    tsys_m = np.zeros(tsys.shape)

    for i in range(shape[0]):
        # select a time range for bandpass estimation 
        if i > time_n+1:
            spec_i = spec[i-time_n-1:i+time_n+1,...]
            tsys_i = tsys[i-time_n-1:i+time_n+1,...]
        else:
            spec_i = spec[0:i+time_n+1,...]
            tsys_i = tsys[0:i+time_n+1,...]
        # mask out the current time
        if i > 1: 
            spec_i = np.delete(spec_i, np.s_[i-1:i+1], 0)
            tsys_i = np.delete(tsys_i, np.s_[i-1:i+1], 0)
        else:
            spec_i = np.delete(spec_i, np.s_[0:i+1], 0)
            tsys_i = np.delete(tsys_i, np.s_[0:i+1], 0)

        spec_m[i,...] = np.median(spec_i, axis=0)
        tsys_m[i,...] = np.median(tsys_i, axis=0)

        del spec_i
        del tsys_i
        gc.collect()

        #spec_m[:,i,:,:,:] = ma.median(spec_i, axis=1)
        #tsys_m[:,i,:,:] = ma.median(tsys_i, axis=1)

    spec = tsys_m[:,:,None,None] * (spec/spec_m) - tsys[:,:,None,None]

    spec = np.ma.array(spec)
    spec[spec==0] = np.ma.masked
    spec[np.isnan(spec)] = np.ma.masked

    Data.data = spec

def flag_by_tsys(data, tsys, sig=3):
    data = ma.array(data)
    tsys = ma.array(tsys)
    tsys[np.isnan(tsys)] = ma.masked
    tsys[tsys==0] = ma.masked
    tsys_mean = ma.mean(tsys, axis=0)
    tsys_std  = ma.std(tsys, axis=0)
    #print tsys_std.T
    tsys[tsys>tsys_mean[None,...]+sig*tsys_std[None,...]] = ma.masked
    tsys[tsys<tsys_mean[None,...]-sig*tsys_std[None,...]] = ma.masked

    data[tsys.mask,:] = ma.masked

    tsys = tsys.filled(0)
    data = data.filled(0)

    return data, tsys


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

