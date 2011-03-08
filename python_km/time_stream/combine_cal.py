#!/usr/bin/python
"""Module that combines cal on and cal off data, subtracts out time means and
generally prepares data for map making."""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single


class CombineCal(base_single.BaseSingle) :
    """Pipeline module that does some final manipulation of data.
    """

    prefix = 'cc_'
    params_init = {
                   'subtract_time_mean' : True, 
                   'average_cal_states' : True,
                   # weights should sum to 1.0 (to remain calibrated).
                   'weights' : (0.5, 0.5)
                   }

    def action(self, Data):
        combine(Data, self.params["weights"], 
                sub_mean=self.params['subtract_time_mean'],
                average_cals=self.params['average_cal_states']), 
        Data.add_history('Combined cal-on and cal-off data.')
        return Data

def combine(Data, weights=(0.5, 0.5), sub_mean=True, average_cals=True) :
    """
    """
    
    if average_cals :
        on_ind = 0
        off_ind = 1
        if (Data.field['CAL'][on_ind] != 'T' or
            Data.field['CAL'][off_ind] != 'F') :
                raise ce.DataError('Cal states not in expected order.')
        newdata = (Data.data[:,:,[0],:]*weights[0] +
                   Data.data[:,:,[1],:]*weights[1])
        Data.set_data(newdata)
        Data.field['CAL'] = sp.array(['A']) # For averaged.

    # Subtract the time median if desired.
    if sub_mean :
        Data.data -= ma.mean(Data.data, 0)

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    CombineCal(str(sys.argv[1])).execute()

