#!/usr/bin/python
"""Module that put data in units of cal temperture and subtracts median."""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single


class CalScale(base_single.BaseSingle) :
    """Pipeline module that performs scales data by the cal on the data."""

    prefix = 'cs_'
    params_init = {'subtract_time_median' : False}

    def action(self, Data):
        scale_by_cal(Data, self.params['subtract_time_median'])
        Data.add_history('Converted to units of noise cal temperture.')
        return Data


def scale_by_cal(Data, sub_med=False) :
    """Puts all data in units of the cal temperature."""
    
    # Here we check the polarizations and cal indicies
    xx_ind = 0
    yy_ind = 3
    xy_inds = [1,2]
    if (Data.field['CRVAL4'][xx_ind] != -5 or
        Data.field['CRVAL4'][yy_ind] != -6 or
        Data.field['CRVAL4'][xy_inds[0]] != -7 or
        Data.field['CRVAL4'][xy_inds[1]] != -8) :
            raise ce.DataError('Polarization types not as expected.')
    on_ind = 0
    off_ind = 1
    if (Data.field['CAL'][on_ind] != 'T' or
        Data.field['CAL'][off_ind] != 'F') :
            raise ce.DataError('Cal states not in expected order.')
    # Find the cal medians and scale by them.
    cal_med_xx = ma.median(Data.data[:,xx_ind,on_ind,:]
                           - Data.data[:,xx_ind,off_ind,:], 0)
    Data.data[:,xx_ind,:,:] /= cal_med_xx
    cal_med_yy = ma.median(Data.data[:,yy_ind,on_ind,:]
                           - Data.data[:,yy_ind,off_ind,:], 0)
    Data.data[:,yy_ind,:,:] /= cal_med_yy
    Data.data[:,xy_inds,:,:] /= ma.sqrt(cal_med_yy*cal_med_xx)

    # Subtract the time median if desired.
    if sub_med :
        Data.data -= ma.median(Data.data, 0)

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    CalScale(str(sys.argv[1])).execute()

