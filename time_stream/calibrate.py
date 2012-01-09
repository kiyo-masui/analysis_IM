#!/usr/bin/python
"""Puts time stream data in units of temperature.
"""

import numpy.ma as ma
import scipy as sp
from scipy import interpolate

import kiyopy.custom_exceptions as ce
import base_single
from core import fitsGBT

class Calibrate(base_single.BaseSingle) :
    """Pipeline module converts data from units of cal temperture to Kelvins.

    This module reads the calibrator temperture from a fits file (as a function
    of polarization and frequency) and multiplies it by the time stream data.
    If the time stream data was in units of calibrator temperture, it will end
    up in units of acctual temperature (K).
    """
    
    prefix = 'cl_'
    params_init = {
                   #TODO: Have a real default (like a test file).
                   'cal_temperature_files' : ('some_file_name.fits',)
                   }

    # Add extra stuff to the constructor.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
        base_single.BaseSingle.__init__(self, parameter_file_or_dict, 
                                        feedback)
        # Read in the calibration file.
        file_names = self.params['cal_temperature_files']
        if len(file_names) > 1 :
            raise NotImplementedError('Can use a single calibration file'
                                      ' at a time.')
        Reader = fitsGBT.Reader(file_names[0], feedback=self.feedback)
        self.CalData = Reader.read([],[])
        if type(self.CalData) is tuple :
            raise ce.DataError('Expected calibration file to have only a'
                               ' single scan and IF.')

    def action(self, Data) :
        """Calls multiply_by_cal.
        
        Reads in cal data file on first iteration and then applys the noise cal
        temperature to the data.
        """
        
        # Scale data by noise cal termperature.
        multiply_by_cal(Data, self.CalData)
        Data.add_history("Put into units of Kelvin.")
        return Data


def multiply_by_cal(Data, CalData) :
    """Function scales data by the noise cal temperature.
    """

    # For now we just assume that the cal and polarizations are arranged in a
    # certain way and then check to make sure we are right.
    calibrate_to_I = False
    if tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6) :
        xx_ind = 0
        yy_ind = 3
        xy_inds = [1,2]
    elif tuple(Data.field['CRVAL4']) == (1, 2, 3, 4) :
        # This is a hack.  Completly temporairy.
        calibrate_to_I = True
    else :
        raise ce.DataError('Polarization types not as expected in data.')

    cal_xx_ind = 0
    cal_yy_ind = 1
    if (CalData.field['CRVAL4'][cal_xx_ind] != -5 or
        CalData.field['CRVAL4'][cal_yy_ind] != -6) :
            raise ce.DataError('Polarization types not as expected in cal.')

    # Cal should only have 1 time, 1 cal state and 2 polarizations.
    if CalData.dims[:3] != (1,2,1) :
        raise ce.DataError('Cal temperature data has wrong dimensions.')

    # Cal state should be special state 'R'.
    if CalData.field['CAL'][0] != 'R' :
        raise ce.DataError("Cal state in cal temperture data should be "
                           "'R'.")

    # Bring the Cal data to the same frequencies as the other data.
    Data.calc_freq()
    CalData.calc_freq()
    if sp.allclose(Data.freq, CalData.freq) :
        cdata = CalData.data
    elif abs(Data.field['CDELT1']) <= abs(CalData.field['CDELT1']) :
        calfunc = interpolate.interp1d(CalData.freq, CalData.data, 
                                  fill_value=sp.nan, bounds_error=False)
        cdata = ma.array(calfunc(Data.freq))
        cdata[sp.logical_not(sp.isfinite(cdata))] = ma.masked
    else :
        nf = len(Data.freq)
        width = abs(Data.field['CDELT1'])
        cdata = ma.empty((1,2,1,nf))
        for find in range(nf) :
            f = Data.freq[find]
            inds, = sp.where(sp.logical_and(CalData.freq >= f - width/2.0,
                                           CalData.freq < f + width/2.0))
            cdata[:,:,:,find] = ma.mean(CalData.data[:,:,:,inds], 3)
    
    if calibrate_to_I :
        Data.data *= (cdata[0,cal_xx_ind,0,:] + cdata[0,cal_yy_ind,0,:])/2.0
    else :
        # Loop over times and cal and scale each polarization appropriately.
        for tind in range(Data.dims[0]) :
            for cind in range(Data.dims[2]) :
                Data.data[tind,xx_ind,cind,:] *= cdata[0,cal_xx_ind,0,:]
                Data.data[tind,yy_ind,cind,:] *= cdata[0,cal_yy_ind,0,:]
                Data.data[tind,xy_inds,cind,:] *= ma.sqrt(
                     cdata[0,cal_yy_ind,0,:] * cdata[0,cal_xx_ind,0,:])


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Calibrate(str(sys.argv[1])).execute()
