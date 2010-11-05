"""Unit tests for calibrate module.
"""

import unittest
import copy

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
from core import fitsGBT, data_block
import calibrate, rebin_freq

test_filename = 'testfile_GBTfits.fits'

class TestMultipyFunction(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_filename, 0)
        self.Data = Reader.read(0,0)
        self.Data.calc_freq()
        calT = self.Data.freq[sp.newaxis,sp.newaxis,sp.newaxis,:]
        calT = sp.concatenate((calT, 5*calT), 1)
        self.CalData = data_block.DataBlock(calT)
        for field in ('CRVAL1', 'CDELT1', 'CRPIX1'):
            self.CalData.set_field(field, self.Data.field[field], 
                self.Data.field_axes[field], self.Data.field_formats[field])
        self.CalData.set_field('CAL', ['R'], ('cal',), '1A')
        self.CalData.set_field('CRVAL4', [-5, -6], ('pol',), '1I')
        self.CalData.verify()

    def test_runs_no_errors(self) :
        calibrate.multiply_by_cal(self.Data, self.CalData)

    
    def test_pol_checking(self) :
        self.Data.field['CRVAL4'][0] = -8
        self.assertRaises(ce.DataError, calibrate.multiply_by_cal, self.Data, 
                          self.CalData)

    def test_cal_checking_cal(self) :
        self.CalData.field['CAL'][0] = 'T'
        self.assertRaises(ce.DataError, calibrate.multiply_by_cal, self.Data, 
                          self.CalData)
    
    def test_pol_checking_cal(self) :
        self.CalData.field['CRVAL4'][1] = -7
        self.assertRaises(ce.DataError, calibrate.multiply_by_cal, self.Data, 
                          self.CalData)
    
    def test_checks_cal_dims(self) :
        self.CalData.dims = (2,2,1,467)
        self.assertRaises(ce.DataError, calibrate.multiply_by_cal, self.Data, 
                          self.CalData)

    def test_scales_right_not_rebined(self) :
        DataCopy = copy.deepcopy(self.Data)
        calibrate.multiply_by_cal(self.Data, self.CalData)
        gains = sp.array([1,sp.sqrt(5),sp.sqrt(5),5])
        gains.shape = (1,4,1,1)
        self.assertTrue(ma.allclose(DataCopy.data*gains*DataCopy.freq, 
                                    self.Data.data))

    def test_scales_rebined(self) :
        self.Data.set_data(ma.ones((10,4,2,600)))
        self.Data.verify()
        rebin_freq.rebin(self.Data, 1.0)
        DataCopy = copy.deepcopy(self.Data)
        DataCopy.calc_freq()
        calibrate.multiply_by_cal(self.Data, self.CalData)
        gains = sp.array([1,sp.sqrt(5),sp.sqrt(5),5])
        gains.shape = (1,4,1,1)
        expected = (DataCopy.data*gains*DataCopy.freq)
        # don't expect this test to be perfect, 4 figures good enough.  Edge
        # bins adversly affected by missing frequencies in CalData.
        self.assertTrue(ma.allclose(expected[:,:,:1:-1], 
                                    self.Data.data[:,:,:1:-1], 4))

    def test_rebinned_cal(self) :
        rebin_freq.rebin(self.CalData, 1.0)
        DataCopy = copy.deepcopy(self.Data)
        DataCopy.calc_freq()
        calibrate.multiply_by_cal(self.Data, self.CalData)
        gains = sp.array([1,sp.sqrt(5),sp.sqrt(5),5])
        gains.shape = (1,4,1,1)
        expected = (DataCopy.data*gains*DataCopy.freq)
        self.assertTrue(ma.count_masked(self.Data.data <
                           sp.size(self.Data.data)*50.0/2000.0))
        expected.mask = sp.array(self.Data.data.mask)
        self.assertTrue(ma.allclose(expected, self.Data.data, 4))


    def tearDown(self) :
        del self.CalData
        del self.Data


if __name__ == '__main__' :
    unittest.main()
