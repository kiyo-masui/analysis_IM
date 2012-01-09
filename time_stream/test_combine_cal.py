"""Unit tests for combine_cal.py"""

import unittest

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
from time_stream import combine_cal, rebin_freq, hanning
from core import data_block, fitsGBT

test_file = 'testdata/testfile_GBTfits.fits'

class TestCombineCal(unittest.TestCase) :
    
    def setUp(self) :
        self.Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = self.Reader.read(1,0)

    def test_cal_checking(self) :
        self.Data.field['CAL'][1] = 'T'
        self.assertRaises(ce.DataError, combine_cal.combine, self.Data)

    def test_subtracts_baseline(self) :
        rebin_freq.rebin(self.Data, 1.0)
        combine_cal.combine(self.Data, sub_mean=True, average_cals=False)
        data = self.Data.data
        self.assertTrue(ma.allclose(ma.mean(data, 0), 0.))

    def test_dims_cal_ave(self) :
        dims = self.Data.dims
        combine_cal.combine(self.Data, sub_mean=False, average_cals=True)
        dims = dims[:2] + (1,) + dims[3:]
        self.assertEqual(self.Data.dims, dims)
        self.Data.verify()


    def tearDown(self) :
        del self.Reader
        del self.Data

if __name__ == '__main__' :
    unittest.main()

