"""Unit tests for cal_scale.py"""

import unittest

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
from time_stream import cal_scale, rebin_freq, hanning
from core import data_block, fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestCalScale(unittest.TestCase) :
    
    def setUp(self) :
        self.Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = self.Reader.read(1,0)

    def test_scale(self) :
        hanning.hanning_smooth(self.Data)
        rebin_freq.rebin(self.Data, 2.)
        cal_scale.scale_by_cal(self.Data)
        data = self.Data.data

        self.assertTrue(ma.allclose(ma.median(data[:,0,0,:] -
                                              data[:,0,1,:], 0), 1.0))
        self.assertTrue(ma.allclose(ma.median(data[:,3,0,:] -
                                              data[:,3,1,:], 0), 1.0))
        #self.assertTrue(ma.allclose(ma.median(
        #                   (data[:,1,0,:] - data[:,1,1,:])**2
        #                   + (data[:,2,0,:] - data[:,2,1,:])**2, 0), 1.0, 
        #                   rtol=0.05))

    def test_pol_checking(self) :
        self.Data.field['CRVAL4'][0] = 1
        self.assertRaises(ce.DataError, cal_scale.scale_by_cal, self.Data)

    def test_cal_checking(self) :
        self.Data.field['CAL'][1] = 'T'
        self.assertRaises(ce.DataError, cal_scale.scale_by_cal, self.Data)

    def test_subtracts_baseline(self) :
        rebin_freq.rebin(self.Data, 1.0)
        cal_scale.scale_by_cal(self.Data, sub_med=True)
        data = self.Data.data
        self.assertTrue(ma.allclose(ma.median(data, 0), 0.))
        # The following fails if you get rid of the rebin line, but in the 7th
        # digit.  Numpy must have only single precision somewhere.
        #self.assertAlmostEqual(ma.median(data[:,0,0,753]), 0.)
        self.assertAlmostEqual(ma.median(data), 0.)

    def test_fave_scale(self) :
        hanning.hanning_smooth(self.Data)
        rebin_freq.rebin(self.Data, 2.)
        cal_scale.scale_by_cal(self.Data, False, True, False)
        data = self.Data.data

        self.assertTrue(ma.allclose(ma.mean(data[:,0,0,:] -
                                              data[:,0,1,:], -1), 1.0))
        self.assertTrue(ma.allclose(ma.mean(data[:,3,0,:] -
                                              data[:,3,1,:], -1), 1.0))

    def test_both_scale(self) :
        data = self.Data.data
        # linear so median = mean
        gt = 0.5*sp.arange(10) + 3.
        gf = sp.cos(0.01*sp.arange(2048)) + 4
        data[:,:,1,:] = 0
        # explicit phase closure: 2 = sqrt(4*1)
        data[:,0,0,:] = gt[:,sp.newaxis]*gf
        data[:,3,0,:] = gt[:,sp.newaxis]*gf*4
        data[:,1,0,:] = gt[:,sp.newaxis]*gf*2
        data[:,2,0,:] = gt[:,sp.newaxis]*gf*2
        cal_scale.scale_by_cal(self.Data, True, True, False)
        self.assertTrue(ma.allclose(data[:,:,0,:], 1.))
        self.assertTrue(ma.allclose(data[:,:,1,:], 0.))

    def tearDown(self) :
        del self.Reader
        del self.Data

if __name__ == '__main__' :
    unittest.main()

