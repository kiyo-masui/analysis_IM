"""Unit tests for cal_scale.py"""

import unittest

import numpy.ma as ma

from time_stream import cal_scale
from core import data_block, fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestCalScale(unittest.TestCase) :
    
    def setUp(self) :
        self.Reader = fitsGBT.Reader(test_file)
        self.Data = self.Reader.read(0,1)

    def test_scale(self) :
        cal_scale.scale_by_cal(self.Data)
        data = self.Data.data
        self.assertTrue(ma.allclose(ma.median(data[:,0,0,:] -
                                              data[:,0,1,:], 0), 1.0))
        self.assertTrue(ma.allclose(ma.median(data[:,3,0,:] -
                                              data[:,3,1,:], 0), 1.0))
        #self.assertTrue(ma.logical
    
    def tearDown(self) :
        del self.Reader
        del self.Data

if __name__ == '__main__' :
    unittest.main()

