"""Unit tests for preprocess module."""

import unittest
import copy

import scipy as sp
import numpy.ma as ma

import preprocess
import core.data_block
import core.fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestHanning(unittest.TestCase) :
    """Since Hanning acctually changes the data, these are only sanity
    tests, far from thurough."""
    
    def setUp(self) :
        Reader = core.fitsGBT.Reader(test_file)
        self.Data = Reader.read(1,1)
    
    def test_hanning_data_changed(self) :
        """Copy the data, see that we did something."""
        self.Data.verify()
        Data_copy = copy.deepcopy(self.Data)
        preprocess.hanning(self.Data)
        self.Data.verify()
        self.assertTrue(not ma.allclose(self.Data.data, Data_copy.data))
        self.assertTrue(sp.allclose(self.Data.field['LST'],
                        Data_copy.field['LST']))

    def tearDown(self) :
        del self.Data


if __name__ == '__main__' :
    unittest.main()

