"""Unit tests for hanning a rebin_freq modules."""

import unittest
import copy

import scipy as sp
import numpy.ma as ma

import sys

import rebin_time
import core.data_block
import core.fitsGBT

test_file = 'testdata/testfile_GBTfits.fits'

class TestFunctions(unittest.TestCase) :
    """Since these operations actually changes the data, these are only sanity
    tests, far from thorough."""
    
    def setUp(self) :
        Reader = core.fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,1)
        self.Data.verify()
        self.Data_copy = copy.deepcopy(self.Data)

    def test_rebin_runs(self) :
        rebin_time.rebin(self.Data, 1.0)
        self.Data.verify()

    def tearDown(self) :
        del self.Data
        del self.Data_copy


if __name__ == '__main__' :
    unittest.main()

