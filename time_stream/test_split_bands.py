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
    
    def setUp(self):
        Reader = core.fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,1)
        self.Data.verify()


if __name__ == '__main__' :
    unittest.main()

