"""Unit tests for hanning a rebin_freq modules."""

import unittest
import copy

import scipy as sp
import numpy.ma as ma

import sys

import split_bands
import core.data_block
import core.fitsGBT

test_file = 'testdata/testfile_guppi_rebinned.fits'

class TestFunctions(unittest.TestCase) :
    
    def setUp(self):
        Reader = core.fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(0, 0)
        self.Data.verify()

    def test_works(self):
        Blocks = split_bands.split(self.Data, 4, 30, 3)
        self.assertEqual(4, len(Blocks))
        self.Data.calc_freq()
        all_freq = self.Data.freq
        for ii in xrange(4):
            this_data = self.Data.data[...,3 + ii * 30:3 + (ii + 1) * 30]
            self.assertTrue(sp.allclose(Blocks[ii].data.filled(-1), 
                                        this_data.filled(-1)))
            this_freq = all_freq[3 + ii * 30:3 + (ii + 1) * 30]
            Blocks[ii].calc_freq()
            self.assertTrue(sp.allclose(this_freq, Blocks[ii].freq))

    def test_raise(self):
        self.assertRaises(ValueError, split_bands.split, self.Data, 2, 60, 23)
    
        


        


if __name__ == '__main__' :
    unittest.main()

