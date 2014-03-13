"""Unit tests for hanning a rebin_freq modules."""

import unittest
import copy
import glob
import os

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


class TestModule(unittest.TestCase) :

    def test_module(self) :
        params = {'sb_n_bands' : 2,
                  'sb_n_bins_band' : 10,
                  'sb_offset' : 3,
                  'sb_input_root' : './testdata/',
                  'sb_file_middles' : ('testfile_guppi',),
                  'sb_input_end' : '_rebinned.fits',
                  'sb_output_root' : './testout_',
                  'sb_output_end' : '_split.fits'
                 }
        split_bands.SplitBands(params, feedback=0).execute()
        files = glob.glob('*testout*')
        self.assertTrue(len(files) > 1)

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)

        


        


if __name__ == '__main__' :
    unittest.main()

