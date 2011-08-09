"""Unit tests for window stitching modules."""

import unittest
import os

import scipy as sp
import numpy.ma as ma
import matplotlib.pyplot as plt

import kiyopy.custom_exceptions as ce
import stitch_windows_crude as swc
from core import data_block, fitsGBT

test_file = 'testdata/testfile_GBTfits.fits'

class TestStitchWindowsCrude(unittest.TestCase) :

    def setUp(self) :
        Reader = fitsGBT.Reader('testfile_GBTfits.fits', feedback=0)
        self.data_blocks = Reader.read(0,())

    def test_runs(self) :
        swc.stitch(self.data_blocks)

    def test_has_required_keys(self) :
        del self.data_blocks[0].field['CRVAL1']
        self.assertRaises(ce.DataError, swc.stitch, self.data_blocks)

    def test_all_keys_same(self) :
        self.data_blocks[0].field['SCAN'] = 3
        self.assertRaises(ce.DataError, swc.stitch, self.data_blocks)

    def test_fails_no_overlap(self) :
        self.data_blocks[0].field['CRVAL1'] = 800e6
        self.assertRaises(ce.DataError, swc.stitch, self.data_blocks)

    def test_stitches(self) :
        self.data_blocks[0].calc_freq()
        min = self.data_blocks[0].freq[-1]
        # Add a random factor to this one and make sure it gets devided out.
        self.data_blocks[0].data[:,:,:,:] = self.data_blocks[0].freq*1.05
        self.data_blocks[1].calc_freq()
        max = self.data_blocks[1].freq[0]
        self.data_blocks[1].data[:,:,:,:] = self.data_blocks[1].freq
        
        NewData = swc.stitch(self.data_blocks)
        tol = abs(NewData.field['CDELT1']/2)
        data_col = NewData.data[2,2,0,:]
        self.assertTrue(abs(data_col[-1] - min) < tol)
        self.assertTrue(abs(data_col[0] - max) < tol)
        self.assertTrue(sp.allclose(sp.sort(-data_col), -data_col))

        # Make sure the stitched frequencies line up frequency axis.
        NewData.calc_freq()
        freq = NewData.freq
        self.assertTrue(sp.alltrue(abs(data_col - freq) < tol))

    def tearDown(self) :
        del self.data_blocks

class TestModule(unittest.TestCase) :
    
    def test_runs(self) :
        # Just run the module with default params, make sure it works.
        Sticher = swc.Stitch(feedback=0)
        Sticher.execute()
        os.remove('testfile_GBTfits.testout.fits')


if __name__ == '__main__' :
    unittest.main()

