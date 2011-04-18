"""Unit tests for reflag.py"""

import unittest
import os
import glob

import scipy as sp
import numpy.ma as ma

from core import fitsGBT
import reflag
import rotate_pol
import combine_cal

test_file = 'testfile_GBTfits.fits'

class TestFlag(unittest.TestCase) :
        
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,0)
        rotate_pol.rotate(self.Data, (1,))
        self.DataSub = Reader.read(1,0)
        rotate_pol.rotate(self.DataSub, (1,))

    def test_runs(self) :
        reflag.flag(self.Data, self.DataSub, thres=3)

    def test_flags(self) :
        self.Data.data[...] = 10
        self.Data.data[4, 0, 0, 235] = 100
        nf = self.Data.dims[-1]
        self.DataSub.data[0::2, :, :, :] = sp.arange(1, nf+1)
        self.DataSub.data[1::2, :, :, :] = 0
        self.DataSub.data[3, 0, 0, 458] = 2*458
        self.DataSub.data[3, 0, 0, 986] = 2*986
        self.DataSub.data[8, 0, 1, 734] = 2*734
        reflag.flag(self.Data, self.DataSub, thres=2.0)
        self.assertTrue(self.Data.data[3, 0, 0, 458] is ma.masked)
        self.assertTrue(self.Data.data[3, 0, 0, 986] is ma.masked)
        self.assertTrue(self.Data.data[8, 0, 1, 734] is ma.masked)
        self.assertEqual(ma.count_masked(self.Data.data), 3)
        self.assertEqual(ma.count_masked(self.DataSub.data), 3)

class TestModule(unittest.TestCase) :

    def test_module(self) :
        params = {'sf_thres' : 3.0,
                  'sf_subtracted_input_root' : './',
                  'sf_subtracted_output_root' : './subtracted_'
                 }
        reflag.ReFlag(params, feedback=0).execute()
        files = glob.glob('*testout*')
        self.assertTrue(len(files) > 1)

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)


if __name__ == '__main__' :
    unittest.main()

