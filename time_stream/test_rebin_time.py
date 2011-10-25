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

    def test_rebin_runs(self):
        rebin_time.rebin(self.Data, 2)
        self.Data.verify()

    def test_time_axis(self):
        self.Data.calc_time()
        t1 = self.Data.time
        rebin_time.rebin(self.Data, 4)
        self.Data.calc_time()
        t2 = self.Data.time
        self.assertEqual(len(t2), 2)
        self.assertAlmostEqual(t2[0], sp.mean(t1[:4]))
        self.assertAlmostEqual(t2[1], sp.mean(t1[4:8]))

    def test_fields(self):
        lst1 = self.Data.field['LST']
        exp1 = self.Data.field['EXPOSURE']
        rebin_time.rebin(self.Data, 4)
        lst2 = self.Data.field['LST']
        exp2 = self.Data.field['EXPOSURE']
        self.assertEqual(len(lst2), 2)
        self.assertAlmostEqual(lst2[0], sp.mean(lst1[:4]))
        self.assertAlmostEqual(lst2[1], sp.mean(lst1[4:8]))
        self.assertEqual(exp2.shape[1], exp1.shape[1])
        self.assertEqual(exp2.shape[0], 2)
        self.assertTrue(sp.allclose(exp2[0,:], sp.mean(exp1[:4,:], 0)))
        self.assertTrue(sp.allclose(exp2[1,:], sp.mean(exp1[4:8,:], 0)))

    def test_data(self):
        nf = self.Data.dims[-1]
        # Rebin by 3's.
        data = sp.arange(10, dtype=float) * 3.0
        data = data[...,None] * sp.arange(4)
        data = data[...,None] * sp.arange(2)
        data = data[...,None] * sp.arange(nf)
        new_data = sp.arange(3, dtype=float) * 9.0 + 3.0
        new_data = new_data[...,None] * sp.arange(4)
        new_data = new_data[...,None] * sp.arange(2)
        new_data = new_data[...,None] * sp.arange(nf)
        self.Data.data[...] = data
        rebin_time.rebin(self.Data, 3)
        self.assertTrue(sp.allclose(self.Data.data, new_data))

    def test_data_masks(self):
        nf = self.Data.dims[-1]
        data = ma.arange(10, dtype=float) * 3.0
        data = data[...,None] * sp.arange(4)
        data = data[...,None] * sp.arange(2)
        data = data[...,None] * sp.arange(nf)
        data[7,1,1,16] = ma.masked
        data[3:8,2,1,15] = ma.masked
        self.Data.data[...] = data
        rebin_time.rebin(self.Data, 3)
        self.assertAlmostEqual(self.Data.data[2,1,1,16],
                               ma.mean(data[6:9,1,1,16]))
        self.assertTrue(self.Data.data[1,2,1,15] is ma.masked)

    def tearDown(self):
        del self.Data


if __name__ == '__main__' :
    unittest.main()

