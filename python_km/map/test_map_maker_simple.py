"""Unit tests for map maker."""

import unittest
import scipy as sp
import numpy.ma as ma

import map_maker_simple as mm

class TestCalcFunctions(unittest.TestCase) :
    
    def test_calc_inds(self) :
        pointing = [4.001, 3.999, -0.999, 8.999]
        centre = 4.
        shape = 10
        inds = mm.calc_inds(pointing, centre, shape)
        self.assertEqual(inds[0], 5)
        self.assertEqual(inds[1], 4)
        self.assertEqual(inds[2], 0)
        self.assertEqual(inds[3], 9)

        shape = 5
        spacing = 2.0
        inds = mm.calc_inds(pointing, centre, shape, spacing)
        self.assertEqual(inds[0], 2)
        self.assertEqual(inds[1], 2)
        self.assertEqual(inds[2], 0)
        self.assertEqual(inds[3], 4)

    def test_cal_bins(self) :
        centre = 5
        shape = 5
        spacing = 2
        bins = mm.calc_bins(centre, shape, spacing, edge='lower')
        self.assertEqual(bins[0], 0)
        self.assertEqual(bins[-1], 8)
        self.assertEqual(bins[2], 4)
        bins = mm.calc_bins(centre, shape, spacing, edge='middle')
        self.assertEqual(bins[0], 1)
        self.assertEqual(bins[-1], 9)
        self.assertEqual(bins[2], 5)
        shape = 10
        spacing = 1
        bins = mm.calc_bins(centre, shape, spacing, edge='upper')
        self.assertEqual(bins[0], 1)
        self.assertEqual(bins[-1], 10)
        self.assertEqual(bins[4], 5)

class TestAddingData(unittest.TestCase) :
    
    def setUp(self) :
        self.shape = (5, 5, 20)
        self.ntime = 15
        self.map = ma.zeros(self.shape, dtyp=float)
        self.counts = sp.zeros(self.shape, dtyp=int)
        self.data = ma.ones((self.ntime, self.shape[-1]))
        self.centre = (90.0,0.0,100.0)
        self.spacing = 1.0
        self.ra = 0.3*sp.arange(self.shape[0]) + self.centre[0]
        self.dec = 0.3*sp.arange(self.shape[1]) + self.centre[1]
        
    def test_adds_data(self) :
        mm.add_data(self.data, self.ra, self.dec, self.map)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1])
        self.assertAlmostEqual(ma.sum(self.counts), self.ntime*self.shape[-1])

    def tearDown(self) :
        del self.data
        del self.map
        del self.counts
        


if __name__ == '__main__' :
    unittest.main()
