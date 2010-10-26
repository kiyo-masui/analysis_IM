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
        bins = mm.calc_bins(centre, shape, spacing, edge='left')
        self.assertEqual(bins[0], 0)
        self.assertEqual(bins[-1], 8)
        self.assertEqual(bins[2], 4)
        bins = mm.calc_bins(centre, shape, spacing, edge='middle')
        self.assertEqual(bins[0], 1)
        self.assertEqual(bins[-1], 9)
        self.assertEqual(bins[2], 5)
        shape = 10
        spacing = 1
        bins = mm.calc_bins(centre, shape, spacing, edge='right')
        self.assertEqual(bins[0], 1)
        self.assertEqual(bins[-1], 10)
        self.assertEqual(bins[4], 5)

    def test_circle(self) :
        centre = 4.6
        shape = 19
        spacing = 1.3
        bins = mm.calc_bins(centre, shape, spacing, edge='left') + 0.1*spacing
        inds = mm.calc_inds(bins, centre, shape, spacing)
        self.assertTrue(sp.allclose(inds, sp.arange(int(shape))))
        spacing = -0.35
        bins = mm.calc_bins(centre, shape, spacing, edge='left') + 0.1*spacing
        inds = mm.calc_inds(bins, centre, shape, spacing)
        self.assertTrue(sp.allclose(inds, sp.arange(int(shape))))

    def test_calc_inds_neg(self) :
        """Needs to work for negitive spacing."""
        pointing = [2.001, 1.999, 0.001, 3.5]
        centre = 2.
        shape = 4
        spacing = -1
        inds = mm.calc_inds(pointing, centre, shape, spacing)
        self.assertEqual(inds[0], 1)
        self.assertEqual(inds[1], 2)
        self.assertEqual(inds[2], 3)
        self.assertEqual(inds[3], 0)


class TestAddingData(unittest.TestCase) :
    
    def setUp(self) :
        self.shape = (5, 5, 20)
        self.ntime = 15
        self.map = ma.zeros(self.shape, dtype=float)
        self.counts = sp.zeros(self.shape, dtype=int)
        self.data = ma.ones((self.ntime, self.shape[-1]))
        self.centre = (90.0,0.0,100.0)
        self.spacing = 1.0
        self.ra = 0.2*(sp.arange(self.ntime)-self.ntime/2.0) + self.centre[0]
        self.dec = 0.2*(sp.arange(self.ntime)-self.ntime/2.0) + self.centre[1]
        self.ra_inds = mm.calc_inds(self.ra, self.centre[0], self.shape[0],
                                    self.spacing)
        self.dec_inds = mm.calc_inds(self.dec, self.centre[1], self.shape[1],
                                    self.spacing)
        
    def test_adds_data(self) :
        mm.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.counts)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1])
        self.assertAlmostEqual(ma.sum(self.counts), self.ntime*self.shape[-1])

    def test_mask_data(self) :
        self.data[5,8] = ma.masked
        mm.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.counts)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1]-1)
        self.assertAlmostEqual(ma.sum(self.counts), self.ntime*self.shape[-1]-1)

    def test_off_map(self) :
        self.ra_inds[3] = 10
        self.dec_inds[6] = -1
        mm.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.counts)
        self.assertAlmostEqual(ma.sum(self.map),
                               self.ntime*self.shape[-1]-2*self.shape[2])
        self.assertAlmostEqual(ma.sum(self.counts), 
                               self.ntime*self.shape[-1]-2*self.shape[2])

    # test a pointing we know is right

    def test_checks_shapes(self) :
        self.assertRaises(ValueError, mm.add_data_2_map, self.data, self.ra,
                          sp.arange(self.ntime+1), self.map, self.counts)
        self.assertRaises(ValueError, mm.add_data_2_map, self.data, self.ra,
                          self.dec, self.map[:,:,0:-1], self.counts)

    def tearDown(self) :
        del self.data
        del self.map
        del self.counts
        


if __name__ == '__main__' :
    unittest.main()
