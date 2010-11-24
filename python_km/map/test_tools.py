"""Unit tests for map tools."""

import unittest
import scipy as sp
import numpy.ma as ma

from core import fitsGBT
import tools

test_file = "testfile_GBTfits.fits"

class TestCalcFunctions(unittest.TestCase) :
    
    def test_calc_inds(self) :
        pointing = [4.001, 3.999, -0.999, 8.999]
        centre = 4.
        shape = 10
        inds = tools.calc_inds(pointing, centre, shape)
        self.assertEqual(inds[0], 5)
        self.assertEqual(inds[1], 4)
        self.assertEqual(inds[2], 0)
        self.assertEqual(inds[3], 9)

        shape = 5
        spacing = 2.0
        inds = tools.calc_inds(pointing, centre, shape, spacing)
        self.assertEqual(inds[0], 2)
        self.assertEqual(inds[1], 2)
        self.assertEqual(inds[2], 0)
        self.assertEqual(inds[3], 4)

    def test_calc_bins(self) :
        centre = 5
        shape = 5
        spacing = 2
        bins = tools.calc_bins(centre, shape, spacing, edge='left')
        self.assertEqual(bins[0], 0)
        self.assertEqual(bins[-1], 8)
        self.assertEqual(bins[2], 4)
        bins = tools.calc_bins(centre, shape, spacing, edge='middle')
        self.assertEqual(bins[0], 1)
        self.assertEqual(bins[-1], 9)
        self.assertEqual(bins[2], 5)
        shape = 10
        spacing = 1
        bins = tools.calc_bins(centre, shape, spacing, edge='right')
        self.assertEqual(bins[0], 1)
        self.assertEqual(bins[-1], 10)
        self.assertEqual(bins[4], 5)

    def test_circle(self) :
        centre = 4.6
        shape = 19
        spacing = 1.3
        bins = tools.calc_bins(centre, shape, spacing, edge='left') + 0.1*spacing
        inds = tools.calc_inds(bins, centre, shape, spacing)
        self.assertTrue(sp.allclose(inds, sp.arange(int(shape))))
        spacing = -0.35
        bins = tools.calc_bins(centre, shape, spacing, edge='left') + 0.1*spacing
        inds = tools.calc_inds(bins, centre, shape, spacing)
        self.assertTrue(sp.allclose(inds, sp.arange(int(shape))))

    def test_calc_inds_neg(self) :
        """Needs to work for negitive spacing."""
        pointing = [2.001, 1.999, 0.001, 3.5]
        centre = 2.
        shape = 4
        spacing = -1
        inds = tools.calc_inds(pointing, centre, shape, spacing)
        self.assertEqual(inds[0], 1)
        self.assertEqual(inds[1], 2)
        self.assertEqual(inds[2], 3)
        self.assertEqual(inds[3], 0)


class TestMapSetUp(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(0,1)

    def test_runs_no_errors(self) :
        Map = tools.set_up_map(self.Data, (0.0, 5.2), (10,4), [0.5,0.1])

    def test_checks_tuple_args(self) :
        self.assertRaises(ValueError, tools.set_up_map, self.Data, 5, (10,10),
                          (0.5,0.5))
        self.assertRaises(ValueError, tools.set_up_map, self.Data, (5, 8), 
                          (10, 10, 11), (0.5,0.5))

    def test_map_dims(self) :
        shape = (20,15)
        Map = tools.set_up_map(self.Data, (0.0, 0.0), shape, [0.5,0.5])
        nf = self.Data.dims[-1]
        self.assertEqual(Map.dims, shape + (nf,))

    def test_gets_axes_right(self) :
        Map = tools.set_up_map(self.Data, (0.0, 5.2), (45,23), (0.4,0.5))
        Map.calc_axes()
        self.Data.calc_freq()
        self.assertTrue(sp.allclose(Map.freq, self.Data.freq))
        self.assertTrue(sp.allclose(Map.long, tools.calc_bins(0.0, 45, 0.4,
                                                            'middle')))
        self.assertTrue(sp.allclose(Map.lat, tools.calc_bins(5.2, 23, 0.5,
                                                            'middle')))

    def tearDown(self) :
        del self.Data

# TODO: Test that the bins calculated here agree with the ones from
# data_map.calc_axes().


if __name__ == '__main__' :
    unittest.main()
