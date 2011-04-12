"""Unit tests for map tools."""

import unittest
import scipy as sp
import numpy.ma as ma

from core import fitsGBT, fits_map
import tools

test_file = "testfile_GBTfits.fits"
test_map = 'testfile_map.fits'
# Known test map parameters.
tm_centre = (325.7, 0.0)
tm_shape = (22, 10)
tm_spacing = (-0.2, 0.2)


class TestCalcFunctions(unittest.TestCase) :
    
    def test_calc_inds(self) :
        pointing = [4.001, 3.999, -0.999, 8.999]
        centre = 4.
        shape = 10
        inds = tools.calc_inds(pointing, centre, shape)
        self.assertEqual(inds[0], 5)
        self.assertEqual(inds[1], 5)
        self.assertEqual(inds[2], 0)
        self.assertEqual(inds[3], 10)

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
        self.assertEqual(bins[0], 0.5)
        self.assertEqual(bins[-1], 9.5)
        self.assertEqual(bins[4], 4.5)

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
        pointing = [2.001, 1.999, 0.001, 3.49]
        centre = 2.
        shape = 4
        spacing = -1
        inds = tools.calc_inds(pointing, centre, shape, spacing)
        self.assertEqual(inds[0], 2)
        self.assertEqual(inds[1], 2)
        self.assertEqual(inds[2], 4)
        self.assertEqual(inds[3], 1)


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


class TestFileVar(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Blocks = Reader.read((),1)

    def test_known(self) :
        dims = self.Blocks[0].dims
        data = ma.ones(dims)*sp.arange(1, dims[-1]+1)
        data[0:-1:2,0,0,:] = 2*sp.arange(1, dims[-1]+1)
        for Data in self.Blocks :
            Data.data = data
        var = tools.calc_time_var_file(self.Blocks, 0, 0)
        self.assertTrue(sp.allclose(var.filled(), 
                        sp.arange(1, dims[-1]+1)**2/4.0))

    def test_known2(self) :
        dims = self.Blocks[0].dims
        data = ma.ones(dims)*sp.arange(1, dims[-1]+1)
        ii = 0.0
        for Data in self.Blocks :
            ii += 1.0
            Data.data = data*ii
        var = tools.calc_time_var_file(self.Blocks, 0, 0)
        self.assertTrue(sp.allclose(var.filled(), 
                        sp.arange(1, dims[-1]+1)**2/4.0))


    def test_known_with_masked(self) :
        dims = self.Blocks[0].dims
        data = ma.ones(dims)*sp.arange(1, dims[-1]+1)
        data[0:-1:2,0,0,:] = 2*sp.arange(1, dims[-1]+1)
        for Data in self.Blocks :
            Data.data = data
        self.Blocks[0].data[2,0,0,43] = ma.masked
        self.Blocks[0].data[7,0,0,43] = ma.masked
        self.Blocks[1].data[4,0,0,43] = ma.masked
        self.Blocks[1].data[3,0,0,43] = ma.masked
        self.Blocks[0].data[:,0,0,103] = ma.masked
        self.Blocks[1].data[:,0,0,103] = ma.masked
        self.Blocks[0].data[:,0,0,554] = ma.masked
        self.Blocks[1].data[1:,0,0,554] = ma.masked
        var = tools.calc_time_var_file(self.Blocks, 0, 0)
        expected = ma.arange(1, dims[-1]+1)**2/4.0
        expected[103] = ma.masked
        expected[554] = ma.masked
        self.assertTrue(sp.allclose(var.filled(-1), expected.filled(-1)))


    def tearDown(self) :
        del self.Blocks


class TestMapPars(unittest.TestCase) :
    
    def setUp(self) :
        self.Map = fits_map.read(test_map, feedback=0)

    def test_test_map(self) :
        centre, shape, spacing = tools.get_map_params(self.Map)
        self.assertAlmostEqual(shape[0], tm_shape[0])
        self.assertAlmostEqual(shape[1], tm_shape[1])
        self.assertAlmostEqual(centre[0], tm_centre[0])
        self.assertAlmostEqual(centre[1], tm_centre[1])
        self.assertAlmostEqual(spacing[0], tm_spacing[0])
        self.assertAlmostEqual(spacing[1], tm_spacing[1])
    
    def tearDown(self) :
        del self.Map


# TODO: Test that the bins calculated here agree with the ones from
# data_map.calc_axes().


if __name__ == '__main__' :
    unittest.main()
