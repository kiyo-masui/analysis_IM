"""unit tests for beam.py."""

import unittest
import copy

import scipy as sp


from core import algebra
import beam
import kiyopy.custom_exceptions as ce


class TestGaussianSetUp(unittest.TestCase):

    def setUp(self) :
        self.width = (5 + sp.arange(6, dtype=float))/10
        self.frequencies = sp.arange(6, dtype=float)

    def test_init(self) :
        self.assertRaises(ValueError, beam.GaussianBeam, sp.arange(4),
                          sp.arange(5))

    def test_beam_function(self) :
        self.Beam = beam.GaussianBeam(self.width, self.frequencies)
        # Check a case calculated by hand (FWHM=1 here).
        self.assertAlmostEqual(self.Beam.beam_function(1.0, 5.0), 0.0551589, 6)
        # Check that the delta r dependance is right.
        norm = self.Beam.beam_function(0, 3.0)
        self.assertAlmostEqual(sp.log(self.Beam.beam_function(0.6, 3.0)/norm), 
                               sp.log(self.Beam.beam_function(1.2, 3.0)/norm)/4)
        # Check that it works if I give an array.
        self.Beam.beam_function(sp.arange(5.0), 5.0)
        # Or 2 arrays of the same length.
        self.Beam.beam_function(sp.arange(5), sp.arange(5))
        # But fails if they are different lengths.
        self.assertRaises(ValueError, self.Beam.beam_function, sp.arange(5),
                          sp.arange(6))

    def test_freq_independant(self) :
        self.Beam = beam.GaussianBeam(1.0)
        self.assertAlmostEqual(self.Beam.beam_function(1.0, 2.0), 0.0551589, 6)
        self.assertAlmostEqual(self.Beam.beam_function(1.0, None), 
                               0.0551589, 6)

    def test_extrapolate(self) :
        self.Beam = beam.GaussianBeam(self.width, self.frequencies,
                                      extrapolate=False)
        self.assertRaises(Exception, self.Beam.beam_function, 1.0, 6)
        self.Beam = beam.GaussianBeam(self.width, self.frequencies,
                                      extrapolate=True)
        # Compare FWHM=1 to FWHM = 1/4.  Just check normalization.
        self.assertAlmostEqual(float(self.Beam.beam_function(0, 5.0)), 
                               float(self.Beam.beam_function(0, -2.5))/4**2)

class TestOperatorsMap(unittest.TestCase) :
    
    def setUp(self):
        arr = sp.zeros((10, 20, 15))
        self.map = algebra.make_vect(arr, ('freq', 'ra', 'dec'))
        self.map.set_axis_info('freq', 700.0, 5.0)
        self.map.set_axis_info('ra', 215.5, 0.075*sp.cos(2*sp.pi/180))
        self.map.set_axis_info('dec', 2.0, 0.075)

    def test_left_apply_map_norm(self) :
        self.Beam = beam.GaussianBeam(0.2, None)
        self.map[3, 6, 8] = 1.0
        self.map[9, 13, 6] = 2.0
        self.map[2, 13, 6] = 1.0
        conv_map = self.Beam.apply(self.map)
        self.assertAlmostEqual(sp.sum(self.map.flat_view()),
                               sp.sum(conv_map.flat_view()), 4)

    def test_left_apply_map_vals(self) :
        self.Beam = beam.GaussianBeam(0.2, None)
        self.map[3, 2, 13] = 1.0
        self.map[3, 19, 13] = 1.0
        conv_map = self.Beam.apply(self.map)
        self.assertNotAlmostEqual(conv_map[3, 1, 13], 0, 3)
        # Check Isotropy
        self.assertAlmostEqual(conv_map[3, 1, 13], conv_map[3, 3, 13])
        self.assertAlmostEqual(conv_map[3, 1, 13], conv_map[3, 2, 12])
        self.assertAlmostEqual(conv_map[3, 1, 13], conv_map[3, 2, 14])
        # Check 0 lag amplitude.
        sig = 0.2/2/sp.sqrt(2*sp.log(2))
        self.assertAlmostEqual(conv_map[3, 2, 13],
                               0.075**2/(2*sp.pi*sig**2))

    def test_checks_axes(self) :
        self.map.axes=('freq', 'lat', 'long')
        self.Beam = beam.GaussianBeam(0.2, None)
        self.assertRaises(ce.DataError, self.Beam.apply, self.map)

    def test_wrap(self) :
        self.Beam = beam.GaussianBeam(0.2, None)
        self.map[5, 0, 0] = 1.0
        conv_map = self.Beam.apply(self.map, wrap=True)
        self.assertAlmostEqual(sp.sum(self.map.flat_view()),
                               sp.sum(conv_map.flat_view()), 4)
        self.assertNotAlmostEqual(conv_map[5, 0, 0], 1.0, 3)
        self.assertNotAlmostEqual(conv_map[5, -1, 0], 0, 3)
        self.assertNotAlmostEqual(conv_map[5, -1, -1], 0, 3)
        self.assertAlmostEqual(conv_map[5, -1, 0], conv_map[5, 0, -1])
        self.assertAlmostEqual(conv_map[5, -1, -1], conv_map[5, 1, 1])

class TestOperatorsMatrix(unittest.TestCase) :

    def setUp(self) :
        arr = sp.zeros((2, 20, 15, 20, 15))
        self.op = algebra.make_mat(arr, col_axes = (0, 3, 4), 
                                   row_axes = (0, 1, 2),
                                   axis_names=('freq','ra','dec','ra','dec'))
        self.op.set_axis_info('freq', 700.0, 5.0)
        self.op.set_axis_info('ra', 215.5, 0.075*sp.cos(2*sp.pi/180))
        self.op.set_axis_info('dec', 2.0, 0.075)
        arr = sp.zeros((2, 20, 15))
        self.map = algebra.make_vect(arr, ('freq', 'ra', 'dec'))
        self.map.set_axis_info('freq', 700.0, 5.0)
        self.map.set_axis_info('ra', 215.5, 0.075*sp.cos(2*sp.pi/180))
        self.map.set_axis_info('dec', 2.0, 0.075)

    def test_apply(self) :
        self.op[0,9,8,6,4] = 2.0
        self.op[1,11,7,2,3] = 2.0
        self.op[1,10,8,2,3] = 1.0
        self.map[1, 11, 7] = 2.0
        self.map[1, 10, 8] = 1.0
        self.Beam = beam.GaussianBeam([0.2, 0.3], [600, 900])
        conv_op = self.Beam.apply(self.op)
        conv_map = self.Beam.apply(self.map)
        self.assertAlmostEqual(sp.sum(self.op), sp.sum(conv_op), 3)
        self.assertTrue(sp.allclose(conv_op[:, :, :, 2, 3], conv_map))

    def test_apply_right(self) :
        self.op[0,9,8,6,4] = 2.0
        self.op[1,11,7,11,8] = 2.0
        self.op[1,11,7,8,6] = 1.0
        self.map[1, 11, 8] = 2.0
        self.map[1, 8, 6] = 1.0
        self.Beam = beam.GaussianBeam([0.2, 0.3], [600, 900])
        conv_op = self.Beam.apply(self.op, right_apply=True)
        conv_map = self.Beam.apply(self.map)
        self.assertAlmostEqual(sp.sum(self.op), sp.sum(conv_op), 3)
        self.assertTrue(sp.allclose(conv_op[:, 11, 7, :, :], conv_map))
    
    def test_checks_axes(self) :
        self.op.axes = ('freq','ra','dec','mode1','mode2')
        self.Beam = beam.GaussianBeam([0.2, 0.3], [600, 900])
        self.assertRaises(ce.DataError, self.Beam.apply, self.op,
                         right_apply=True)
        conv_op = self.Beam.apply(self.op)
        self.assertEqual(conv_op.axes, ('freq','ra','dec','mode1','mode2'))

if __name__ == '__main__' :
    unittest.main()



