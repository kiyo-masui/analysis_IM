"""Unit tests for utils module."""

import unittest

import scipy as sp

import utils

class TestElAz2RaDec_LST(unittest.TestCase) :
    
    def test_known_cases(self) :
        el = [90, 21, 60]
        az = [12, 180, 90]
        lst = [0, 26.*86400./360., 0]
        lat = [13, 0, 0]
        
        ra, dec = utils.elaz2radec_lst(el, az, lst, lat)

        self.assertAlmostEqual(ra[0], 0)
        self.assertAlmostEqual(dec[0], lat[0])
        
        self.assertAlmostEqual(ra[1], 26, 4)
        self.assertAlmostEqual(dec[1], -90+el[1])
        
        self.assertAlmostEqual(ra[2], 360-30)
        self.assertAlmostEqual(dec[2], 0)

class TestElAz2RaDec(unittest.TestCase) :
    
    def test_works(self) :
        UT = '2000-01-01T12:00:00.00'
        el = 90.
        az = 50.5
        ra, dec = utils.elaz2radecGBT(el, az, UT)
        self.assertAlmostEqual(dec, 38.43312, 1) #GBT Latitude
        el = 38.43312
        az = 0.
        ra, dec = utils.elaz2radecGBT(el, az, UT)
        self.assertAlmostEqual(dec, 90, 1)
        el = 90 - 38.43312
        az = 180.
        ra, dec = utils.elaz2radecGBT(el, az, UT)
        self.assertAlmostEqual(dec, 0, 1)

class TestMakeMapGrid(unittest.TestCase) :
    
    def test_basic(self) :
        centre = [80., 20.]
        shape = (11, 11)
        spacing = 1.
        grid_ra, grid_dec = utils.mk_map_grid(centre, shape, spacing)
        self.assertEqual(sp.shape(grid_ra), shape)
        self.assertEqual(sp.shape(grid_dec), shape)
        self.assertTrue(sp.allclose(grid_ra[:,5], centre[0]))
        self.assertTrue(sp.allclose(grid_dec[5,:], centre[1]))
        self.assertAlmostEqual(sp.amax(grid_dec), 25)
        self.assertAlmostEqual(sp.amin(grid_dec), 15)
        self.assertTrue(sp.allclose(grid_ra[:,7], centre[0] + 
                        2/sp.cos(centre[1]*sp.pi/180.)))

class TestPolInt2String(unittest.TestCase):

    def test_all_functionality(self):
        self.assertEqual(utils.polint2str(1), 'I')
        self.assertEqual(utils.polint2str(2), 'Q')
        self.assertEqual(utils.polint2str(3), 'U')
        self.assertEqual(utils.polint2str(4), 'V')
        self.assertEqual(utils.polint2str(-5), 'XX')
        self.assertEqual(utils.polint2str(-6), 'YY')
        self.assertEqual(utils.polint2str(-7), 'XY')
        self.assertEqual(utils.polint2str(-8), 'YX')
        self.assertEqual(utils.polint2str(-1), 'RR')
        self.assertEqual(utils.polint2str(-2), 'LL')
        self.assertEqual(utils.polint2str(-3), 'RL')
        self.assertEqual(utils.polint2str(-4), 'LR')
        self.assertRaises(ValueError, utils.polint2str, 7)

if __name__ == '__main__' :
    unittest.main()
