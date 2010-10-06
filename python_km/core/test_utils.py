"""Unit tests for utils module."""

import unittest

import utils

class TestElAz2RaDec(unittest.TestCase) :
    
    def test_known_cases(self) :
        el = [90, 21, 60]
        az = [12, 180, 90]
        lst = [0, 26.*86400./360., 0]
        lat = [13, 0, 0]
        
        ra, dec = utils.elaz2radec(el, az, lst, lat)

        self.assertAlmostEqual(ra[0], 0)
        self.assertAlmostEqual(dec[0], lat[0])
        
        self.assertAlmostEqual(ra[1], 26, 4)
        self.assertAlmostEqual(dec[1], -90+el[1])
        
        self.assertAlmostEqual(ra[2], 360-30)
        self.assertAlmostEqual(dec[2], 0)


if __name__ == '__main__' :
    unittest.main()
