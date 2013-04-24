"""Unit tests for utils module."""

import unittest

import scipy as sp
import numpy as np
from numpy import random
from scipy import linalg, special

import misc as utils

import matplotlib.pyplot as plt

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
        c_az, c_el = utils.radec2azelGBT(ra, dec, UT)
        self.assertAlmostEqual(c_el, el, 1)
        
        el = 38.43312
        az = 0.
        ra, dec = utils.elaz2radecGBT(el, az, UT)
        self.assertAlmostEqual(dec, 90, 1)
        c_az, c_el = utils.radec2azelGBT(ra, dec, UT)
        self.assertAlmostEqual((c_az + 180) % 360, (az + 180) % 360, 1)
        self.assertAlmostEqual(c_el, el, 1)
        
        el = 90 - 38.43312
        az = 180.
        ra, dec = utils.elaz2radecGBT(el, az, UT)
        self.assertAlmostEqual(dec, 0, 1)
        c_az, c_el = utils.radec2azelGBT(ra, dec, UT)
        self.assertAlmostEqual(c_az, az, 1)
        self.assertAlmostEqual(c_el, el, 1)

class TestAzEl2P_GBT(unittest.TestCase):

    def nieve_conversion(self, az, el):
        """Basic conversion. Ignores precesions and only works for 
        -90 < p < 90."""
        
        az_r = az * np.pi / 180
        el_r = el * np.pi / 180
        lat_r = (38. + 25./60 + 59.23/360) * np.pi / 180

        sin_dec = np.sin(lat_r) * np.sin(el_r)
        sin_dec += np.cos(lat_r) * np.cos(el_r) * np.cos(az_r)
        dec_r = np.arcsin(sin_dec)

        sin_p = -np.sin(az_r) * np.cos(lat_r) / np.cos(dec_r)
        return np.arcsin(sin_p) * 180. / np.pi

    def test_known_cases(self):
        UT = '2000-01-01T12:00:00.00'
        self.assertAlmostEqual(utils.azel2pGBT(0, 60, UT) % 360, 180, 1)
        self.assertAlmostEqual(utils.azel2pGBT(90, 89.99, UT) % 360, 270, 1)
        self.assertAlmostEqual(utils.azel2pGBT(270, 89.99, UT) % 360, 90, 1)

    def test_simple_conversion(self):
        UT = '2000-01-01T12:00:00.00'
        self.assertAlmostEqual(utils.azel2pGBT(120, 88, UT) % 360, 
                               self.nieve_conversion(120, 88) % 360, 0)
        self.assertAlmostEqual(utils.azel2pGBT(15, 5, UT) % 360, 
                               self.nieve_conversion(15, 5) % 360, 0)
        self.assertAlmostEqual(utils.azel2pGBT(210, 30, UT) % 360, 
                               self.nieve_conversion(210, 30) % 360, 0)


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

class TestAmpFit(unittest.TestCase):

    def test_uncorrelated_noscatter(self):
        data = sp.arange(10, dtype=float)
        theory = data/2.0
        C = sp.identity(10)
        out = utils.ampfit(data, C, theory)
        a, s = out['amp'], out['error']
        self.assertAlmostEqual(a, 2)

    def test_uncorrelated_noscatter_error(self):
        data = sp.arange(10, dtype=float) + 1.0
        amp = 5.0
        theory = data/amp
        C = sp.diag((sp.arange(10, dtype=float) + 1.0)**2)
        out = utils.ampfit(data, C, theory)
        a, s = out['amp'], out['error']
        self.assertAlmostEqual(a, amp)
        self.assertAlmostEqual(s/a, 1.0/sp.sqrt(10))

    def test_uncorrelated_scatter_error(self):
        n = 1000
        data = sp.arange(n, dtype=float) + 1.0
        amp = 5.0
        theory = data/amp
        C = sp.diag((sp.arange(n, dtype=float) + 1.0)**2)
        data += random.normal(size=n)*data
        out = utils.ampfit(data, C, theory)
        a, s = out['amp'], out['error']
        self.assertTrue(sp.allclose(a, amp, rtol=5.0/sp.sqrt(n), atol=0))
        # Expect the next line to fail 1/100 trials.
        self.assertFalse(sp.allclose(a, amp, rtol=0.01/sp.sqrt(n), atol=0))
        self.assertAlmostEqual(s, amp/sp.sqrt(1000))

    def test_correlated_scatter(self) :
        n = 50
        r = (sp.arange(n, dtype=float) + 10.0*n)/10.0*n
        data = sp.sin(sp.arange(n)) * r 
        amp = 25.0
        theory = data/amp
        # Generate correlated matrix.
        C = random.rand(n, n) # [0, 1) 
        # Raise to high power to make values near 1 rare.
        C = (C**10) * 0.2
        C = (C + C.T)/2.0
        C += sp.identity(n)
        C *= r[:, None]/2.0
        C *= r[None, :]/2.0
        # Generate random numbers in diagonal frame.
        h, R = linalg.eigh(C)
        self.assertTrue(sp.alltrue(h>0))
        rand_vals = random.normal(size=n)*sp.sqrt(h)
        # Rotate back.
        data += sp.dot(R.T, rand_vals)
        out = utils.ampfit(data, C, theory)
        a, s = out['amp'], out['error']
        self.assertTrue(sp.allclose(a, amp, atol=5.0*s, rtol=0))
        # Expect the next line to fail 1/100 trials.
        self.assertFalse(sp.allclose(a, amp, atol=0.01*s, rtol=0))

class TestFloatTime(unittest.TestCase):

    def test_near_epoch(self):
        UT = "2000-01-01T00:00:54.51"
        self.assertAlmostEqual(54.51, utils.time2float(UT))
        self.assertEqual(utils.float2time(54.51), UT)

    def test_full_circle(self):
        seconds = 451732642.56
        UT = utils.float2time(seconds)
        seconds2 = utils.time2float(UT)
        self.assertEqual(UT, utils.float2time(seconds2))
        self.assertAlmostEqual(seconds, seconds2)

    def test_full_circle_vectorized(self):
        seconds = sp.array([[34252672.12, 623421.65], [23464.1, 644656784.56]])
        UT = utils.float2time(seconds)
        seconds2 = utils.time2float(UT)
        self.assertTrue(sp.all(UT == utils.float2time(seconds2)))
        self.assertTrue(sp.allclose(seconds, seconds2))

class Test_OrthoPoly(unittest.TestCase):

    def test_flat_even_spaced(self):
        # Use uniform weight, should just get the lagendre polynomials.
        m = 20
        n = 10
        x = sp.arange(m, dtype=float)
        window = 1.
        polys = utils.ortho_poly(x, n, window)
        # The first one should just be the mean.
        self.assertTrue(sp.allclose(polys[0,:], 1./sp.sqrt(m)))
        # The second one should be a slope.
        expected = x - sp.mean(x)
        expected = expected / sp.sqrt(sp.sum(expected**2))
        self.assertTrue(sp.allclose(polys[1,:], expected))
        # Check that they are all orthonormal.
        self.check_ortho_norm(polys, 1.)

    def test_uneven(self):
        # Use uniform weight, should just get the lagendre polynomials.
        # Get lots of polynomials to test the numerical stability.
        m = 300
        n = 10
        x = sp.log(sp.arange(m, dtype=float)/2 + 0.5)
        window = sp.sin(x)**2
        polys = utils.ortho_poly(x, n, window)
        self.check_ortho_norm(polys, window)
        #plt.plot(x, window, '.')
        #plt.plot(x, polys[0,:])
        #plt.plot(x, polys[1,:])
        #plt.plot(x, polys[2,:])
        #plt.plot(x, polys[3,:])
        #plt.plot(x, polys[6,:])
        #plt.show()

    def test_multiD(self):
        # Use uniform weight, should just get the lagendre polynomials.
        m = 40
        n = 10
        x = sp.log(sp.arange(m, dtype=float)/2 + 0.5)
        window = sp.empty((m, 2), dtype=float)
        window[:,0] = sp.sin(x)**2
        window[:,1] = sp.cos(x)**2
        x.shape = (m, 1)
        polys = utils.ortho_poly(x, n, window, axis=0)
        self.check_ortho_norm(polys, window, axis=0)

    def test_multiD2(self):
        # Use uniform weight, should just get the lagendre polynomials.
        m = 40
        n = 10
        x = sp.log(sp.arange(m, dtype=float)/2 + 0.5)
        window = sp.empty((2, m), dtype=float)
        window[0,:] = sp.sin(x)**2
        window[1,:] = sp.cos(x)**2
        x.shape = (1, m)
        polys = utils.ortho_poly(x, n, window, axis=1)
        self.check_ortho_norm(polys, window, axis=1)

    def check_ortho_norm(self, polys, window=1., axis=-1):
        # Always check that they are all orthonormal.
        n = polys.shape[0]
        for ii in range(n):
            for jj in range(n):
                prod = sp.sum(window * polys[ii,:] * polys[jj,:], axis)
                if ii == jj:
                    self.assertTrue(sp.alltrue(abs(prod - 1.) < 1e-8))
                else:
                    self.assertTrue(sp.alltrue(abs(prod) < 1e-8))

class TestOrthoPoly2D(unittest.TestCase):
    """Unit tests for the 2D orthonormal polynomial generator."""

    def test_raise_errors_shape(self):
        # Not 2D.
        self.assertRaises(ValueError, utils.ortho_poly_2D, np.zeros(5), 
                          np.zeros(5), 4, np.zeros(5))
        # Not broadcastable.
        self.assertRaises(ValueError, utils.ortho_poly_2D, np.zeros((5, 5)), 
                          np.zeros(3), 4, np.zeros(5))
        self.assertRaises(ValueError, utils.ortho_poly_2D, np.zeros((5, 5)), 
                          np.zeros((5, 1)), 4, np.zeros((5, 6)))
        self.assertRaises(ValueError, utils.ortho_poly_2D, np.zeros((5, 5)), 
                          np.zeros((1, 5)), 4, np.zeros((5, 5, 1)))
        # Make sure the output is the right shape.
        x = np.arange(5)
        x.shape = (5, 1)
        y = np.arange(7)
        y.shape = (1, 7)
        #out = utils.ortho_poly_2D(x, y, 3, 1.)
        #self.assertEquals(out.shape, (3, 3, 5, 7))

    def test_gets_legendre(self):
        # Define problem size.
        m = 100
        n = 8
        x = np.arange(m)
        x.shape = (m, 1)
        y = np.arange(m)
        y.shape = (1, m)
        # Normalized domain.
        z = (np.arange(m, dtype=float) + 0.5) * 2. / m - 1.
        # Generate exact Legendre.
        legendre = np.zeros((n, n, m, m), dtype=float)
        for ii in range(n):
            for jj in range(n - ii):
                ljj = special.legendre(jj)
                lii = special.legendre(ii)
                this_leg = lii(z)
                this_leg.shape = (m, 1)
                this_leg = this_leg * ljj(z)
                # Normalize.
                this_leg /= np.sqrt(np.sum(this_leg**2))
                legendre[ii,jj,...] = this_leg
        # Generate the same with Gram-Schmidt.
        #out = utils.ortho_poly_2D(x, y, n, 1.)
        # Check that they are the same to within the a few times the number of
        # points. Absolute toerance ~1/m**2 because 1/m is the magnitude of 
        # the elements and I want an precision of 1/m.
        #plt.imshow(out[1, 2])
        #plt.colorbar()
        #plt.figure()
        #plt.imshow(legendre[1, 2])
        #plt.colorbar()
        #plt.figure()
        #plt.imshow(out[1, 2] - legendre[1, 2])
        #plt.colorbar()
        #plt.show()
        #self.assertTrue(np.allclose(out, legendre, atol=5./m**2))



if __name__ == '__main__' :
    unittest.main()
