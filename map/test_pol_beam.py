import unittest
from math import sqrt, exp, log

import numpy as np

import pol_beam

class TestSimple(unittest.TestCase):

    def setUp(self):
        """Set up constructs a beam of known parameters."""
    
        n_f = 6
        self.n_f = n_f
        self.freq = 10. + np.arange(n_f)
        Beam = pol_beam.SimpleBeam(self.freq)
        self.width = 0.9 + (np.arange(n_f) / 5. / n_f)
        self.sigma = self.width / (2 * np.sqrt(2 * np.log(2)))
        Beam.set_width(self.width)
        self.Beam = Beam
            
    def test_0_mode(self):
        coef = np.ones((self.n_f, 1, 1, 4))
        self.Beam.set_coefficients(coef)
        # Test that centre beam is unity.
        self.assertTrue(np.allclose(self.Beam.get_skewer(0, 0), 1))
        # Check that we get 1/2 at the FWHM.
        s = self.Beam.get_skewer(self.width/2, 0)
        self.assertTrue(np.allclose(s, 0.5))
        s = self.Beam.get_skewer(-self.width/2/sqrt(2), self.width/2/sqrt(2))
        self.assertTrue(np.allclose(s, 0.5))
        # Check the normalization.
        n_grid = 30
        beam = self.Beam.get_full(n_grid, 5)
        norm = np.sum(np.sum(beam**2, 2), 2) * (5. / n_grid)**2
        self.assertTrue(np.allclose(norm, (self.sigma**2 * np.pi)[:,None]))

    def test_ortho(self):
        coef = np.zeros((self.n_f, 5, 5, 4))
        # Set each polarization to a different polynomial.
        coef[:,0,0,0] = 1.
        coef[:,0,1,1] = 1.
        coef[:,1,1,2] = 1.
        coef[:,2,0,3] = 1.
        self.Beam.set_coefficients(coef)
        # Get grid representation of beam.
        n_grid = 20
        beam = self.Beam.get_full(n_grid, 4)
        # Take all 16 * n_f inner products.
        prods = np.sum(np.sum(beam[:,:,None,:,:] * beam[:,None,:,:,:], 3), 3)
        prods *= (4. / n_grid)**2
        # Take out the normalization.
        prods /= (self.sigma**2 * np.pi)[:,None,None]
        self.assertTrue(np.allclose(prods, np.eye(4)[None,:,:]))

class TestHermiteBasis(unittest.TestCase):

    def setUp(self):
        n_f = 6
        self.n_f = n_f
        self.freq = 10. + np.arange(n_f)
        self.width = 0.9 + (np.arange(n_f) / 5. / n_f)
        self.sigma = self.width / (2 * np.sqrt(2 * np.log(2)))
        self.center = np.arange(n_f)[:,None] * 0.2 / n_f * np.array([1, 2])
        self.Basis = pol_beam.HermiteBasis(self.freq, self.center, self.width)

    def test_center(self):
        # Test that the center of the basis is unity for the zero mode.
        center_skewer = self.Basis.eval_basis((0, 0), self.center[:,0,None], 
                                              self.center[:,1,None])
        self.assertTrue(np.allclose(center_skewer, 1.))
        # Test that it is zero for all odd modes.
        skewer = self.Basis.eval_basis((1, 0), self.center[:,0,None], 0)
        self.assertTrue(np.allclose(skewer, 0.))
        skewer = self.Basis.eval_basis((4, 7), 1, self.center[:,1,None])
        self.assertTrue(np.allclose(skewer, 0.))

    def test_orthonormal(self):
        # Modes we will test for orthonormality.
        mode_list = [(0, 0), (1, 0), (0, 2), (1, 2), (2, 5), (6, 8)]
        n_modes = len(mode_list)
        n_side = 35
        size = 6.
        mode_maps = np.empty((self.n_f, n_modes, n_side, n_side))
        for ii in range(n_modes):
            mode_maps[:,ii,:,:] = self.Basis.get_basis_grid(mode_list[ii],
                                                            n_side, size)
        # Take inner produce of all mode pairs.
        prods = np.sum(np.sum(mode_maps[:,:,None,:,:]
                              * mode_maps[:,None,:,:,:], -1), -1)
        # Normalize the integral.
        prods *= (size / n_side)**2
        # Gaussian normalizations
        prods /= (self.sigma**2 * np.pi)[:,None,None]
        # Check orthonormality.
        self.assertTrue(np.allclose(prods, np.eye(n_modes)))



if __name__ == '__main__' :
    unittest.main()

