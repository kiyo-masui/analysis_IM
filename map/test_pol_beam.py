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



if __name__ == '__main__' :
    unittest.main()
