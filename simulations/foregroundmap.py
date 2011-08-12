
import numpy as np
import scipy.linalg as la

import gaussianfield

from maps import *


def matrix_root_manynull(mat, threshold = 1e-16):
    """Square root a matrix.

    An inefficient alternative to the Cholesky decomposition for a
    matrix with a large dynamic range in eigenvalues. Numerical
    roundoff causes Cholesky to fail as if the matrix were not
    positive semi-definite. This does an explicit eigen-decomposition,
    setting small and negative eigenvalue to zero.

    Parameters
    ==========
    mat - ndarray
        An N x N matrix to decompose.
    threshold : scalar, optional
        Set any eigenvalues a factor `threshold` smaller than the
        largest eigenvalue to zero.

    Returns
    =======
    root : ndarray
        The decomposed matrix. This is truncated to the number of
        non-zero eigen values.
    num_pos : integer
        The number of positive eigenvalues.
    """
    evals, evecs = la.eigh(mat)

    evals[np.where(evals < evals.max() * threshold)] = 0.0
    num_pos = len(np.flatnonzero(evals))

    return (evecs[:,-num_pos:] * evals[np.newaxis,-num_pos:]**0.5), num_pos

    

class ForegroundMap(Map3d):
    r"""Simulate foregrounds with a seperable angular and frequency
    covariance.

    Used to simulate foregrounds that can be modeled by angular
    covariances of the form:
    .. math:: C_l(\nu,\nu') = A_l B(\nu, \nu')
    """
    
    _weight_gen = False

    def angular_powerspectrum(self, karray):
        pass


    def frequency_covariance(self, nu1, nu2):
        pass


    def generate_weight(self, regen = False):
        r"""Pregenerate the k weights array.

        Parameters
        ==========
        regen : boolean, optional
            If True, force regeneration of the weights, to be used if
            parameters have been changed,
        """
        
        if self._weight_gen and not regen:
            return
        
        f1, f2 = np.meshgrid(self.nu_pixels, self.nu_pixels)

        ch = self.frequency_covariance(f1, f2)

        self._freq_weight, self._num_corr_freq = matrix_root_manynull(ch)

        rf = gaussianfield.RandomFieldA2.like_map(self)
        #rf.widthx = self.widthx
        #rf.widthy = self.widthy
        #rf.numx = self.numx
        #rf.numy = self.numy
        
        rf.powerspectrum = self.angular_powerspectrum

        self._ang_field = rf

        self._weight_gen = True

    def getfield(self):

        self.generate_weight()

        aff = np.fft.rfftn(self._ang_field.getfield())

        s2 = (self._num_corr_freq,) + aff.shape

        norm = np.tensordot(self._freq_weight, np.random.standard_normal(s2), axes = (1,0))

        return np.fft.irfft(np.fft.ifft(norm * aff[np.newaxis,:,:], axis = 1), axis=2)

        
        
        
        
