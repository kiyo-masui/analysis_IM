"""Unit tests for dirty map maker."""

import os
import glob

import unittest
import scipy as sp
import numpy.ma as ma
import numpy.random as rand
from scipy import linalg

import clean_map
import dirty_map
import core.algebra as al
import _cholesky as _c

class TestSolver(unittest.TestCase):

    def setUp(self):
        # Make a positive definite noise matrix, clean map, and dirty_map.
        self.nra = 10
        self.ndec = 5
        self.nf = 20
        self.shape = (self.nf, self.nra, self.ndec)
        self.size = self.nra * self.ndec * self.nf
        # Clean map.
        clean_map = sp.empty(self.shape, dtype=float)
        clean_map = al.make_vect(clean_map, axis_names=('freq', 'ra', 'dec'))
        clean_map[...] = sp.sin(sp.arange(self.nf))[:,None,None]
        clean_map *= sp.cos(sp.arange(self.nra))[:,None]
        clean_map *= sp.cos(sp.arange(self.ndec))
        # Noise inverse matrix.
        noise_inv = sp.empty(self.shape * 2, dtype=float)
        noise_inv = al.make_mat(noise_inv, axis_names=('freq', 'ra', 'dec')*2,
                                row_axes=(0, 1, 2), col_axes=(3, 4, 5))
        rand_mat = rand.randn(*((self.size,) * 2))
        information_factor = 1.e6  # K**-2
        rand_mat = sp.dot(rand_mat, rand_mat.transpose()) * information_factor
        noise_inv.flat[...] = rand_mat.flat
        # Dirty map.
        dirty_map = al.partial_dot(noise_inv, clean_map)
        # Store in self.
        self.clean_map = clean_map
        self.noise_inv = noise_inv
        self.dirty_map = dirty_map

    def eig(self):
        noise_evalsinv, noise_evects = linalg.eigh(sp.reshape(self.noise_inv, 
                                            (self.size, self.size)))
        self.noise_evalsinv = al.make_mat(noise_evalsinv, axis_names=('mode',),
                                     row_axes=(0,), col_axes=(0,))
        self.noise_evects = al.make_mat(sp.reshape(noise_evects, 
                                           (self.shape + (self.size,))),
                                   axis_names=('freq', 'ra', 'dec', 'mode'),
                                   row_axes=(0, 1, 2), col_axes=(3,))

        
    def test_tri_copy(self):
        self.noise_inv.shape = (self.size, self.size)
        tri_noise_inv = _c.up_tri_copy(self.noise_inv)
        for ii in range(self.size):
            self.assertTrue(sp.allclose(tri_noise_inv[ii,ii:],
                                        self.noise_inv[ii,ii:]))
        
    def test_tri_copy_memmap(self):
        self.noise_inv.shape = (self.size, self.size)
        noise_mem = al.open_memmap("testout.npy", mode='w+',
                                   shape=self.noise_inv.shape)
        noise_mem[...] = self.noise_inv
        tri_noise_inv = _c.up_tri_copy(noise_mem)
        for ii in range(self.size):
            self.assertTrue(sp.allclose(tri_noise_inv[ii,ii:],
                                        self.noise_inv[ii,ii:]))

    def test_cholesky(self):
        self.noise_inv.shape = (self.size, self.size)
        # Do a normal cholesky and compare.
        lapack_chol = linalg.cholesky(self.noise_inv)
        _c.call_cholesky(self.noise_inv)
        for ii in range(self.size):
            self.assertTrue(sp.allclose(self.noise_inv[ii,ii:],
                                        lapack_chol[ii,ii:]))

    def test_solve(self):
        noise_cpy = self.noise_inv.copy()
        new_clean_map, noise_diag, chol = clean_map.solve(self.noise_inv,
                                                          self.dirty_map, True)
        self.assertTrue(sp.allclose(self.noise_inv, noise_cpy))
        self.assertTrue(sp.allclose(new_clean_map, self.clean_map))
        # Test the noise diagonal.
        self.noise_inv.shape = (self.size, self.size)
        noise = linalg.inv(self.noise_inv)
        new_noise_diag = noise.flat[::self.size + 1]
        new_noise_diag.shape = self.shape
        self.assertTrue(sp.allclose(noise_diag, new_noise_diag))

    def test_solve_eig(self):
        # Eigen decomposition of the noise inverse.
        self.eig()
        # now to actually test.
        noise_evects_copy = self.noise_evects.copy()
        new_clean_map, noise_diag = clean_map.solve_from_eig(
            self.noise_evalsinv, self.noise_evects, self.dirty_map, True,
            feedback=0)
        self.assertTrue(sp.allclose(new_clean_map, self.clean_map))
        self.assertTrue(sp.allclose(self.noise_evects, noise_evects_copy))
        expected_noise_diag = linalg.inv(sp.reshape(self.noise_inv, 
                                                    (self.size, self.size)))
        expected_noise_diag = expected_noise_diag.flat[::self.size + 1]
        self.assertTrue(sp.allclose(expected_noise_diag, noise_diag.flat))

    def test_solve_eig_bad_ind(self):
        # Set all the information in one pixel to nil.
        self.noise_inv[17,3,1,...] = 0
        self.noise_inv[...,17,3,1] = 0
        self.dirty_map = al.partial_dot(self.noise_inv, self.clean_map)
        self.eig()
        new_clean_map, noise_diag = clean_map.solve_from_eig(
            self.noise_evalsinv, self.noise_evects, self.dirty_map, True,
            feedback=0)
        self.clean_map[17,3,1] = 0
        self.assertTrue(sp.allclose(new_clean_map, self.clean_map))


    def tearDown(self):
        files = glob.glob("testout*")
        for file in files:
            os.remove(file)


if __name__ == '__main__' :
    unittest.main()
