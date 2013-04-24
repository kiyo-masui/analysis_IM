"""Unit tests for hdf5 IO for algebra.py."""

import unittest
import os

import numpy as np
import h5py

import algebra


class TestH5(unittest.TestCase):
    
    def setUp(self):
        info = { 'a' : 42, 'b' : 'a string', 'c' : 3.14159 }
        data = np.arange(80).reshape((2,8,5))
        self.Mat = algebra.info_array(data, info)

    def test_save(self):
        F = h5py.File('temp_out.hdf5', 'w')
        algebra.save_h5(F, 'Mat', self.Mat)
        # Check that the object is there.
        self.assertTrue('Mat' in F.keys())
        # Make sure it has the same data.
        self.assertTrue(np.all(F['Mat'][...] == self.Mat))
        # Check that the meta data is all there, but in string form.
        self.assertEqual(F['Mat'].attrs['a'], '42')
        self.assertEqual(F['Mat'].attrs['b'], "'a string'")

    def test_save_deep(self):
        F = h5py.File('temp_out.hdf5', 'w')
        algebra.save_h5(F, 'here/there/Mat', self.Mat)
        self.assertTrue('here' in F.keys())
        self.assertTrue(np.all(F['here/there/Mat'][...] == self.Mat))

    def test_load(self):
        F = h5py.File('temp_out.hdf5', 'w')
        algebra.save_h5(F, 'Mat', self.Mat)
        del F
        F = h5py.File('temp_out.hdf5', 'r')
        new_Mat = algebra.load_h5(F, 'Mat')
        self.assertTrue(np.all(new_Mat == self.Mat))
        self.assertTrue(new_Mat.info == self.Mat.info)

    def tearDown(self):
        if 'temp_out.hdf5' in os.listdir('./'):
            os.remove('temp_out.hdf5')


if __name__ == '__main__' :
    unittest.main()



