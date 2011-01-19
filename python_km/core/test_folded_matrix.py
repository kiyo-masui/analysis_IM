"""Unit tests for folded_matrix.py."""

import os
import unittest

import scipy as sp

import folded_matrix as fm
import kiyopy.custom_exceptions as ce


class TestInit(unittest.TestCase) :
    
    def test_init_from_scratch(self) :
        Mat = fm.FoldedMatrix((5, 4, 8))
        self.assertEqual(Mat.shape, (20, 40))
        self.assertEqual(Mat.data.shape, (5, 4, 8))

    def test_init_memap(self) :
        Mat = fm.FoldedMatrix((10, 2, 4), 'temp.npz')
        self.assertTrue(isinstance(Mat.data, sp.memmap))
        self.assertEqual(Mat.shape, (20, 40))
        self.assertEqual(Mat.data.shape, (10, 2, 4))
        os.remove('temp.npz')

    def test_init_array(self) :
        a = sp.empty((15, 2, 3))
        Mat = fm.FoldedMatrix(a)
        del a
        self.assertEqual(Mat.shape, (30, 45))
        self.assertEqual(Mat.data.shape, (15, 2, 3))

    def test_init_2D(self) :
        # Initilize empty.
        Mat = fm.FoldedMatrix((4, 8))
        self.assertEqual(Mat.shape, (4, 8))
        self.assertEqual(Mat.data.shape, (1, 4, 8))
        # And from an array.
        a = sp.empty((15, 7))
        Mat = fm.FoldedMatrix(a)
        del a
        self.assertEqual(Mat.shape, (15, 7))
        self.assertEqual(Mat.data.shape, (1, 15, 7))

if __name__ == '__main__' :
    unittest.main()
