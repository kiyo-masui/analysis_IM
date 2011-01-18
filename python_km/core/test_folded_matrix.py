"""Unit tests for folded_matrix.py."""

import os
import unittest

import scipy as sp

import folded_matrix as fm
import kiyopy.custom_exceptions as ce


class TestInit(unittest.TestCase) :
    
    def test_init_from_scratch(self) :
        Mat = fm.FoldedMatrix((20, 40), (16, 5, 5))
        self.assertEqual(Mat.shape, (20,40))
        self.assertEqual(Mat.data.shape, (16, 5, 5))

    def test_init_memap(self) :
        Mat = fm.FoldedMatrix((20, 40), (16, 5, 5), 'temp.npz')
        self.assertTrue(isinstance(Mat.data, sp.memmap))
        self.assertEqual(Mat.shape, (20,40))
        self.assertEqual(Mat.data.shape, (16, 5, 5))
        os.remove('temp.npz')

    def test_init_array(self) :
        a = sp.empty((16, 5, 5))
        Mat = fm.FoldedMatrix((20, 40), a, 'temp.npz')
        del a
        self.assertEqual(Mat.shape, (20,40))
        self.assertEqual(Mat.data.shape, (16, 5, 5))

    def test_demands_right_dimensions(self) :
        self.assertRaises(ValueError, fm.FoldedMatrix, (20, 40, 1), (16, 5, 5))
        self.assertRaises(ValueError, fm.FoldedMatrix, (20, 40), (16, 55))
        self.assertRaises(ValueError, fm.FoldedMatrix, (20, 40),
                          sp.ones((16,55)))
        

#class Test_
        

if __name__ == '__main__' :
    unittest.main()
