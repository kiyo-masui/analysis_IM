"""Unit tests for algebra.py."""

import os
import unittest

import scipy as sp
import numpy.lib.format as npfor

import algebra
import kiyopy.custom_exceptions as ce

class TestMatTypes(unittest.TestCase) :
    
    def test_from_memory(self) :
        # Works if constructed from array.
        data = sp.empty((5, 6, 6))
        data[:] = 4.0
        Mat = algebra.info_array(data, {'a' : 'b'})
        self.assertEqual(Mat.shape, (5, 6, 6))
        self.assertEqual(Mat.info['a'], 'b')
        self.assertTrue(sp.allclose(Mat, 4.0))
        # Check that the dictionary is there and empty if uninitialized.
        Mat4 = algebra.info_array(data)
        self.assertEqual(Mat4.info, {})
        # Works if constructed from a slice.
        Mat2 = Mat[1:2, :, :]
        self.assertEqual(Mat2.shape, (1, 6, 6))
        self.assertEqual(Mat2.info['a'], 'b')
        # Works if constructed from a view.
        Mat3 = data.view(algebra.info_array)
        self.assertEqual(Mat3.shape, (5, 6, 6))
        Mat3.info['a'] = 'b'

    def test_from_memmap(self) :
        # Works if constructed from array.
        data = npfor.open_memmap('temp.npy', mode='w+', shape=(4,3,3))
        data[:] = 5.0
        Mat = algebra.info_memmap(data, {'a': 'b'})
        Mat.flush()
        self.assertEqual(Mat.shape, (4, 3, 3))
        self.assertEqual(Mat.info['a'], 'b')
        self.assertTrue(sp.allclose(Mat, 5.0))
        self.assertTrue(isinstance(Mat,  sp.memmap))
        del Mat
        os.remove('temp.npy')

    def test_fails_if_not_memmap(self) :
        data = sp.empty((5, 7, 7))
        self.assertRaises(TypeError, algebra.info_memmap, data)

    def test_assert_info(self) :
        """Test the assert_info function."""
        # info_memaps should pass.
        data = npfor.open_memmap('temp.npy', mode='w+', shape=(4,3,3))
        data[:] = 5.0
        Mat = algebra.info_memmap(data)
        algebra.assert_info(Mat)
        del Mat
        os.remove('temp.npy')
        # info_arrays should pass.
        data = sp.empty((5, 6, 6))
        data[:] = 4.0
        Mat = algebra.info_array(data)
        algebra.assert_info(Mat)
        # arrays should fail.
        self.assertRaises(TypeError, algebra.assert_info, data)


class TestLoadSave(unittest.TestCase) :
    
    def setUp(self) :
        info = { 'a' : 42, 'b' : 'a string', 'c' : 3.14159 }
        data = sp.arange(80).reshape((2,8,5))
        self.Mat = algebra.info_array(data, info)

    def test_runs_filename(self) :
        algebra.save('temp.npy', self.Mat)
        self.assertTrue('temp.npy' in os.listdir('./'))
        self.assertTrue('temp.npy.meta' in os.listdir('./'))

    def test_runs_fid(self) :
        fid = open('temp.npy', 'w')
        algebra.save(fid, self.Mat)
        self.assertTrue('temp.npy' in os.listdir('./'))
        self.assertTrue('temp.npy.meta' in os.listdir('./'))

    def test_checks_dict_executable(self) :
        self.Mat.info['un_readable_object'] = sp.arange(100000)
        self.assertRaises(ce.DataError, algebra.save, 'temp.npy', self.Mat)

    def test_loads(self) :
        algebra.save('temp.npy', self.Mat)
        Loaded = algebra.load('temp.npy')
        self.assertTrue(isinstance(Loaded, sp.ndarray))
        self.assertTrue(isinstance(Loaded, algebra.info_array))
        self.assertTrue(sp.allclose(Loaded, self.Mat))
        self.assertEqual(Loaded.info['a'], self.Mat.info['a'])
        self.assertEqual(Loaded.info['b'], self.Mat.info['b'])
        self.assertAlmostEqual(Loaded.info['c'], self.Mat.info['c'])

    def test_memmap_read(self) :
        algebra.save('temp.npy', self.Mat)
        marray = algebra.open_memmap('temp.npy', mode="r")
        self.assertTrue(isinstance(marray, sp.ndarray))
        self.assertTrue(isinstance(marray, algebra.info_memmap))
        self.assertTrue(sp.allclose(marray, self.Mat))
        self.assertEqual(marray.info['a'], self.Mat.info['a'])

    def test_memmap_write(self) :
        marray = algebra.open_memmap('temp.npy', mode="w+",
                                     shape=self.Mat.shape,
                                     dtype=self.Mat.dtype,
                                     )
        marray[:] = self.Mat[:]
        marray.info = self.Mat.info
        del marray
        marray = None
        Loaded = algebra.load('temp.npy')
        self.assertTrue(isinstance(Loaded, sp.ndarray))
        self.assertTrue(isinstance(Loaded, algebra.info_array))
        self.assertTrue(sp.allclose(Loaded, self.Mat))
        self.assertEqual(Loaded.info['a'], self.Mat.info['a'])


    def tearDown(self) :
        del self.Mat
        if 'temp.npy' in os.listdir('./') :
            os.remove('temp.npy')
        if 'temp.npy.meta' in os.listdir('./') :
            os.remove('temp.npy.meta')


class TestMakeAlgebra(unittest.TestCase) :
    
    def setUp(self) :
        data = sp.arange(60, dtype=float)
        data.shape = (2, 5, 6)
        self.array = algebra.info_array(data)

    def test_make_vector(self) :
        algebra.make_vect(self.array)
        self.assertEqual(self.array.info['type'], 'vect')
        self.assertEqual(self.array.info['axes'], (None, None, None))

    def test_make_vector_named(self) :
        algebra.make_vect(self.array, ('freq', None, 'ra'))
        self.assertEqual(self.array.info['axes'], ('freq', None, 'ra'))

    def test_make_vector_raises(self) :
        self.assertRaises(TypeError, algebra.make_vect, self.array, (1,2,3))
        self.assertRaises(ValueError, algebra.make_vect, self.array, ('a','b'))

    def test_make_mat(self) :
        algebra.make_mat(self.array, row_axes=(0,1), col_axes=(0,2))
        self.assertEqual(self.array.info['type'], 'mat')
        self.assertEqual(self.array.info['axes'], (None, None, None))
        self.assertEqual(self.array.info['rows'], (0,1))
        self.assertEqual(self.array.info['cols'], (0,2))

    def test_make_mat_raises(self) :
        self.assertRaises(ValueError, algebra.make_mat, self.array, (1,), (2,))
        self.assertRaises(ValueError, algebra.make_mat, self.array, (0,1,'a'), 
                          (2,))
        self.assertRaises(ValueError, algebra.make_mat, self.array, (0,1), 
                          (3,))


class TestMatUtils(unittest.TestCase) :
    
    def test_mat_shape(self) :
        mat = sp.empty((3, 8, 6))
        self.assertEqual(algebra.mat_shape(mat), (24, 18))
        memmat = sp.memmap('temp.npy', shape=(3, 8, 6), dtype=sp.float32,
                           mode='w+')
        self.assertEqual(algebra.mat_shape(memmat), (24, 18))
        del memmat
        os.remove('temp.npy')
        mat2 = list(mat)
        self.assertRaises(TypeError, algebra.mat_shape, mat2)
        mat3 = sp.empty((7,8))
        self.assertRaises(TypeError, algebra.mat_shape, mat3)
       

if __name__ == '__main__' :
    unittest.main()
