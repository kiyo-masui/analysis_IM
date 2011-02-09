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
        self.array = algebra.vect_array(self.array)
        self.assertEqual(self.array.info['type'], 'vect')
        self.assertEqual(self.array.info['axes'], (None, None, None))
        self.assertTrue(isinstance(self.array, algebra.vect))
        self.assertTrue(isinstance(self.array, algebra.info_array))

    def test_make_vector_named(self) :
        self.array = algebra.vect_array(self.array, ('freq', None, 'ra'))
        self.assertEqual(self.array.info['axes'], ('freq', None, 'ra'))
        self.assertEqual(self.array.axes, ('freq', None, 'ra'))

    def test_make_vector_raises(self) :
        self.assertRaises(TypeError, algebra.vect_array, self.array, (1,2,3))
        self.assertRaises(ValueError, algebra.vect_array, self.array, 
                          ('a','b'))

    def test_make_mat(self) :
        self.array = algebra.mat_array(self.array, row_axes=(0,1), 
                                       col_axes=(0,2))
        self.assertEqual(self.array.info['type'], 'mat')
        self.assertEqual(self.array.info['axes'], (None, None, None))
        self.assertEqual(self.array.info['rows'], (0,1))
        self.assertEqual(self.array.info['cols'], (0,2))
        self.assertEqual(self.array.cols, (0,2))

    def test_make_mat_raises(self) :
        self.assertRaises(ValueError, algebra.mat_array, self.array, (1,), (2,))
        self.assertRaises(ValueError, algebra.mat_array, self.array, (0,1,'a'), 
                          (2,))
        self.assertRaises(ValueError, algebra.mat_array, self.array, (0,1), 
                          (3,))


class TestMatUtils(unittest.TestCase) :
    
    def setUp(self) :
        data = sp.arange(30)
        data.shape = (5, 2, 3)
        self.vect = algebra.info_array(data)
        self.vect = algebra.vect_array(self.vect, 
                                       axis_names=('freq', 'a', 'b'))
        data = sp.arange(120)
        self.mat = algebra.info_array(data)
        self.mat.shape = (5, 4, 6)
        self.mat = algebra.mat_array(self.mat, row_axes=(0,1), col_axes=(0,2), 
                                     axis_names=('freq', 'mode', 'mode'))

    def test_mat_shape(self) :
        # Check basic functionality.
        self.assertEqual(self.vect.mat_shape, (30,))
        self.assertEqual(self.mat.mat_shape, (20, 30))
        # Check that it fails if I mess with stuff.
        self.vect.info['axes'] = (None, 5, 'a')
        self.assertRaises(TypeError, self.vect.__getattr__, 'mat_shape')
        self.vect.info['axes'] = (None, 'a')
        self.assertRaises(ValueError, self.vect.__getattr__, 'mat_shape')
        self.mat.info['rows'] = (0,)
        self.assertRaises(ValueError, self.mat.__getattr__, 'mat_shape')
        self.mat.info['rows'] = (0,1,4)
        self.assertRaises(ValueError, self.mat.__getattr__, 'mat_shape')

    def test_assert_axes_ordered(self) :
        algebra.assert_axes_ordered(self.mat)
        self.mat.info['cols'] = (2,)
        algebra.assert_axes_ordered(self.mat)
        self.mat.info['cols'] = (0,2)
        self.mat.info['rows'] = (1,0)
        self.assertRaises(NotImplementedError, algebra.assert_axes_ordered, 
                          self.mat)
        self.mat.info['rows'] = (0,2)
        self.mat.info['cols'] = (0,1)
        self.assertRaises(NotImplementedError, algebra.assert_axes_ordered, 
                          self.mat)

    def test_expand_mat(self) :
        expanded = algebra.expand_mat(self.mat)
        self.assertEqual(expanded.shape, algebra.get_shape(self.mat))
        self.assertTrue(sp.allclose(self.mat[0,:,:], expanded[0:4, 0:6]))

    def test_dot_mat_vect(self) :
        self.mat.shape = (5, 4, 2, 3)
        self.mat.info['cols'] = (0, 2, 3)
        self.mat.info['axes'] = ('freq', 'mode', 'a', 'b')
        prod = algebra.dot(self.mat, self.vect)
        self.assertEqual(algebra.get_shape(prod), (20,))
        self.assertEqual(prod.shape, (5, 4))
        self.assertTrue(sp.allclose(prod.flatten(),
                        sp.dot(algebra.expand_mat(self.mat),
                        self.vect.flatten())))
        #self.assertEqual(prod.info['axes'], ('freq', 'mode'))

if __name__ == '__main__' :
    unittest.main()
