"""Unit tests for psrfits_to_sdfits.py."""

import unittest
import sys

import scipy as sp
import numpy.random as rand

import psrfits_to_sdfits as p2s
import kiyopy.custom_exceptions as ce


class TestFormatData(unittest.TestCase) :

    def setUp(self) :
        self.ntime = 5
        self.npol = 4
        self.nfreq = 10
        self.good_data = sp.empty((self.ntime, self.npol, self.nfreq),
                                  dtype=int)
        self.good_data[:,:,:] = sp.reshape(sp.arange(self.ntime*self.nfreq),
                                           (self.ntime, 1, self.nfreq))
        self.good_data[:,0,:] +=  100
        self.good_data[:,1:,:] -= self.ntime*self.nfreq//2
        self.raw_data = sp.empty((self.ntime, self.npol, self.nfreq),
                                  dtype=sp.uint8)
        self.raw_data[:,0,:] = self.good_data[:,0,:]
        self.raw_data.dtype = sp.int8
        self.raw_data[:,1:,:] = self.good_data[:,1:,:]
        self.raw_data.dtype = sp.uint8
        self.raw_data = self.raw_data.flatten()
    
    def test_runs(self) :
        p2s.format_data(self.raw_data, self.ntime, self.npol, self.nfreq)

    def test_requires_uint(self) :
        self.assertRaises(TypeError, p2s.format_data, self.good_data,
                          self.ntime, self.npol, self.nfreq)

    def test_right_answer(self):
        reformated = p2s.format_data(self.raw_data, self.ntime, self.npol,
                                     self.nfreq)
        self.assertTrue(sp.allclose(reformated, self.good_data))

class TestFoldOnCal(unittest.TestCase) :

    def setUp(self) :
        self.ntime = 2048
        self.nfreq = 10
        self.data = sp.zeros((self.ntime, 4, self.nfreq))
        self.n_bins_cal = 64
        self.offset = 23 # Make it less than n_bins_cal/2
        self.level = 0.001*(5*self.nfreq + sp.arange(self.nfreq))
        self.data[:,:,:] += rand.normal(loc=0, scale=self.level/10,
                                        size=(self.ntime, 4, self.nfreq))
        for ii in range(self.ntime//self.n_bins_cal) :
            s = slice(self.offset + ii*self.n_bins_cal, self.offset +
                      (2*ii+1)*self.n_bins_cal//2)
            self.data[s, :, :] += self.level
        self.t_slice = slice(self.offset, sys.maxint, self.n_bins_cal//2)
        self.t_vals = 0.5*self.level + rand.normal(loc=0, scale=self.level/10,
                            size=(2*self.ntime//self.n_bins_cal, 4, self.nfreq))

    def test_runs(self) : 
        p2s.get_cal_mask(self.data, self.n_bins_cal)

    def test_right_answer_basic(self) :
        first_ind_on, n_blank = p2s.get_cal_mask(self.data, self.n_bins_cal)
        self.assertEqual(first_ind_on, self.offset+1)
        self.assertEqual(n_blank, 2)

    def test_right_answer_partial(self) :
        self.data[self.t_slice, :, :] = self.t_vals
        first_ind_on, n_blank = p2s.get_cal_mask(self.data, self.n_bins_cal)
        self.assertEqual(first_ind_on, self.offset+1)
        self.assertEqual(n_blank, 1)

    def test_checks_cal_per(self) :
        self.assertRaises(ValueError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal + 1)

    def test_fails_to_many_transitions(self) :
        self.data[self.t_slice, :, :] = self.t_vals

        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal*2)
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal//2)

    def test_fails_any_nan(self) :
        self.data[self.t_slice,0,:] = float('nan')
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fails_offs_in_ons(self) :
        self.data[self.t_slice, :, :] = self.t_vals
        
        s = slice(self.offset+7, sys.maxint, self.n_bins_cal)
        self.data[s, :, :] *= 0
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fails_late_on(self) :
        self.data[self.t_slice, :, :] = 0
        s = slice(self.offset+1, sys.maxint, self.n_bins_cal)
        self.data[s, :, :] = 0
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fails_to_many_semi_bins(self) :
        self.data[self.t_slice, :, :] = self.t_vals

        s = slice(self.offset+7, sys.maxint, self.n_bins_cal)
        self.data[s, :, :] *= 0.7
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

class TestSeparateCal(unittest.TestCase) :
    """Unlike the tests for get_cal_mask, these tests are tihgtly controled
    with no noise so we can detect deviations from expected."""
    
    def setUp(self) :
        self.ntime = 2048
        self.nfreq = 10
        self.data = sp.zeros((self.ntime, 4, self.nfreq))
        self.n_bins_cal = 64
        self.offset = 10

    def post_setup(self) :
        if self.offset > self.n_bins_cal//2 :
            last_on_start = (self.offset + self.n_bins_cal//2)% self.n_bins_cal
            self.data[:last_on_start, :, :] = 1
        for ii in range(self.ntime//self.n_bins_cal) :
            s = slice(self.offset + ii*self.n_bins_cal, self.offset +
                      (2*ii+1)*self.n_bins_cal//2)
            self.data[s, :, :] = 1
        self.t_slice_on = slice(self.offset, sys.maxint, self.n_bins_cal)
        self.t_slice_off = slice((self.offset +
                                  self.n_bins_cal//2)%self.n_bins_cal,
                                 sys.maxint, self.n_bins_cal)

    def check_answer(self) :
        data = p2s.separate_cal(self.data, self.n_bins_cal)
        self.assertTrue(sp.allclose(data[:,:,0,:], 1.0))
        self.assertTrue(sp.allclose(data[:,:,1,:], 0.0))

    def test_works_no_transition(self) :
        self.post_setup()
        self.check_answer()

    def test_works_transition(self) :
        self.post_setup()
        self.data[self.t_slice_off, :, :] = 0.3
        self.data[self.t_slice_on, :, :] = 0.5
        self.check_answer()
    
    # Move the offset to the the second half and make sure it works.

    def test_works_no_transition_late(self) :
        self.offset = 57
        self.post_setup()
        self.check_answer()

    def test_works_transition_late(self) :
        self.offset = 57
        self.post_setup()
        self.data[self.t_slice_off, :, :] = 0.3
        self.data[self.t_slice_on, :, :] = 0.5
        self.check_answer()

    # Test offset = 63
    def test_works_no_transition__1(self) :
        self.offset = 63
        self.post_setup()
        self.check_answer()

    def test_works_transition__1(self) :
        self.offset = 63
        self.post_setup()
        self.data[self.t_slice_off, :, :] = 0.3
        self.data[self.t_slice_on, :, :] = 0.5
        self.check_answer()
    
    # Test offset = 32
    def test_works_no_transition_32(self) :
        self.offset = 32
        self.post_setup()
        self.check_answer()

    def test_works_transition_32(self) :
        self.offset = 32
        self.post_setup()
        self.data[self.t_slice_off, :, :] = 0.3
        self.data[self.t_slice_on, :, :] = 0.5
        self.check_answer()

    # Test offset = 0
    def test_works_no_transition_0(self) :
        self.offset = 0
        self.post_setup()
        self.check_answer()

    def test_works_transition_0(self) :
        self.offset = 0
        self.post_setup()
        self.data[self.t_slice_off, :, :] = 0.3
        self.data[self.t_slice_on, :, :] = 0.5
        self.check_answer()




    # test offset = 32 and 63!!!



if __name__ == '__main__' :
        unittest.main()
