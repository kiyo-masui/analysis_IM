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

    def setUp(self):
        self.ntime = 2048
        self.nfreq = 10
        self.data = sp.zeros((self.ntime, 4, self.nfreq))
        self.n_bins_cal = 64
        # Set channel dependant gain.
        self.level = 0.1*(self.nfreq + sp.arange(self.nfreq))
        # Add noise.
        self.data[:,:,:] += (0.1 * self.level
                             * rand.randn(self.ntime, 4, self.nfreq))
        # Add DC level.
        self.dc = 10 * self.level
        self.data += self.dc
        # First can transition.
        self.first_trans = rand.randint(0, self.n_bins_cal // 2)
        # The following randomly assigns self.neg to -1 or 1.
        self.neg = 0
        while not self.neg: self.neg = rand.randint(-1, 2)
        # First upward edge:
        if self.neg == 1:
            self.offset = self.first_trans
        else:
            self.offset = self.first_trans + self.n_bins_cal // 2
            self.data[:,0,:] += self.level
        for ii in range(self.ntime//self.n_bins_cal) :
            s = slice(self.first_trans + ii*self.n_bins_cal, self.first_trans +
                      (2*ii+1)*self.n_bins_cal//2)
            self.data[s, 0, :] += self.neg * self.level
        # Transition values and locations.
        self.t_slice = slice(self.first_trans, sys.maxint, self.n_bins_cal//2)
        self.t_vals = 0.5 + 0.1 * rand.randn(2*self.ntime//self.n_bins_cal,
                                             self.nfreq)
        self.t_vals *= - self.level

    def test_runs(self) : 
        p2s.get_cal_mask(self.data, self.n_bins_cal)

    def test_right_answer_basic(self) :
        first_ind_on, n_blank = p2s.get_cal_mask(self.data, self.n_bins_cal)
        self.assertEqual(first_ind_on, (self.offset + 1) % self.n_bins_cal)
        self.assertEqual(n_blank, 2)

    def test_right_answer_partial(self) :
        self.data[self.t_slice, 0, :] += self.t_vals
        first_ind_on, n_blank = p2s.get_cal_mask(self.data, self.n_bins_cal)
        self.assertEqual(first_ind_on, (self.offset + 1) % self.n_bins_cal)
        self.assertEqual(n_blank, 1)

    def test_checks_cal_per(self) :
        self.assertRaises(ValueError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal + 1)

    def test_fails_to_many_transitions(self) :
        self.data[self.t_slice, 0, :] += self.t_vals

        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal*2)
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal//2)

    def test_fails_any_nan(self) :
        self.data[self.t_slice,0,:] = float('nan')
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fails_offs_in_ons(self) :
        self.data[self.t_slice, 0, :] += self.t_vals
        
        s = slice((self.offset + 7) % self.n_bins_cal, sys.maxint,
                  self.n_bins_cal)
        self.data[s, :, :] = self.dc
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fails_late_on(self) :
        self.data[self.t_slice, 0, :] = self.dc
        s = slice(self.offset+1, sys.maxint, self.n_bins_cal)
        self.data[s, :, :] = self.dc
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fails_to_many_semi_bins(self) :
        self.data[self.t_slice, 0, :] += self.t_vals

        s = slice((self.offset + 7) % self.n_bins_cal, sys.maxint,
                  self.n_bins_cal)
        self.data[s, :, :] = self.dc + self.level * 0.7
        self.assertRaises(ce.DataError, p2s.get_cal_mask, self.data,
                          self.n_bins_cal)

    def test_fast_flagger(self):
        for ii in range(self.ntime * self.nfreq * 4 // self.n_bins_cal // 10):
        #for ii in range(3):
            i_f = rand.randint(0, self.nfreq)
            i_t = rand.randint(0, self.ntime)
            i_p = rand.randint(0, 4)
            self.data[i_t,i_p,i_f] += self.level[i_f] * 5 
        data, weights = p2s.separate_cal(self.data, self.n_bins_cal, flag=10)
        right_answer = sp.empty((4, 2, self.nfreq))
        right_answer[...] = self.dc
        right_answer[0,0,:] += self.level
        self.assertTrue(sp.allclose(data, right_answer, atol=self.level / 10))
        self.assertTrue(sp.all(weights <= 1.))
        kept_fraction = 1. - 4./self.n_bins_cal - (4./self.n_bins_cal/10) 
        self.assertTrue(sp.allclose(sp.mean(weights), kept_fraction, rtol=1e-3))
            


class TestSeparateCal(unittest.TestCase) :
    """Unlike the tests for get_cal_mask, these tests are tightly controled
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
        data = self.data.copy()
        outdata, weights = p2s.separate_cal(data, self.n_bins_cal, flag=-1)
        self.assertTrue(sp.allclose(outdata[:,:,0,:], 1.0))
        self.assertTrue(sp.allclose(outdata[:,:,1,:], 0.0))
        data = self.data.copy()
        outdata, weights = p2s.separate_cal(data, self.n_bins_cal, flag=10)
        self.assertTrue(sp.allclose(outdata[:,:,0,:], 1.0))
        self.assertTrue(sp.allclose(outdata[:,:,1,:], 0.0))

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


if __name__ == '__main__' :
        unittest.main()
