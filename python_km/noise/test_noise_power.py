"""Unit tests for noise_power module."""

import unittest

import scipy as sp
import scipy.fftpack as fft
from numpy import random
import scipy.linalg as linalg

import noise_power as np
from core import fitsGBT
import kiyopy.custom_exceptions as ce

testfile = "testfile_GBTfits.fits"

class TestWindowedPower(unittest.TestCase) :

    def setUp(self) :
        self.n = 10000
        self.mode1 = 677 # Choose the 677th mode.
        self.k1 = (2.0*sp.pi*self.mode1)/self.n
        self.amp1 = 3.45
        self.wave1 = self.amp1*sp.sin(self.k1*(sp.arange(self.n) + 6.45))

    def test_no_window(self) :
        window = sp.ones_like(self.wave1)
        power = np.windowed_power(self.wave1, window)
        self.assertAlmostEqual(power[self.mode1]/self.amp1**2/self.n*4, 1)
        self.assertTrue(sp.allclose(power[:self.mode1], 0, 
                                    atol=self.amp1**2*self.n/1e15))
        self.assertTrue(sp.allclose(power[self.mode1+1:], 0, 
                                    atol=self.amp1**2*self.n/1e15))
        # With no window, we have a quick way to the answer
        quick_power = (abs(fft.fft(self.wave1))**2)[:self.n//2]/self.n
        self.assertTrue(sp.allclose(power, quick_power))

    def test_non_zero_window(self) :
        window = sp.ones_like(self.wave1)
        window += 0.5*sp.sin(sp.arange(self.n)/15.0)
        self.wave1 *= window
        power = np.windowed_power(self.wave1, window)
        self.assertAlmostEqual(power[self.mode1]/self.amp1**2/self.n*4, 1.0, 3)
        self.assertTrue(sp.allclose(power[:self.mode1], 0, 
                                    atol=self.amp1**2*self.n/1e3))
        self.assertTrue(sp.allclose(power[self.mode1+1:], 0, 
                                    atol=self.amp1**2*self.n/1e3))

    def test_window_with_zeros(self) :
        window = sp.ones_like(self.wave1)
        window += sp.sin(sp.arange(self.n)/50.0)
        self.wave1 *= window
        power = np.windowed_power(self.wave1, window)
        self.assertAlmostEqual(power[self.mode1]/self.amp1**2/self.n*4, 1.0, 3)
        self.assertTrue(sp.allclose(power[:self.mode1], 0, 
                                    atol=self.amp1**2*self.n/1e3))
        self.assertTrue(sp.allclose(power[self.mode1+1:], 0, 
                                    atol=self.amp1**2*self.n/1e3))

    def test_different_windows(self) :
        window1 = sp.ones_like(self.wave1)
        window1 += sp.sin(sp.arange(self.n)/40.0)
        window2 = 2*sp.ones_like(self.wave1)
        window2 += 2*sp.sin(sp.arange(self.n)/62.0)
        window2[window2<0.5] = 0
        wave1 = self.wave1*window1
        wave2 = self.wave1*window2
        power = np.windowed_power(wave1, window1, wave2, window2)
        self.assertAlmostEqual(power[self.mode1]/self.amp1**2/self.n*4, 1.0, 3)
        self.assertTrue(sp.allclose(power[:self.mode1], 0, 
                                    atol=self.amp1**2*self.n/1e3))
        self.assertTrue(sp.allclose(power[self.mode1+1:], 0, 
                                    atol=self.amp1**2*self.n/1e3))

    def test_statistical_no_window(self) :
        n_trials = 1000
        n_points = 200
        window = sp.ones(n_points, dtype=float)
        power = sp.zeros(n_points//2)
        for ii in range(n_trials) :
            wave = self.amp1*random.randn(n_points)
            power += np.windowed_power(wave, window)
        power /= n_trials
        self.assertTrue(sp.allclose(power/self.amp1**2, 1.0,
                                    atol=4.0*(2.0/sp.sqrt(n_trials))))
        # Expect this to fail ~1/100 times.
        self.assertFalse(sp.allclose(power/self.amp1**2, 1.0,
                                    atol=0.01*(2.0/sp.sqrt(n_trials))))

    def test_statistical_different_windows(self) :
        n_trials = 1000
        n_points = 200
        window1 = sp.ones(n_points, dtype=float)
        window1 += sp.sin(sp.arange(n_points)/10.0)
        window2 = 2*sp.ones(n_points, dtype=float)
        window2 += 2*sp.sin(sp.arange(n_points)/22.0)
        window2[window2<0.5] = 0
        window1[window1<0.5] = 0
        power = sp.zeros(n_points//2)
        for ii in range(n_trials) :
            wave = self.amp1*random.randn(n_points)
            power += np.windowed_power(wave*window1, window1, wave*window2,
                                       window2)
        power /= n_trials
        self.assertTrue(sp.allclose(power/self.amp1**2, 1.0,
                                    atol=6.0*(2.0/sp.sqrt(n_trials))))
        # Expect this to fail ~1/100 times.
        self.assertFalse(sp.allclose(power/self.amp1**2, 1.0,
                                    atol=0.01*(2.0/sp.sqrt(n_trials))))

class TestMakeMasked(unittest.TestCase) :
    
    def setUp(self) :
        self.Reader = fitsGBT.Reader(testfile, feedback=0)
        self.Blocks = self.Reader.read(None, 0)

    def test_runs(self) :
        data, mask, dt = np.make_masked_time_stream(self.Blocks)
        self.assertAlmostEqual(dt, 1.0, 2)
        self.assertTrue(data.shape[0] > 60)
        self.assertTrue(data.shape[0] < 120)

    def test_fails_overlapping(self) :
        Blocks = self.Reader.read(0, None)
        self.assertRaises(ce.DataError, np.make_masked_time_stream,
                          Blocks)

    def test_shortend(self) :
        data, mask, dt = np.make_masked_time_stream(self.Blocks, 65)
        self.assertEqual(data.shape[0], 65)
        ref_data, ref_mask, dt = np.make_masked_time_stream(self.Blocks)
        self.assertTrue(sp.allclose(ref_data[:65,...], data))



if __name__ == '__main__' :
    unittest.main()

