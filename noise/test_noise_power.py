"""Unit tests for noise_power module."""

import unittest

import scipy as sp
import scipy.fftpack as fft
from numpy import random
import scipy.linalg as linalg
import matplotlib.pyplot as plt

import noise_power as npow
from core import fitsGBT
import kiyopy.custom_exceptions as ce

testfile = "testdata/testfile_GBTfits.fits"

class TestWindowedPower(unittest.TestCase) :

    def setUp(self) :
        self.n = 10000
        self.mode1 = 677 # Choose the 677th mode.
        self.k1 = (2.0*sp.pi*self.mode1)/self.n
        self.amp1 = 3.45
        self.wave1 = self.amp1*sp.sin(self.k1*(sp.arange(self.n) + 6.45))

    def test_no_window(self) :
        window = sp.ones_like(self.wave1)
        power = npow.windowed_power(self.wave1, window)
        power = npow.prune_power(power)
        self.assertAlmostEqual(power[self.mode1]/self.amp1**2/self.n*4, 1)
        self.assertTrue(sp.allclose(power[:self.mode1], 0, 
                                    atol=self.amp1**2*self.n/1e15))
        self.assertTrue(sp.allclose(power[self.mode1+1:], 0, 
                                    atol=self.amp1**2*self.n/1e15))
        # With no window, we have a quick way to the answer
        quick_power = npow.calculate_power(self.wave1)
        self.assertTrue(sp.allclose(power, quick_power[:self.n//2]))

    def test_convolve_normalization(self) :
        window = sp.ones_like(self.wave1)
        power = npow.calculate_power(self.wave1)
        window_power = npow.calculate_power(window)
        # Convolving with the window function should do nothing for all ones.
        convolved_power = npow.convolve_power(power, window_power)
        self.assertTrue(sp.allclose(power, convolved_power))

    def test_convolve(self) :
        window = sp.ones_like(self.wave1)
        window += sp.sin(sp.arange(self.n)/50.0)
        self.wave1 *= window
        power = npow.calculate_power(self.wave1)
        window_power = npow.calculate_power(window)
        deconvolved_power = npow.deconvolve_power(power, window_power)
        reconvolved_power = npow.convolve_power(deconvolved_power, window_power)
        self.assertTrue(sp.allclose(power, reconvolved_power))

    def test_non_zero_window(self) :
        window = sp.ones_like(self.wave1)
        window += 0.5*sp.sin(sp.arange(self.n)/15.0)
        self.wave1 *= window
        power = npow.windowed_power(self.wave1, window)
        power = npow.prune_power(power)
        self.assertAlmostEqual(power[self.mode1]/self.amp1**2/self.n*4, 1.0, 3)
        self.assertTrue(sp.allclose(power[:self.mode1], 0, 
                                    atol=self.amp1**2*self.n/1e3))
        self.assertTrue(sp.allclose(power[self.mode1+1:], 0, 
                                    atol=self.amp1**2*self.n/1e3))

    def test_window_with_zeros(self) :
        window = sp.ones_like(self.wave1)
        window += sp.sin(sp.arange(self.n)/50.0)
        self.wave1 *= window
        power = npow.windowed_power(self.wave1, window)
        power = npow.prune_power(power)
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
        power = npow.windowed_power(wave1, window1, wave2, window2)
        power = npow.prune_power(power)
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
            p = npow.windowed_power(wave, window)
            power += npow.prune_power(p).real
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
            p = npow.windowed_power(wave*window1, window1, wave*window2,
                                       window2)
            power += npow.prune_power(p).real
        power /= n_trials
        self.assertTrue(sp.allclose(power/self.amp1**2, 1.0,
                                    atol=6.0*(2.0/sp.sqrt(n_trials))))
        # Expect this to fail ~1/100 times.
        self.assertFalse(sp.allclose(power/self.amp1**2, 1.0,
                                    atol=0.01*(2.0/sp.sqrt(n_trials))))

    def test_statistical_physical_units(self) :
        n_trials = 1000
        n_points = 200
        dt = 0.001
        window = sp.ones(n_points, dtype=float)
        power = sp.zeros(n_points//2)
        for ii in range(n_trials) :
            wave = self.amp1*random.randn(n_points)
            power += npow.prune_power(npow.calculate_power(wave)).real
        power /= n_trials
        power = npow.make_power_physical_units(power, dt)
        freqs = npow.ps_freq_axis(dt, n_points)
        df = abs(sp.mean(sp.diff(freqs)))
        # The integral of the power spectrum should be the variance. Factor of
        # 2 get the negitive frequencies.
        integrated_power = sp.sum(power) * df * 2
        self.assertTrue(sp.allclose(integrated_power/self.amp1**2, 1.0,
                                    atol=4.0*(2.0/sp.sqrt(n_trials*n_points))))

class TestMakeMasked(unittest.TestCase) :
    
    def setUp(self) :
        self.Reader = fitsGBT.Reader(testfile, feedback=0)
        self.Blocks = self.Reader.read(None, 0)

    def test_runs(self) :
        data, mask, dt = npow.make_masked_time_stream(self.Blocks)
        self.assertAlmostEqual(dt, 1.0, 2)
        self.assertTrue(data.shape[0] > 60)
        self.assertTrue(data.shape[0] < 120)

    def test_fails_overlapping(self) :
        Blocks = self.Reader.read(0, None)
        self.assertRaises(ce.DataError, npow.make_masked_time_stream,
                          Blocks)

    def test_shortend(self) :
        data, mask, dt = npow.make_masked_time_stream(self.Blocks, 65)
        self.assertEqual(data.shape[0], 65)
        ref_data, ref_mask, dt = npow.make_masked_time_stream(self.Blocks)
        self.assertTrue(sp.allclose(ref_data[:65,...], data))

class TestUtils(unittest.TestCase):
    
    def test_overf_correlation_white(self):
        n = 1000
        dt = 0.01
        BW = 1.0/2/dt
        amplitude = 0.5 / BW / 2   # K^2/Hz
        corr = npow.calculate_overf_correlation(amplitude, 0, 1, dt, n)
        self.assertAlmostEqual(corr[0], amplitude*BW*2)
        self.assertTrue(sp.allclose(corr[1:], 0))

    def test_overf_correlation_bw(self):
        n1 = 50000
        dt1 = 0.001
        corr1 = npow.calculate_overf_correlation(10, -1.3, 1, dt1, n1)
        n2 = 5000
        dt2 = 0.01
        corr2 = npow.calculate_overf_correlation(10, -1.3, 1, dt2, n2)
        # Asside from finer sampling, only the short time scale stuff should
        # differ.
        self.assertTrue(sp.allclose(sp.diff(corr1[0::10])[10:], 
                                    sp.diff(corr2)[10:], rtol=1e-2))
    
    def test_overf_correlation_total_time(self):
        n1 = 50000
        dt1 = 0.01
        BW1 = 1.0/2/dt1
        corr1 = npow.calculate_overf_correlation(10.0 / BW1 / 2, -1.3, 1.0, 
                                                      dt1, n1)
        n2 = 5000
        dt2 = 0.01
        BW2 = 1.0/2/dt2
        corr2 = npow.calculate_overf_correlation(10.0 / BW2 / 2, -1.3, 1.0,
                                                    dt2, n2)
        # The overlapping time scale parts should all be about the same.
        self.assertTrue(sp.allclose(sp.diff(corr1)[:n2//2], 
                                    sp.diff(corr2)[:n2//2], rtol=1e-2))
    
    def test_overf_correlation_statistical(self):
        # Parameters
        n = 10000
        n_corr = 50
        n_trials = 100
        dt = 0.64  # s
        f0 = 0.01  # Hz
        index = -1.3
        amp = 2.34  # K**2/Hz
        BW = 1.0/2/dt  # Hz
        corr = sp.zeros(n_corr, dtype=float)
        for ii in range(n_trials):
            # Generate noise.
            noise_data = npow.generate_overf_noise(amp, index, f0, dt, n)
            # Calculate the correlation function.
            corr[0] += sp.mean(noise_data **2)
            for jj in range(1, n_corr):
                this_corr = sp.mean(noise_data[:-jj] * noise_data[jj:])
                corr[jj] += this_corr
        corr /= n_trials
        corr_theory = npow.calculate_overf_correlation(amp, index, f0, dt,
                                                       n)[:n_corr]
        self.assertTrue(sp.allclose(sp.diff(corr), sp.diff(corr_theory),
                                    rtol=0.1))




if __name__ == '__main__' :
    unittest.main()

