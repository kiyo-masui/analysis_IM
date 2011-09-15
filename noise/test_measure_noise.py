"""Unit tests for reflag.py"""

import unittest
import os
import glob

import scipy as sp
import scipy.signal as sig
import numpy.random as rand
import scipy.fftpack as fft

from core import data_block, utils
import measure_noise as mn
import noise_power

class TestModule(unittest.TestCase) :
    
    def test_module(self) :
        params = {"mn_output_filename" : "noise_params.shelf",
                  "mn_output_root" : "./testout_"
                 }
        mn.Measure(params, feedback=0).execute()
        files = glob.glob('*testout*')
        self.assertTrue(len(files) > 1)

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)

nf = 5
nt = 1000
dt = 0.13
bw = 1.0/dt/2
nb = 5
class TestMeasures(unittest.TestCase) :
    
    def setUp(self) :
        start_time = 87942.34
        self.time = sp.arange(nt*nb*1.5)*dt + start_time
        self.nt_total = len(self.time)
        self.data = rand.normal(size=(self.nt_total, 1, 1, nf))

    def make_blocks(self) :
        Blocks = []
        for ii in range(nb):
            start_ind = ii*1.4*nt
            data = self.data[start_ind:start_ind + nt, ...]
            Data = data_block.DataBlock(data)
            time_sec = self.time[start_ind:start_ind + nt]
            time_strings = utils.float2time(time_sec)
            Data.set_field("DATE-OBS", time_strings, axis_names=('time',))
            Blocks.append(Data)
        return Blocks
    
    def test_get_variance(self) :
        self.data *= sp.arange(1, nf+1)
        Blocks = self.make_blocks()
        self.assertTrue(sp.allclose(mn.get_var(Blocks)
                                    / sp.arange(1, nf+1)**2, 1.0, atol=0.1))

    def test_get_overf(self) :
        # Give every channel a different thermal noise floor.
        thermal_norm = 1.0 + 1.0/nf*sp.arange(nf)  # K**2/Hz
        self.data *= sp.sqrt(thermal_norm * bw)
        # Now make a 1/f like noise component
        ntime = self.data.shape[0]
        index = -1.45
        amp = 2.3  # K**2/Hz
        f_0 = 0.5 #  Hz
        correlated_overf = noise_power.generate_overf_noise(amp, index, f_0,
                                                            dt, ntime)
        # Add 1/f to all channels (perfectly correlated).
        self.data += correlated_overf[:,None,None,None]
        Blocks = self.make_blocks()
        # Measure the 1/f and is if we get back what we put in.
        thermal, overf = mn.get_correlated_overf(Blocks, f_0)
        self.assertTrue(sp.allclose(thermal, thermal_norm, rtol=0.15))
        self.assertTrue(sp.allclose(overf[0], amp, rtol=0.3))
        self.assertTrue(sp.allclose(overf[1], index, atol=0.1))

if __name__ == '__main__' :
    unittest.main()
