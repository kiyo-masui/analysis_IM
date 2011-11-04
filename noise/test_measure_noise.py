"""Unit tests for reflag.py"""

import unittest
import os
import glob
import shelve

import scipy as sp
import scipy.signal as sig
import numpy.random as rand
import scipy.fftpack as fft
import numpy.ma as ma

from core import data_block, utils
import measure_noise as mn
import noise_power

class TestModule(unittest.TestCase) :
    
    def test_module(self) :
        params = {"mn_output_filename" : "noise_params.shelve",
                  "mn_output_root" : "./testout_",
                  "mn_file_middles" : ("testfile_guppi_combined",)
                 }
        mn.Measure(params, feedback=0).execute()
        files = glob.glob('*testout*')
        self.assertTrue(len(files) > 1)
        pdata = shelve.open("./testout_noise_params.shelve")
        

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)

nf = 5
nt = 1000
dt = 0.13
bw = 1.0/dt/2
nb = 8
class TestMeasures(unittest.TestCase) :
    
    def setUp(self) :
        start_time = 87942.34
        self.time = sp.arange(nt*nb*1.6)*dt + start_time
        self.nt_total = len(self.time)
        self.data = rand.normal(size=(self.nt_total, 2, 1, nf))

    def make_blocks(self) :
        Blocks = []
        for ii in range(nb):
            start_ind = ii*1.4*nt
            data = self.data[start_ind:start_ind + nt, ...]
            Data = data_block.DataBlock(data)
            time_sec = self.time[start_ind:start_ind + nt]
            time_strings = utils.float2time(time_sec)
            Data.set_field("DATE-OBS", time_strings, axis_names=('time',))
            pols = [1, 3,]
            Data.set_field("CRVAL4", pols, axis_names=('pol',))
            cals = ['T']
            Data.set_field("CAL", cals, axis_names=('cal',))
            Blocks.append(Data)
        return Blocks
    
    def test_get_variance(self) :
        self.data *= sp.arange(1, nf+1)
        Blocks = self.make_blocks()
        # Mask a channel out completly.
        for Data in Blocks:
            Data.data[:,:,:,3] = ma.masked
        parameters = mn.measure_noise_parameters(Blocks, ['channel_var'])
        for ii in [1,3]:
            self.assertTrue(parameters.has_key(ii))
        for p in parameters.itervalues():
            variance = p['channel_var']
            right_ans = sp.arange(1, nf+1)**2
            right_ans[3] = 0
            self.assertTrue(sp.allclose(variance,
                                        right_ans, rtol=0.1))

    def test_get_overf(self) :
        # Give every channel a different thermal noise floor.
        thermal_norm = 1.0 + 1.0/nf*sp.arange(nf)  # K**2/Hz
        self.data *= sp.sqrt(thermal_norm * bw)
        # Now make a 1/f like noise component
        ntime = self.data.shape[0]
        index = -1.45
        amp = 0.9  # K**2/Hz
        f_0 = 1.0 #  Hz
        correlated_overf = noise_power.generate_overf_noise(amp, index, f_0,
                                                            dt, ntime)
        # Add 1/f to all channels (perfectly correlated).
        self.data += correlated_overf[:,None,None,None]
        Blocks = self.make_blocks()
        # Mask a channel.
        for Data in Blocks:
            Data.data[:,:,:,1] = ma.masked
        thermal_norm_masked = thermal_norm.copy()
        thermal_norm_masked[1] = 0
        # Measure the 1/f and see if we get back what we put in.
        parameters = mn.measure_noise_parameters(Blocks, ['mean_over_f'])
        for ii in [1,3]:
            self.assertTrue(parameters.has_key(ii))
        for p in parameters.itervalues():
            over_f_params = p['mean_over_f'][1]
            thermal = p['mean_over_f'][0]
            self.assertTrue(sp.allclose(thermal, thermal_norm_masked, 
                                        rtol=0.15, atol=0.01))
            self.assertTrue(sp.allclose(over_f_params[0], amp, rtol=0.3))
            self.assertTrue(sp.allclose(over_f_params[1], index, atol=0.1))

if __name__ == '__main__' :
    unittest.main()
