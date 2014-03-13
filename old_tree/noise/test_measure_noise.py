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
import scipy.signal as sig
import matplotlib.pyplot as plt

from core import data_block
import utils.misc as utils
import measure_noise as mn
import noise_power
from map.constants import T_infinity, T_small


class TestModule(unittest.TestCase) :
#class TestModule(object) :
    
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


class TestMeasures(unittest.TestCase) :
#class TestMeasures(object) :
    
    def setUp(self) :
        nf = 20
        nt = 1000
        dt = 0.13
        bw = 1.0/dt/2
        nb = 8
        self.nf = nf
        self.nt = nt
        self.dt = dt
        self.bw = bw
        self.nb = nb
        start_time = 87942.34
        self.time = sp.arange(nt*nb*1.6)*dt + start_time
        self.nt_total = len(self.time)
        self.data = rand.normal(size=(self.nt_total, 2, 1, nf))

    def make_blocks(self) :
        nf = self.nf
        nt = self.nt
        dt = self.dt
        bw = self.bw
        nb = self.nb
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
        nf = self.nf
        nt = self.nt
        dt = self.dt
        bw = self.bw
        nb = self.nb
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
            right_ans[3] = T_infinity**2
            self.assertTrue(sp.allclose(variance,
                                        right_ans, rtol=0.1))

    def test_gets_thermal_no_correlated(self):
        """Just check that the part of the freq_modes code that measures
        thermal is sane."""
        self.data *= sp.sqrt(self.bw * 2)
        Blocks = self.make_blocks()
        # Mask a channel out completly.
        for Data in Blocks:
            Data.data[:,:,:,3] = ma.masked
        parameters = mn.measure_noise_parameters(Blocks,
                                                ['freq_modes_over_f_0'])
        right_ans = sp.ones(self.nf)
        right_ans[3] = T_infinity**2
        for p in parameters.itervalues():
            thermal = p['freq_modes_over_f_0']['thermal']
            self.assertTrue(sp.allclose(thermal,
                                        right_ans, rtol=0.3))
            mean_thermal = sp.mean(thermal[right_ans==1])
            self.assertTrue(sp.allclose(mean_thermal, 1, rtol=0.05))

    def test_gets_thermal_with_correlated(self):
        """Checks that the part of the freq_modes code that compensates the
        thermal for mode subtraction works."""
        self.data *= sp.sqrt(self.bw * 2)  # Makes thermal unity.
        # Need to add something correlated so the modes arn't just the
        # channels.
        correlated_overf = noise_power.generate_overf_noise(1, 
                                    -2, 0.5, self.dt, self.data.shape[0])
        correlated_overf += (rand.normal(size=(self.data.shape[0],))
                             * sp.sqrt((self.bw * 2) * 0.3))
        self.data += correlated_overf[:,None,None,None] / sp.sqrt(self.nf)
        Blocks = self.make_blocks()
        # Mask a channel out completly.
        for Data in Blocks:
            Data.data[:,:,:,3] = ma.masked
        model = 'freq_modes_over_f_4'  # Take out 20% of the thermal power.
        parameters = mn.measure_noise_parameters(Blocks,
                                                [model])
        right_ans = sp.ones(self.nf)
        right_ans[3] = T_infinity**2
        for p in parameters.itervalues():
            pars = p[model]
            thermal = pars['thermal']
            self.assertTrue(sp.allclose(thermal,
                                        right_ans, rtol=0.3))
            mean_thermal = sp.mean(thermal[right_ans==1])
            self.assertTrue(sp.allclose(mean_thermal, 1, rtol=0.05))
            self.assertTrue(sp.allclose(pars['over_f_mode_0']['thermal'],
                                        0.3, atol=0.1))

    def test_gets_modes(self):
        nf = self.nf
        nt = self.nt
        dt = self.dt
        bw = self.bw
        nb = self.nb
        # Give every channel a different thermal noise floor.
        thermal_norm = 1.0 + 1.0/nf*sp.arange(nf)  # K**2/Hz
        self.data *= sp.sqrt(thermal_norm * bw * 2)
        n_time = self.data.shape[0]
        # Now make a 1/f like noise component in a few frequency modes.
        n_modes = 3
        index = -0.8 * (2.0 - sp.arange(n_modes, dtype=float)/n_modes)
        amp = 1.2 * (3.**(n_modes - 1.
                          - sp.arange(n_modes, dtype=float))) # K**2/Hz
        f_0 = 1.0 # Hz
        modes = sp.empty((n_modes, nf))
        for ii in range(n_modes):
            correlated_overf = noise_power.generate_overf_noise(amp[ii], 
                                        index[ii], f_0, dt, n_time)
            # Generate the frequency mode.  They should all be orthonormal.
            mode = sp.sin(2.*sp.pi*(ii + 1)*sp.arange(nf, dtype=float)/nf
                          + 6.4 * (ii + 3))
            mode /= sp.sqrt(sp.sum(mode**2))
            modes[ii] = mode
            self.data += correlated_overf[:,None,None,None] * mode
        # Add a subdominant general 1/f noise to all channels. Only do this for
        # stokes I, for V add nothing (to check that it measures 0).
        general_amp = 0.01
        general_index = -0.9
        general_cross_over = f_0 * general_amp**(-1./general_index)
        for ii in range(nf):
            tmp_a = general_amp * thermal_norm[ii]
            correlated_overf = noise_power.generate_overf_noise(tmp_a, 
                                        general_index, f_0, dt, n_time)
            self.data[:,0,:,ii] += correlated_overf[:,None]
        # Now put the data into the form of the real data.
        Blocks = self.make_blocks()
        # Measure all the noise parameters.
        model_name = 'freq_modes_over_f_' + str(n_modes)
        parameters = mn.measure_noise_parameters(Blocks, [model_name])
        for pol, pol_params in parameters.iteritems():
            
            for ii in range(n_modes):
                mode_noise = pol_params[model_name]['over_f_mode_' + str(ii)]
                self.assertTrue(sp.allclose(mode_noise['amplitude'], amp[ii],
                                            rtol=0.5))
                self.assertTrue(sp.allclose(mode_noise['index'], index[ii],
                                            atol=0.2))
                thermal_proj = sp.sum(thermal_norm * modes[ii,:]**2)
                #self.assertTrue(sp.allclose(mode_noise['thermal'], thermal_proj,
                #                            rtol=0.5))
                self.assertTrue(abs(sp.dot(mode_noise['mode'], modes[ii,:]))
                                > 0.90)
            thermal = pol_params[model_name]['thermal']
            loss = float(nf - n_modes) / nf
            self.assertTrue(sp.allclose(thermal, thermal_norm * loss,
                                        rtol=0.4))
            measured_general_ind = pol_params[model_name]['all_channel_index']
            measured_corner = pol_params[model_name]['all_channel_corner_f']
            if pol == 1:
                self.assertTrue(sp.allclose(measured_general_ind,
                                            general_index, atol=0.4))
                #Only need logarithmic accuracy on the corner.
                self.assertTrue(sp.allclose(sp.log(measured_corner),
                                            sp.log(general_cross_over),
                                            atol=1.5))
            elif pol == 3:
                self.assertTrue(measured_corner < 4. / dt / nt / nb)

    def test_gets_modes_by_scan(self):
        nf = self.nf
        nt = self.nt
        dt = self.dt
        bw = self.bw
        nb = self.nb
        # Give every channel a different thermal noise floor.
        thermal_norm = 1.0 + 1.0/nf*sp.arange(nf)  # K**2/Hz
        self.data *= sp.sqrt(thermal_norm * bw * 2)
        n_time = self.data.shape[0]
        # Now make a 1/f like noise component in a few frequency modes.
        n_modes = 3
        index = -0.8 * (2.0 - sp.arange(n_modes, dtype=float)/n_modes)
        amp = 1.2 * (3.**(n_modes - 1.
                          - sp.arange(n_modes, dtype=float))) # K**2/Hz
        f_0 = 1.0 # Hz
        modes = sp.empty((n_modes, nf))
        for ii in range(n_modes):
            correlated_overf = noise_power.generate_overf_noise(amp[ii], 
                                        index[ii], f_0, dt, n_time)
            # Generate the frequency mode.  They should all be orthonormal.
            mode = sp.sin(2.*sp.pi*(ii + 1)*sp.arange(nf, dtype=float)/nf
                          + 6.4 * (ii + 3))
            mode /= sp.sqrt(sp.sum(mode**2))
            modes[ii] = mode
            self.data += correlated_overf[:,None,None,None] * mode
        # Add a subdominant general 1/f noise to all channels.
        general_amp = 0.1
        general_index = -0.9
        general_cross_over = f_0 * general_amp**(-1./general_index)
        for ii in range(nf):
            tmp_a = general_amp * thermal_norm[ii]
            correlated_overf = noise_power.generate_overf_noise(tmp_a, 
                                        general_index, f_0, dt, n_time)
            self.data[:,0,:,ii] += correlated_overf[:,None]
        # Now put the data into the form of the real data.
        Blocks = self.make_blocks()
        # Measure all the noise parameters.
        model_name = 'freq_modes_over_f_' + str(n_modes)
        parameters = mn.measure_noise_parameters(Blocks, [model_name],
                                                 split_scans=True)
        for pol, pol_params in parameters.iteritems():
            
            for ii in range(n_modes):
                mode_noise = pol_params[model_name]['over_f_mode_' + str(ii)]
                self.assertTrue(sp.allclose(mode_noise['amplitude'], amp[ii],
                                            rtol=0.5))
                self.assertTrue(sp.allclose(mode_noise['index'], index[ii],
                                            atol=0.2))
                #thermal_proj = sp.sum(thermal_norm * modes[ii,:]**2)
                #self.assertTrue(sp.allclose(mode_noise['thermal'], thermal_proj,
                #                            rtol=0.5))
                self.assertTrue(abs(sp.dot(mode_noise['mode'], modes[ii,:]))
                                > 0.90)
            thermal = pol_params[model_name]['thermal']
            loss = float(nf - n_modes) / nf
            self.assertTrue(sp.allclose(thermal, thermal_norm * loss,
                                        rtol=0.4))
            measured_general_ind = pol_params[model_name]['all_channel_index']
            measured_corner = pol_params[model_name]['all_channel_corner_f']
            if pol == 1:
                self.assertTrue(sp.allclose(measured_general_ind,
                                            general_index, atol=0.4))
                #Only need logarithmic accuracy on the corner.
                self.assertTrue(sp.allclose(sp.log(measured_corner),
                                            sp.log(general_cross_over),
                                            atol=1.5))
            elif pol == 3:
                self.assertTrue(measured_corner < 4. / dt / nt)

    def test_all_masked(self):
        Blocks = self.make_blocks()
        for Data in Blocks:
            Data.data[...] = ma.masked
        model_name = 'freq_modes_over_f_' + str(3)
        parameters = mn.measure_noise_parameters(Blocks, [model_name])
        for pol, pol_params in parameters.iteritems():
            for ii in range(3):
                mode_noise = pol_params[model_name]['over_f_mode_' + str(ii)]
                self.assertTrue(sp.allclose(mode_noise['amplitude'], 0))
                self.assertTrue(sp.allclose(mode_noise['index'], 0))
                self.assertTrue(sp.allclose(mode_noise['thermal'],
                                            T_infinity**2))
            thermal = pol_params[model_name]['thermal']
            self.assertTrue(sp.allclose(thermal, T_infinity**2))
            measured_general_ind = pol_params[model_name]['all_channel_index']
            measured_corner = pol_params[model_name]['all_channel_corner_f']
            self.assertTrue(sp.allclose(measured_general_ind, 0))
            self.assertTrue(sp.allclose(measured_corner, 0))

    def test_few_masked(self):
        Blocks = self.make_blocks()
        for Data in Blocks:
            Data.data[:,:,:,6] = ma.masked
            Data.data[:,:,:,13] = ma.masked
        model_name = 'freq_modes_over_f_' + str(3)
        parameters = mn.measure_noise_parameters(Blocks, [model_name])
        for pol, pol_params in parameters.iteritems():
            for ii in range(3):
                mode_noise = pol_params[model_name]['over_f_mode_' + str(ii)]
                #self.assertTrue(sp.allclose(mode_noise['thermal'], self.dt,
                #                            rtol=0.5))
                self.assertTrue(sp.allclose(mode_noise['mode'][6], 0,
                                            atol=1.e-10))
                self.assertTrue(sp.allclose(mode_noise['mode'][13], 0,
                                            atol=1.e-10))
            expected = sp.ones(self.nf, dtype=float) * self.dt
            expected[6] = T_infinity**2
            expected[13] = T_infinity**2
            thermal = pol_params[model_name]['thermal']
            # weak test since the above modes may have favoured a few channels.
            #print thermal, expected
            self.assertTrue(sp.allclose(thermal, expected, rtol=0.8))
            measured_general_ind = pol_params[model_name]['all_channel_index']
            measured_corner = pol_params[model_name]['all_channel_corner_f']
            self.assertTrue(measured_corner < 4. / self.dt / self.nt / self.nb)


class TestFunctions(unittest.TestCase):

    def test_fit_over_f_plus_const(self):
        dt = 0.13
        n_time = 10000
        amp = 0.67 # K**2/Hz
        index = -1.3
        f_0 = 1.0
        thermal = 2.7 # K**2/Hz
        BW = 1./dt/2
        window = sig.get_window('hanning', n_time)
        n_spec = 10
        p = 0
        for ii in range(n_spec):
            time_stream = noise_power.generate_overf_noise(amp, index, f_0,
                                                                dt, n_time)
            time_stream += rand.normal(size=n_time) * sp.sqrt(thermal * BW * 2)
            time_stream -= sp.mean(time_stream)
            time_stream *= window
            p += noise_power.calculate_power(time_stream)
        p /= n_spec
        p = noise_power.make_power_physical_units(p, dt)
        w = noise_power.calculate_power(window)
        w_norm = sp.mean(w).real
        #w /= w_norm
        p = noise_power.prune_power(p).real
        #p /= w_norm
        f = noise_power.ps_freq_axis(dt, n_time)
        p = p[1:]
        f = f[1:]
        amp_m, index_m, f0_m, thermal_m = mn.fit_overf_const(p, w, f)
        self.assertTrue(sp.allclose(amp_m, amp, atol=0.2))
        self.assertTrue(sp.allclose(index_m, index, atol=0.1))
        self.assertTrue(sp.allclose(thermal_m, thermal, atol=0.1))


if __name__ == '__main__' :
    unittest.main()
