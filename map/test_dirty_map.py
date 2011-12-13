"""Unit tests for dirty map maker."""

import os
import glob
import time as time_module

import unittest
import scipy as sp
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import linalg
from scipy import interpolate
from numpy import random

import core.algebra as al
import dirty_map as dirty_map
from noise import noise_power
from core import fitsGBT

# Default data parameters used for most of the tests.
nscans_d = 6
nt_scan_d = 50
nt_d = nt_scan_d * nscans_d
nf_d = 5
scan_size_d = 5.0

nra_d = 20
ndec_d = 20
map_size_d = 6.0

dt_d = 0.1
BW_d = 1./2./dt_d


class DataMaker(object):
    """Class that generates test data."""

    def __init__(self, nscans=nscans_d, nt_scan=nt_scan_d, nf=nf_d, nra=nra_d,
                 ndec=ndec_d, scan_size=scan_size_d, map_size=map_size_d,
                 thermal=0, correlated_noise=None, random_mean=0,
                 random_slope=0, add_ground=False, freq_mode_noise=None,
                 universal_over_f=None):
        # Store parameters.
        self.nscans = nscans
        self.nt_scan = nt_scan
        self.nf = nf
        self.nra = nra
        self.ndec = ndec
        self.scan_size = float(scan_size)
        self.map_size = float(map_size)
        self.add_ground = add_ground
        # Noise parameters.
        self.thermal = thermal
        self.correlated_noise = correlated_noise
        self.freq_mode_noise = freq_mode_noise
        self.random_mean = random_mean
        self.random_slope = random_slope
        self.universal_over_f = universal_over_f
        # Initialize a counter.
        self.data_set_number = 0
        self.nt = nt_scan * nscans
        self.dt = dt_d

    def ground_az_el(self, az, el):
        time_stream = sp.zeros((self.nf, self.nt))
        time_stream += sp.sin(az*0.3)
        time_stream += sp.cos(el*0.2)
        time_stream *= -sp.arange(self.nf, dtype=float)[:,None]/self.nf + 2.0
        return time_stream

    def sky_ra_dec(self, ra, dec):
        ra_dec_shape = (ra * dec).shape
        sl = (slice(None),) + (None,) * (ra * dec).ndim
        time_stream = sp.zeros((self.nf,) + ra_dec_shape)
        time_stream += sp.sin(ra*0.9 + 3.0)*sp.cos(dec*1.3 + 2.0) # The sky map.
        time_stream *= sp.cos(sp.arange(self.nf, dtype=float) + 5.0)[sl] + 1.2
        return time_stream

    def get_time(self, number):
        time = (sp.arange(self.nt_scan) + sp.arange(self.nscans)[:,None] 
                * 1.5 * self.nt_scan)
        time = time.flat[...] + 2.0 * self.nt*number
        time = time * self.dt + 123456.0
        return time

    def _get_scanning(self):
        nt = self.nt
        nt_scan = self.nt_scan
        nscans = self.nscans
        scanning = sp.zeros(self.nt)
        for ii in range(0, nscans, 2):
            scanning[ii*nt_scan:(ii+1)*nt_scan] = sp.arange(-0.5*nt_scan,
                                                            0.5*nt_scan, 
                                                            dtype=float)
        for ii in range(1, nscans, 2):
            scanning[ii*nt_scan:(ii+1)*nt_scan] = -sp.arange(-0.5*nt_scan,
                                                             0.5*nt_scan, 
                                                             dtype=float)
        scanning *= self.scan_size / nt_scan
        return scanning

    def get_az_el(self, number):
        # Elevation always constant with current scan strategy.
        el = sp.zeros(self.nt) + (13.13 * number) % 40 + 10
        # Az scans back and forth.
        az = self._get_scanning()
        az /= sp.cos(el[0]*sp.pi/180.0)
        az += (26.64 * number) % 290

        return az, el

    def get_ra_dec(self, number):
        # Position angle of the telescope controls how much ra vs dec we scan.
        # Number between -70 and 70 degrees.
        # Using pi as an irrational number so we don't get repeats.
        pa = -70 + sp.pi*2*(number+1) % 140
        # The angles in the drift direction.
        drifting = (sp.arange(-self.nt/2.0, self.nt/2.0, dtype=float)
                    / self.nt * self.scan_size)
        scanning =  self._get_scanning()
        # Use `pa` to rotate these to ra-dec.
        ra = sp.sin(pa*sp.pi/180)*drifting + sp.cos(pa*sp.pi/180)*scanning
        dec = sp.cos(pa*sp.pi/180)*drifting - sp.sin(pa*sp.pi/180)*scanning
        # Add maps centres and the some offset jitter.
        ra += 21 + (((sp.pi * (number+3)) % 1.0) - 0.5) * 0.3
        dec += 0 + (((sp.pi * (number+7)) % 1.0) - 0.5) * 0.3
        return ra, dec

    def get_map(self):
        map = sp.zeros((self.nf, self.nra, self.ndec))
        map = al.make_vect(map, ('freq', 'ra', 'dec'))
        map.set_axis_info('freq', 800e6, 1e6)
        map.set_axis_info('ra', 21, self.map_size/self.nra)
        map.set_axis_info('dec', 0, self.map_size/self.ndec)
        return map

    def get_noise(self, number):
        thermal = sp.zeros(self.nf) + self.thermal
        noise_data = random.randn(self.nf, self.nt) * sp.sqrt(thermal)[:,None]
        t = self.get_time(number)
        # Noisy mean mode.
        noise_data += random.randn(self.nf, 1) * sp.sqrt(self.random_mean)
        # Noisy slope mode.
        slope_noise = random.randn(self.nf, 1) * sp.sqrt(self.random_slope)
        slope_mode = (sp.sqrt(self.nt) * (t - sp.mean(t)) 
                       / sp.sqrt(sp.sum((t - sp.mean(t))**2)))
        noise_data += slope_noise * slope_mode
        if not self.universal_over_f is None:
            amps = sp.empty(self.nf)
            amps[:] = self.thermal * 2.0 * self.dt
            index = self.universal_over_f[0]
            f0 = self.universal_over_f[1]
            for ii in range(self.nf):
                over_f = noise_power.generate_overf_noise(amps[ii], index, f0,
                                                          self.dt, self.nt*2)
                ov_f_interpolator = interpolate.interp1d(t[0] 
                                    + sp.arange(self.nt*2)*self.dt, over_f)
                noise_data[ii,:] += ov_f_interpolator(t)
        if not self.correlated_noise is None:
            amp = self.correlated_noise[0]
            index = self.correlated_noise[1]
            f0 = self.correlated_noise[2]
            over_f = noise_power.generate_overf_noise(amp, index, f0,
                                                      self.dt, self.nt*2)
            ov_f_interpolator = interpolate.interp1d(t[0] 
                                    + sp.arange(self.nt*2)*self.dt, over_f)
            noise_data += ov_f_interpolator(t)
        if not self.freq_mode_noise is None:
            for ii, p in enumerate(self.freq_mode_noise):
                mode = make_frequency_mode(ii, self.nf)
                amp = p[0]
                index = p[1]
                f0 = p[2]
                thermal = p[3]
                over_f = noise_power.generate_overf_noise(amp, index, f0,
                                                      self.dt, self.nt*2)
                ov_f_interpolator = interpolate.interp1d(t[0] 
                                    + sp.arange(self.nt*2)*self.dt, over_f)
                mode_noise = ov_f_interpolator(t)
                thermal_amp = sp.sqrt(thermal / 2. / self.dt)
                mode_noise += random.randn(self.nt) * thermal_amp
                noise_data += mode[:,None] * mode_noise
        return noise_data

    def get_all(self, number=None):
        if number == None:
            number = self.data_set_number
            self.data_set_number += 1
        time = self.get_time(number)
        az, el = self.get_az_el(number)
        ra, dec = self.get_ra_dec(number)
        time_stream = sp.zeros((self.nf, self.nt))
        time_stream += self.get_noise(number)
        time_stream += self.sky_ra_dec(ra, dec)
        if self.add_ground:
            time_stream += self.ground_az_el(az, el)
        time_stream = al.make_vect(time_stream, ('freq', 'time'))
        # Make of some mask_inds.
        mask_inds = (sp.empty(4, dtype=int), sp.empty(4, dtype=int))
        for jj in range(4):
            mask_inds[0][jj] = (13 * ((number + 1)*(jj + 2))) % self.nf
            mask_inds[1][jj] = (53 * ((number + 3)*(jj + 4))) % self.nt

        return time_stream, ra, dec, az, el, time, mask_inds

    def get_all_trimmed(self, number=None):
        time_stream, ra, dec, az, el, time, mask_inds = self.get_all(number)
        map = self.get_map()
        time_stream, inds = dirty_map.trim_time_stream(time_stream, (ra, dec),
                (min(map.get_axis('ra')), min(map.get_axis('dec'))),
                (max(map.get_axis('ra')), max(map.get_axis('dec'))))
        ra = ra[inds]
        dec = dec[inds]
        time = time[inds]
        az = az[inds]
        el = el[inds]
        new_mask_inds = ([], [])
        for jj in range(len(mask_inds[0])):
            if mask_inds[1][jj] in inds:
                new_time_ind = sp.where(inds == mask_inds[1][jj])[0][0]
                new_mask_inds[0].append(mask_inds[0][jj])
                new_mask_inds[1].append(new_time_ind)
        mask_inds = (sp.asarray(new_mask_inds[0], dtype=int),
                     sp.asarray(new_mask_inds[1], dtype=int))
        return time_stream, ra, dec, az, el, time, mask_inds


def make_frequency_mode(mode_num, nf):
    """Deterministically return a normal vector of length `nf`.

    Different modes are orthoganal for `mode_num` < `nf`/2.
    """

    k = 2. * sp.pi * float(mode_num + 1) / nf
    phase = (mode_num + 6) * sp.sqrt(2)
    mode = sp.cos(k * sp.arange(nf) + phase)
    mode /= sp.sqrt(sp.dot(mode, mode))
    return mode


class TestClasses(unittest.TestCase):

    def setUp(self):
        self.DM = DataMaker()
        self.map = self.DM.get_map()
    
    def test_get_matrix(self):
        time_stream, ra, dec, az, el, time, mask_inds = \
                self.DM.get_all_trimmed()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), self.map, 'nearest')
        pointing_matrix = P.get_matrix()
        self.assertTrue(pointing_matrix.axes == ('time', 'ra', 'dec'))
        self.assertEqual(pointing_matrix.shape[0], ra.shape[0])
        self.assertEqual(self.map.shape[1], pointing_matrix.shape[1])
        self.assertEqual(self.map.shape[2], pointing_matrix.shape[2])

    def test_apply_time_axis(self):
        # The function that this is testing hasn't been implemented. Test
        # disabled.
        time_stream, ra, dec, az, el, time, mask_inds = \
                self.DM.get_all_trimmed()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), self.map, 'linear')
        pointing_matrix = P.get_matrix()
        gridded_data_mult = al.partial_dot(time_stream,
                                           pointing_matrix)
        gridded_data_fast = P.apply_to_time_axis(time_stream)
        self.assertTrue(sp.allclose(gridded_data_mult, gridded_data_fast))

    def test_build_noise(self):
        map = self.map
        time_stream, ra, dec, az, el, time, mask_inds = \
                                               self.DM.get_all_trimmed()
        nt = len(time)
        Noise = dirty_map.Noise(time_stream, time)
        thermal_noise_levels = sp.zeros((nf_d)) + 0.04  # Kelvin**2
        Noise.add_thermal(thermal_noise_levels)
        Noise.add_mask(mask_inds)
        self.assertTrue(sp.alltrue(Noise.diagonal[mask_inds] > 10))
        Noise.deweight_time_mean()
        Noise.deweight_time_slope()
        Noise.add_correlated_over_f(0.01, -1.2, 0.1)
        Noise.finalize()
        #### Test the full inverse.
        # Frist get a full representation of the noise matrix
        tmp_mat = sp.zeros((nf_d, nt, nf_d, nt))
        tmp_mat.flat[::nt*nf_d + 1] += Noise.diagonal.flat
        for jj in xrange(Noise.time_modes.shape[0]):
            tmp_mat += (Noise.time_mode_noise[jj,:,None,:,None]
                        * Noise.time_modes[jj,None,:,None,None]
                        * Noise.time_modes[jj,None,None,None,:])
        for jj in xrange(Noise.freq_modes.shape[0]):
            tmp_mat +=  (Noise.freq_mode_noise[jj,None,:,None,:]
                         * Noise.freq_modes[jj,:,None,None,None]
                         * Noise.freq_modes[jj,None,None,:,None])
        tmp_mat.shape = (nt*nf_d, nt*nf_d)
        # Check that the matrix I built for testing is indeed symetric.
        self.assertTrue(sp.allclose(tmp_mat, tmp_mat.transpose()))
        noise_inv = Noise.get_inverse()
        noise_inv.shape = (nt*nf_d, nt*nf_d)
        # Check that the production matrix is symetric.
        self.assertTrue(sp.allclose(noise_inv, noise_inv.transpose()))
        tmp_eye = sp.dot(tmp_mat, noise_inv)
        #print tmp_eye
        noise_inv.shape = (nf_d, nt, nf_d, nt)
        self.assertTrue(sp.allclose(tmp_eye, sp.identity(nt*nf_d)))
        #### Test the noise weighting of the data.
        noise_weighted_data = Noise.weight_time_stream(time_stream)
        self.assertTrue(sp.allclose(noise_weighted_data, al.dot(noise_inv,
                                                                time_stream)))
        #### Test making noise in map space.
        # First make the noise matrix by brute force.
        P = dirty_map.Pointing(("ra", "dec"), (ra, dec), map, 'nearest')
        P_mat = P.get_matrix()
        tmp_map_noise_inv = al.partial_dot(noise_inv,
                                           P_mat)
        tmp_map_noise_inv = al.partial_dot(P_mat.mat_transpose(), 
                                           tmp_map_noise_inv)
        # I mess up the meta data by doing this, but rotate the axes so they
        # are in the desired order.
        tmp_map_noise_inv = sp.rollaxis(tmp_map_noise_inv, 2, 0)
        # Now use fast methods.
        map_noise_inv = sp.zeros((nf_d, nra_d, ndec_d, nf_d, nra_d, ndec_d),
                                 dtype=float)
        map_noise_inv = al.make_mat(map_noise_inv, axis_names=('freq', 'ra', 
            'dec', 'freq', 'ra', 'dec'), row_axes=(0, 1, 2), 
            col_axes=(3, 4, 5))
        start = time_module.clock()
        for ii in xrange(nf_d):
            for jj in xrange(nra_d):
                P.noise_to_map_domain(Noise, ii, jj,
                                      map_noise_inv[ii,jj,:,:,:,:])
        stop = time_module.clock()
        #print "Constructing map noise took %5.2f seconds." % (stop - start)
        self.assertTrue(sp.allclose(map_noise_inv, tmp_map_noise_inv))

    def test_mean_mode_equivalent(self):
        """Test 2 equivalent ways to deweight mean, see if they agree."""
        self.DM.nf = 1
        time_stream, ra, dec, az, el, time, mask_inds = \
                                               self.DM.get_all_trimmed()
        nt = len(time)
        # Frist Way.
        Noise1 = dirty_map.Noise(time_stream, time)
        thermal_noise_levels = sp.zeros((1,)) + 0.04  # Kelvin**2
        Noise1.add_thermal(thermal_noise_levels)
        Noise1.add_mask(mask_inds)
        Noise1.deweight_time_mean()
        Noise1.deweight_time_slope()
        Noise1.add_correlated_over_f(0.01, -1.2, 0.1)
        Noise1.finalize()
        N1 = Noise1.get_inverse()
        # Second Way.
        Noise2 = dirty_map.Noise(time_stream, time)
        thermal_noise_levels = sp.zeros((1,)) + 0.04  # Kelvin**2
        Noise2.add_thermal(thermal_noise_levels)
        Noise2.add_mask(mask_inds)
        Noise2.deweight_time_slope()
        Noise2.add_correlated_over_f(0.01, -1.2, 0.1)
        Noise2.freq_mode_noise += dirty_map.T_large
        Noise2.finalize()
        N2 = Noise2.get_inverse()
        N2_m = N2.view()
        N2_m.shape = (nt, nt)
        self.assertTrue(sp.allclose(N2, N1))

    def deactivated_test_numerical_stability(self):
        nf = 20
        nt = 150
        dt = 0.3
        BW = 1. / dt / 2.
        time_stream = sp.zeros((nf, nt))
        time_stream = al.make_vect(time_stream, axis_names=("freq", "time"))
        time = dt * (sp.arange(nt) + 50)
        N = dirty_map.Noise(time_stream, time)
        # Thermal.
        thermal = sp.zeros(nf, dtype=float) + 0.0003 * BW * 2.
        N.add_thermal(thermal)  # K**2
        # Mask.
        mask_f = []
        mask_t = []
        # Mask out a whole time.
        mask_f += range(nf)
        mask_t += [74] * nf
        # Mask out a channel.
        mask_f += [13] * nt
        mask_t += range(nt)
        # Mask out a few more points.
        mask_f += [7, 15, 9, 7]
        mask_t += [87, 123, 54, 130]
        N.add_mask((sp.array(mask_f, dtype=int), sp.array(mask_t, dtype=int)))
        # Time mean and slope.
        #N.deweight_time_mean()
        #N.deweight_time_slope()
        # Correlated modes.
        f_0 = 1.
        amps = [0.01, 0.001, 0.0001]
        index = [-2.5, -1.7, -1.2]
        for ii in range(3):
            mode = sp.arange(nf, dtype=float)
            mode += 5 * nf * ii
            mode *= ii * 2. * sp.pi / nf
            mode = sp.cos(mode)
            mode /= sp.sqrt(sp.sum(mode**2))
            N.add_over_f_freq_mode(amps[ii], index[ii], f_0, 
                                   0.0003 * BW * 2., mode)
        # All freq modes over_f.
        N.add_all_chan_low(thermal, -0.9, 0.01)
        N.finalize()
        N_mat = N.get_inverse()
        N_mat.shape = (nf * nt,) * 2
        e, v = linalg.eigh(N_mat)
        print sp.amax(e) / sp.amin(e)


    def test_uncoupled_channels(self):
        time_stream, ra, dec, az, el, time, mask_inds = \
                                               self.DM.get_all_trimmed()
        nt = len(time)
        map = self.map
        # Calculate two noise matrices.  Niether have noise couplig terms, but
        # the second is calculated in the most general way.
        # The first one.
        P = dirty_map.Pointing(("ra", "dec"), (ra, dec), map, 'nearest')
        Noise1 = dirty_map.Noise(time_stream, time)
        thermal_noise_level = 0.4  # Kelvin**2
        Noise1.add_thermal(thermal_noise_level)
        Noise1.add_mask(mask_inds)
        Noise1.deweight_time_mean()
        Noise1.deweight_time_slope()
        Noise1.finalize(frequency_correlations=False)
        # The second one.
        Noise2 = dirty_map.Noise(time_stream, time)
        Noise2.add_thermal(thermal_noise_level)
        Noise2.add_mask(mask_inds)
        Noise2.deweight_time_mean()
        Noise2.deweight_time_slope()
        Noise2.finalize(frequency_correlations=True)
        # Check that they agree.
        N_inv = Noise1.get_inverse()
        self.assertEqual(N_inv.shape, (nf_d, nt, nt))
        N_inv2 = Noise2.get_inverse()
        for ii in range(nf_d):
            self.assertTrue(sp.allclose(N_inv[ii,...], N_inv2[ii,:,ii,:]))
        # Check that Noise1 is right.
        # Calculate the noise on a slice and make sure it multiplies to give
        # the identity.
        N = sp.zeros((nf_d, nt * nt))
        N[:,::nt + 1] = Noise1.diagonal
        N.shape = (nf_d, nt, nt)
        for ii in range(Noise1.time_mode_noise.shape[0]):
            N += (Noise1.time_mode_noise[ii,...].flat[::nf_d + 1][:,None,None]
                  * Noise1.time_modes[ii,None,:,None]
                  * Noise1.time_modes[ii,None,None,:])
        eye = sp.eye(nt)
        for ii in range(nf_d):
            tmp_eye = sp.dot(N[ii,...], N_inv[ii,...])
            self.assertTrue(sp.allclose(tmp_eye, eye))
        # Now test that weighting the time stream is the same for both
        # algorithms.
        t1 = Noise1.weight_time_stream(time_stream)
        t2 = Noise2.weight_time_stream(time_stream)
        self.assertTrue(sp.allclose(t1, t2))
        # Now transform to the map domain and make sure they still match.
        # Noise1.
        map_noise_inv1 = sp.zeros((nf_d, nra_d, ndec_d, nra_d, ndec_d),
                                 dtype=float)
        map_noise_inv1 = al.make_mat(map_noise_inv1, axis_names=('freq', 'ra', 
            'dec', 'ra', 'dec'), row_axes=(0, 1, 2), 
            col_axes=(0, 3, 4))
        for ii in xrange(nf_d):
            P.noise_channel_to_map(Noise1, ii, map_noise_inv1[ii,...])
        # Noise2.
        map_noise_inv2 = sp.zeros((nf_d, nra_d, ndec_d, nf_d, nra_d, ndec_d),
                                 dtype=float)
        map_noise_inv2 = al.make_mat(map_noise_inv2, axis_names=('freq', 'ra', 
            'dec', 'freq', 'ra', 'dec'), row_axes=(0, 1, 2), 
            col_axes=(3, 4, 5))
        for ii in xrange(nf_d):
            for jj in xrange(nra_d):
                P.noise_to_map_domain(Noise2, ii, jj,
                                      map_noise_inv2[ii,jj,:,:,:,:])
        # Check them.
        for ii in xrange(nf_d):
            self.assertTrue(sp.allclose(map_noise_inv1[ii,...],
                                        map_noise_inv2[ii,:,:,ii,:,:]))

    def deactivate_test_profile(self):
        """Not an actual test, this is for profiling."""
        
        nf = 32
        nra = 32
        ndec = 32
        DM = DataMaker(nscans=10, nt_scan=200, nf=nf, nra=nra,
                 ndec=ndec, scan_size=5.0, map_size=7.0,
                 add_noise=False, add_ground=False)
        map = DM.get_map()
        time_stream, ra, dec, az, el, time, mask_inds = DM.get_all_trimmed()
        P = dirty_map.Pointing(("ra", "dec"), (ra, dec), map, 'linear')
        Noise = dirty_map.Noise(time_stream, time)
        thermal_noise_levels = sp.zeros((nf)) + 0.04  # Kelvin**2
        Noise.add_thermal(thermal_noise_levels)
        Noise.add_mask(mask_inds)
        self.assertTrue(sp.alltrue(Noise.diagonal[mask_inds] > 10))
        Noise.deweight_time_mean()
        Noise.deweight_time_slope()
        Noise.add_correlated_over_f(0.01, -1.2, 0.1)
        start = time_module.clock()
        Noise.finalize()
        stop = time_module.clock()
        print "Finalizing noise took %5.2f seconds." % (stop - start)
        # Do the profiling.
        map_noise_inv = sp.zeros((nf, nra, ndec, nf, nra, ndec),
                                 dtype=float)
        print "Frequency ind:",
        start = time_module.clock()
        for ii in xrange(1):
            print ii,
            for jj in xrange(nra):
                P.noise_to_map_domain(Noise, ii, jj, 
                                      map_noise_inv[ii,jj,:,:,:,:])
        stop = time_module.clock()
        print
        print "Constructing map noise took %5.2f seconds." % (stop - start)


class TestEngine(unittest.TestCase):
#class TestEngine(object):

    def test_over_f(self):
        """Test with over f dominated noise."""

        nf = 2
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        # These parameters.ers are teetering on the edge of numerical stability.
        thermal_var = 0.00005 * BW_d * 2
        over_f_pars = (0.0001, -0.95, 5.0)
        
        Data = DataMaker(nscans=6, nt_scan=50, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=0., random_slope=0.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 20)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > -3)
        self.assertTrue(nor < 3)

    def test_over_f_fails(self):
        """Test with over f dominated noise and make sure that chis squared is
        wrong if the input noise is wrong."""

        nf = 2
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        # These parameters are teetering on the edge of numerical stability.
        thermal_var = 0.00004 * BW_d * 2
        over_f_pars = (0.0001, -0.95, 5.0)
        over_f_noise_model = (0.00005, -0.95, 5.0)
        
        Data = DataMaker(nscans=6, nt_scan=50, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=0., random_slope=0.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 20,
                                                over_f=over_f_noise_model)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > 0)

    def test_normal(self): 
        """Test under 'normal' conditions."""
        
        nf = 20
        nra = 5
        ndec = 5
        map_size = 0.6
        scan_size = 0.7
        thermal_var = 0.00005 * BW_d * 2
        over_f_pars = [(.0003, -1.4, 1., 0.0001), (.0001, -1.1, 1., 0.0002),
                       (.00003, -0.8, 1., 0.0003)]
        universal_pars = (-1.5, 0.1)
        
        Data = DataMaker(nscans=6, nt_scan=20, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=None,
                 random_mean=10., random_slope=10.,
                 add_ground=False, freq_mode_noise=over_f_pars,
                 universal_over_f=universal_pars)
        new_map, map_diffs, nor = self.make_map(Data, 80)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > -3)
        self.assertTrue(nor < 3)

    def test_all_freq(self): 
        """Test under all freq dominated conditions."""
        
        nf = 20
        nra = 5
        ndec = 5
        map_size = 0.6
        scan_size = 0.5
        thermal_var = 0.00005 * BW_d * 2 * (5. + sp.arange(nf)/nf) / 5.
        over_f_pars = [(.00003, -1.4, 1., 0.00001)]
        # Spectra needs to be pretty steep since we cut it off 1 octave after
        # the corner.
        universal_pars = (-1.8, 0.2)
        
        Data = DataMaker(nscans=6, nt_scan=20, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=None,
                 random_mean=10., random_slope=10.,
                 add_ground=False, freq_mode_noise=over_f_pars,
                 universal_over_f=universal_pars)
        new_map, map_diffs, nor = self.make_map(Data, 80)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > -3)
        self.assertTrue(nor < 3)

    def test_all_freq_fails(self): 
        """Test under that the chi squared is wrong for bad all_freq model.
        """
        
        nf = 20
        nra = 5
        ndec = 5
        map_size = 0.6
        scan_size = 0.7
        thermal_var = 0.00005 * BW_d * 2 * (5. + sp.arange(nf)/nf) / 5.
        over_f_pars = [(.00003, -1.4, 1., 0.00001)]
        universal_pars = (-1.6, 0.25)
        universal_pars_model = (-1.9, 0.15)
        
        Data = DataMaker(nscans=6, nt_scan=20, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=None,
                 random_mean=10., random_slope=10.,
                 add_ground=False, freq_mode_noise=over_f_pars,
                 universal_over_f=universal_pars)
        new_map, map_diffs, nor = self.make_map(Data, 80, 
                                         all_channel=universal_pars_model)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > 1)

    def test_freq_modes_fails(self):
        nf = 10
        nra = 5
        ndec = 5
        map_size = 0.6
        scan_size = 0.7
        thermal_var = 0.00005 * BW_d * 2
        over_f_pars = [(.003, -1.4, 1., 0.0001), (.001, -1.1, 1., 0.0002),
                       (.0003, -0.8, 1., 0.0003)]
        over_f_noise = [(.009, -1.7, 1., 0.0001), (.006, -1.1, 1., 0.0002),
                       (.002, -0.5, 1., 0.0003)]
        universal_pars = (-0.7, 0.)
        
        Data = DataMaker(nscans=6, nt_scan=20, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=None,
                 random_mean=10., random_slope=10.,
                 add_ground=False, freq_mode_noise=over_f_noise,
                 universal_over_f=universal_pars)
        new_map, map_diffs, nor = self.make_map(Data, 80, 
                                                over_f_modes=over_f_pars)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > 1)
        
    def test_freq_modes_fails_thermal(self):
        nf = 10
        nra = 5
        ndec = 5
        map_size = 0.6
        scan_size = 0.7
        thermal_var = 0.00005 * BW_d * 2
        over_f_pars = [(.0003, -1.4, 1., 0.001), (.0001, -1.1, 1., 0.002),
                       (.00003, -0.8, 1., 0.003)]
        over_f_noise = [(.0003, -1.4, 1., 0.005), (.0001, -1.1, 1., 0.008),
                       (.00003, -0.5, 1., 0.015)]
        universal_pars = (-0.7, 0.)
        
        Data = DataMaker(nscans=6, nt_scan=20, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=None,
                 random_mean=10., random_slope=10.,
                 add_ground=False, freq_mode_noise=over_f_noise,
                 universal_over_f=universal_pars)
        new_map, map_diffs, nor = self.make_map(Data, 80, 
                                                over_f_modes=over_f_pars)
        expected_map = new_map - map_diffs
        # Very weak test just to make sure we are in the ball park.
        self.assertTrue(nor > -1)

    def test_mean(self):
        """Test the mean subtration."""
        
        nf = 3
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        thermal_var = 0.0005 * BW_d * 2
        over_f_pars = (.0000001, -0.95, 0.1)
        
        Data = DataMaker(nscans=2, nt_scan=10, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=100., random_slope=0.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 40)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > -3)
        self.assertTrue(nor < 3)

    def test_mean_fails(self):
        """Test that mean subtraction breaks down."""
        
        nf = 3
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        thermal_var = 0.0005 * BW_d * 2
        over_f_pars = (.000000001, -0.95, 0.1)
        
        Data = DataMaker(nscans=2, nt_scan=10, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=100000., random_slope=0.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 40)
        expected_map = new_map - map_diffs
        # Pretty weak test since this only affects a small number of modes.
        self.assertTrue(nor > 2)

    def test_slope(self):
        """Test slope subtraction."""
        
        nf = 3
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        thermal_var = 0.0005 * BW_d * 2
        over_f_pars = (.0000001, -0.95, 0.1)
        
        Data = DataMaker(nscans=2, nt_scan=10, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=0., random_slope=100.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 40)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > -3)
        self.assertTrue(nor < 3)

    def deactivated_test_profile(self):
        """Profiles the map maker under realistic conditions."""
        
        nf = 32
        nra = 32
        ndec = 32
        map_size = 5.
        scan_size = 5.
        thermal_var = 0.01
        over_f_pars = (.01, -0.95, 0.1)
        
        Data = DataMaker(nscans=8, nt_scan=150, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=10., random_slope=10.,
                 add_ground=False)
        t, t2 = self.make_map(Data, 1, True)
        print "Took %7.2f core seconds and %7.2f real seconds." % (t, t2)

    def test_slope_fails(self):
        """Test that slope subtraction fails."""

        nf = 3
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        thermal_var = 0.0005 * BW_d * 2
        over_f_pars = (.000000001, -0.95, 0.1)
        
        Data = DataMaker(nscans=2, nt_scan=10, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=0., random_slope=100000.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 40)
        expected_map = new_map - map_diffs
        # Pretty weak test since this only affects a small number of modes.
        self.assertTrue(nor > 1)

    def test_untouched_pixels(self):

        nf = 2
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 2.0
        thermal_var = 0.0005 * BW_d * 2
        over_f_pars = (.001, -0.95, 0.1)
        
        # Small number of densly sampled scans.  The pixels that are hit are
        # hit multiple time such that the interpolation degeneracy doesn't kick
        # in.  Our untouched pixel fixer only fixes pixels that are completly
        # untouched, not degenerate. 
        Data = DataMaker(nscans=3, nt_scan=400, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=10., random_slope=10.,
                 add_ground=False)
        # This will fail if the noise matrix is singular.
        new_map, map_diffs, nor = self.make_map(Data, 6)
        expected_map = new_map - map_diffs

    def make_map(self, Data, n_data=10, time=False, thermal=None, over_f=None,
                 over_f_modes=None, all_channel=None):
        """This does all the map making for the input parameters.  Contains no
        tests."""
        
        nf = Data.nf
        nra = Data.nra
        ndec = Data.ndec
        map_size = Data.map_size
        scan_size = Data.scan_size
        # Figure out the thermal values.
        if thermal is None:
            thermal_var = Data.thermal
        else:
            thermal_var = thermal
        # Figure out what the noise frequeny modes are.
        if over_f is None:
            over_f_pars = Data.correlated_noise
        else:
            over_f_pars = over_f
        if over_f_modes is None:
            over_f_mode_pars = Data.freq_mode_noise
        else:
            over_f_mode_pars = over_f_modes
        # Figure out what to use all the all channel noise.
        if all_channel is None:
            all_channel_pars = Data.universal_over_f
        else:
            all_channel_pars = all_channel
        # Figure out what noise model will be used.
        if (not over_f is None) or (over_f_modes is None and not
                                    Data.correlated_noise is None):
            n_modes = 0
            noise_model = 'mean'
        else :
            n_modes = len(over_f_mode_pars)
            noise_model = 'measured'

        params = {
            'dm_file_middles' : ['a'] * n_data,
            'dm_map_shape' : (nra, ndec),
            'dm_field_centre' : (21., 0.),
            'dm_pixel_spacing' : map_size/nra,
            'dm_time_block' : 'file',
            'dm_frequency_correlations' : noise_model,
            'dm_number_frequency_modes' : n_modes,
            'dm_deweight_time_mean' : True,
            'dm_deweight_time_slope' : True
                  }
        Maker = dirty_map.DirtyMapMaker(params, feedback=0)
        Maker.pol = 0
        Maker.band = 0
        # Replace the data reader with mock data generator.
        def iterate_data(file_middles):
            for ii in xrange(n_data):
                yield Data.get_all_trimmed(ii)
        Maker.iterate_data = iterate_data
        # Replace the noise parameter reader.
        def get_noise_parameter(MakerClass, parameter_name):
            if parameter_name == "thermal":
                thermal = sp.empty(nf)
                thermal[...] = thermal_var
                return thermal
            elif parameter_name == "mean_over_f":
                return over_f_pars
            elif parameter_name[:12] == "over_f_mode_":
                mode_number = int(parameter_name[12:])
                mode = make_frequency_mode(mode_number, nf)
                return over_f_mode_pars[mode_number] + (mode,)
            elif parameter_name == "all_channel":
                thermal = sp.empty(nf)
                thermal[...] = thermal_var / BW_d / 2.
                return (thermal_var,) + all_channel_pars
            else :
                raise ValueError()
        Maker.get_noise_parameter = lambda p_name : get_noise_parameter(Maker,
                                                                        p_name)
        Maker.uncorrelated_channels = False
        Maker.n_chan = nf
        Maker.n_ra = nra
        Maker.n_dec = ndec
        Maker.n_processes = 8
        # Initialize the map maker's map and inverse covariance.
        map = Data.get_map()
        Maker.map = map
        # I use ones instead of zeros to make sure it gets overwritten.
        cov_inv = sp.ones((nf, nra, ndec, nf, nra, ndec),
                           dtype=float)
        cov_inv = al.make_mat(cov_inv, axis_names=('freq', 'ra', 'dec',
            'freq', 'ra', 'dec'), row_axes=(0, 1, 2), col_axes=(3, 4, 5))
        cov_inv[...] = 0
        Maker.cov_inv = cov_inv
        # Run the engine of the map maker.
        start = time_module.clock()
        real_start = time_module.time()
        Maker.make_map()
        t = time_module.clock() - start
        t2 = time_module.time() - real_start
        if time:
            return t, t2
        # Clean the map.
        cov_inv_m = cov_inv.view()
        cov_inv_m.shape = (nf * nra * ndec, nf * nra * ndec)
        map_v = map.view()
        map_v.shape = (nf * nra * ndec,)
        v, U = linalg.eigh(cov_inv_m)
        self.assertTrue(sp.all(v > 0))
        new_map = sp.dot(U.T, map_v)
        new_map /= v
        new_map = sp.dot(U, new_map)
        new_map.shape = (nf, nra, ndec)
        # Get the expected map.
        map_ra = map.get_axis('ra')[:,None]
        map_dec = map.get_axis('dec')[None,:]
        expected_map = Data.sky_ra_dec(map_ra, map_dec)
        # Figure out the differences (that aren't in the slice mean mode).
        map_diffs = new_map - expected_map
        diff_v = map_diffs.view()
        diff_v.shape = (nf * nra * ndec,)
        chi_sq = sp.dot(sp.dot(diff_v, cov_inv_m), diff_v)
        n_deg_f = nf * nra * ndec
        # Transformation to make it normally distributed.
        nor = (chi_sq/n_deg_f)**(1.0/3) - (1.0 - 2.0/9/n_deg_f)
        nor /= sp.sqrt(2.0/9.0/n_deg_f)
        return new_map, map_diffs, nor


class TestPreprocessor(unittest.TestCase):

    def setUp(self):
        Reader = fitsGBT.Reader("./testdata/testfile_guppi_combined.fits",
                                feedback=0)
        self.Blocks = Reader.read((), 1)
        Data = self.Blocks[0]
        Data.calc_freq()
        Maker = dirty_map.DirtyMapMaker({}, feedback=0)
        n_chan = Data.dims[-1]
        Maker.n_chan = n_chan
        Maker.pols = (1, 2, 3, 4)
        Maker.pol_ind = 0
        Maker.band_centres = (Data.freq[Data.dims[-1]//2],)
        Maker.band_ind = 0
        map = sp.zeros((Data.dims[-1], 32, 15))
        map = al.make_vect(map, ('freq', 'ra', 'dec'))
        map.set_axis_info('freq', Data.freq[Data.dims[-1]//2],
                          Data.field['CRVAL1'])
        map.set_axis_info('ra', 218, 0.075)
        map.set_axis_info('dec', 2, 0.075)
        Maker.map = map
        self.Maker = Maker
        # The variances of each channel.
        self.norms = (sp.arange(1., 2., 0.25)[:,None]
                      * (sp.arange(1., 2., 1./n_chan)[None,:]))
        for Data in self.Blocks:
            Data.data[...] = random.randn(*Data.data.shape)
            Data.data *= sp.sqrt(self.norms[:,None,:])
            Data.data += 50.

    def test_basic(self):
        time_stream, ra, dec, az, el, time, mask_inds = \
                self.Maker.preprocess_data(self.Blocks)
        self.assertTrue(sp.allclose(sp.mean(time_stream, 1), 0, atol=0.2))
        self.assertTrue(sp.allclose(sp.var(time_stream, 1), self.norms[0,:], 
                                    rtol=0.4))
        self.assertTrue(sp.allclose(self.Maker.channel_vars,
                                    self.norms[0,:], rtol=0.4))

    def test_masked_channel(self):
        for Data in self.Blocks:
            Data.data[:,:,:,13] = ma.masked
            Data.data[:,:,:,17] = 0.
        time_stream, ra, dec, az, el, time, mask_inds = \
                self.Maker.preprocess_data(self.Blocks)
        self.assertTrue(sp.allclose(self.Maker.channel_vars[[13,17]],
                                    dirty_map.T_infinity))

    def test_masked_inds(self):
        # Since the preprocessor moves stuff around, put a high value beside
        # all the masked ends so we can find it.
        self.Blocks[0].data[11,0,0,12] = ma.masked
        self.Blocks[0].data[11,0,0,13] = 100
        self.Blocks[0].data[6,0,0,2] = ma.masked
        self.Blocks[0].data[6,0,0,3] = 100
        self.Blocks[0].data[37,0,0,21] = ma.masked
        self.Blocks[0].data[37,0,0,22] = 100
        self.Blocks[0].data[73,0,0,6] = ma.masked
        self.Blocks[0].data[73,0,0,7] = 100
        self.Blocks[0].data[134,0,0,9] = ma.masked
        self.Blocks[0].data[134,0,0,10] = 100
        time_stream, ra, dec, az, el, time, mask_inds = \
                self.Maker.preprocess_data(self.Blocks)
        self.assertTrue(len(mask_inds[0]) > 0)
        for ii in range(len(mask_inds[0])):
            t_ind = mask_inds[1][ii]
            f_ind = mask_inds[0][ii]
            self.assertTrue(time_stream[f_ind + 1, t_ind] > 40)


class TestModuleIO(unittest.TestCase):
#class TestModuleIO(object):
    
    def setUp(self):

        nra = 16
        ndec = 8
        self.params = {
            'dm_input_root' : './testdata/',
            'dm_file_middles' : ('testfile_guppi_rotated',),
            'dm_input_end' : '.fits',
            'dm_output_root' : './testoutput_',
            'dm_noise_parameter_file' :
                    './testdata/testfile_guppi_noise_parameters.shelve',
            'dm_scans' : (),
            'dm_IFs' : (),
            'dm_polarizations' : ('I',),
            'dm_map_shape' : (nra, ndec),
            'dm_field_centre' : (218., 2.),
            'dm_pixel_spacing' : 0.3,
            'dm_time_block' : 'file',
            'dm_frequency_correlations' : 'measured',
            'dm_number_frequency_modes' : 4,
            'dm_deweight_time_mean' : True,
            'dm_deweight_time_slope' : True
                  }
    
    def test_normal(self):
        params = self.params
        dirty_map.DirtyMapMaker(params, feedback=0).execute(4)
        files = glob.glob('*testout*')
        # 17 files = 2 pols * 2 bands * 2 (map and noise) * 2 (.npy and
        # .npy.meta) + 1 (parameter file)
        self.assertTrue(len(files) == 9)

    def test_separate_scans(self) :
        params = self.params
        params['dm_time_block'] = 'scan'
        dirty_map.DirtyMapMaker(params, feedback=0).execute(4)
        files = glob.glob('*testout*')
        # 17 files = 2 pols * 2 bands * 2 (map and noise) * 2 (.npy and
        # .npy.meta) + 1 (parameter file)
        self.assertTrue(len(files) == 9)
    
    def test_independant_channels(self) :
        params = self.params
        params["dm_frequency_correlations"] = "None"
        dirty_map.DirtyMapMaker(params, feedback=0).execute(4)
        files = glob.glob('*testout*')
        # 17 files = 2 pols * 2 bands * 2 (map and noise) * 2 (.npy and
        # .npy.meta) + 1 (parameter file)
        self.assertTrue(len(files) == 9)

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)

if __name__ == '__main__' :
    unittest.main()
