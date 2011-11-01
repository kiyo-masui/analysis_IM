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
import new_dirty_map as dirty_map
from noise import noise_power

# Default data parameters used for most of the tests.
nscans_d = 6
nt_scan_d = 50
nt_d = nt_scan_d * nscans_d
nf_d = 5
scan_size_d = 5.0

nra_d = 20
ndec_d = 20
map_size_d = 6.0


class DataMaker(object):
    """Class that generates test data."""

    def __init__(self, nscans=nscans_d, nt_scan=nt_scan_d, nf=nf_d, nra=nra_d,
                 ndec=ndec_d, scan_size=scan_size_d, map_size=map_size_d,
                 thermal=0, correlated_noise=None, random_mean=0,
                 random_slope=0, add_ground=False):
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
        self.random_mean = random_mean
        self.random_slope = random_slope
        # Initialize a counter.
        self.data_set_number = 0
        self.nt = nt_scan * nscans
        self.dt = 0.2

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
        ra += 21 + (((sp.pi * (number+3)) % 1.0) - 0.5) * 1.0
        dec += 0 + (((sp.pi * (number+7)) % 1.0) - 0.5) * 1.0
        return ra, dec

    def get_map(self):
        map = sp.zeros((self.nf, self.nra, self.ndec))
        map = al.make_vect(map, ('freq', 'ra', 'dec'))
        map.set_axis_info('freq', 800e6, 1e6)
        map.set_axis_info('ra', 21, self.map_size/self.nra)
        map.set_axis_info('dec', 0, self.map_size/self.ndec)
        return map

    def get_noise(self, number):
        noise_data = random.randn(self.nf, self.nt) * sp.sqrt(self.thermal)
        t = self.get_time(number)
        # Noisy mean mode.
        noise_data += random.randn(self.nf, 1) * sp.sqrt(self.random_mean)
        # Noisy slope mode.
        slope_noise = random.randn(self.nf, 1) * sp.sqrt(self.random_slope)
        slope_mode = (sp.sqrt(self.nt) * (t - sp.mean(t)) 
                       / sp.sqrt(sp.sum((t - sp.mean(t))**2)))
        noise_data += slope_noise * slope_mode
        if not self.correlated_noise is None:
            amp = self.correlated_noise[0]
            index = self.correlated_noise[1]
            f0 = self.correlated_noise[2]
            over_f = noise_power.generate_overf_noise(amp, index, f0,
                                                      self.dt, self.nt*2)
            ov_f_interpolator = interpolate.interp1d(t[0] 
                                    + sp.arange(self.nt*2)*self.dt, over_f)
            noise_data += ov_f_interpolator(t)
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
        Noise2.freq_mode_noise += dirty_map.T_infinity
        Noise2.finalize()
        N2 = Noise2.get_inverse()
        N2_m = N2.view()
        N2_m.shape = (nt, nt)
        self.assertTrue(sp.allclose(N2, N1))

    def test_uncoupled_channels(self):
        time_stream, ra, dec, az, el, time, mask_inds = \
                                               self.DM.get_all_trimmed()
        nt = len(time)
        Noise = dirty_map.Noise(time_stream, time)
        thermal_noise_level = 0.04  # Kelvin**2
        Noise.add_thermal(thermal_noise_level)
        Noise.add_mask(mask_inds)
        Noise.deweight_time_mean()
        Noise.deweight_time_slope()
        Noise.add_correlated_over_f(0.01, -1.2, 0.1)
        self.assertRaises(NotImplementedError, Noise.finalize,
                          frequency_correlations=False)
        return
        Noise.finalize(frequency_correlations=False)
        N_inv = Noise.get_inverse()
        self.assertEqual(N_inv.shape, (nf_d, nt, nt))
        # Calculate the noise on a slice and make sure it multiplies to give
        # the identity.
        N = sp.zeros(nf_d, nt * nt)
        N[:,::nt + 1] = thermal_noise_level
        for ii in range(Noise.time_mode_noise.shape[0]):
            N += (Noise.time_mode_noise[ii,...].flat[::nf_d + 1,None,None]
                  * Noise.time_modes[ii,None,:,None]
                  * Noise.time_modes[ii,None,None,:])
        for ii in xrange(nf_d):
            for jj in xrange(Noise.freq_modes.shape[0]):
                N[ii,...] += (Noise.freq_mode_noise[jj,...] 
                      * Noise.freq_modes[jj,ii]**2)
        eye = sp.eye(nt, nt)
        for ii in range(nf_d):
            tmp_eye = sp.dot(N, N_inv[ii,...])
            self.assertTrue(sp.allclose(tmp_eye, eye))

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

    def test_over_f(self):
        """Test with over f dominated noise."""

        nf = 1
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        # These parameters are teetering on the edge of numerical stability.
        thermal_var = 0.0001
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

    def test_normal(self):
        """Test under 'normal' conditions."""
        
        nf = 3
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        thermal_var = 0.001
        over_f_pars = (.001, -0.95, 0.1)
        
        Data = DataMaker(nscans=6, nt_scan=50, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=10., random_slope=10.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 20)
        expected_map = new_map - map_diffs
        # These statistical tests fail 1% of the time under normal conditions.
        self.assertTrue(nor > -3)
        self.assertTrue(nor < 3)

    def test_mean(self):
        """Test the mean subtration."""
        
        nf = 3
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 1.3
        thermal_var = 0.001
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
        thermal_var = 0.001
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
        thermal_var = 0.001
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
        thermal_var = 0.001
        over_f_pars = (.000000001, -0.95, 0.1)
        
        Data = DataMaker(nscans=2, nt_scan=10, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 thermal=thermal_var, correlated_noise=over_f_pars,
                 random_mean=0., random_slope=100000.,
                 add_ground=False)
        new_map, map_diffs, nor = self.make_map(Data, 40)
        expected_map = new_map - map_diffs
        # Pretty weak test since this only affects a small number of modes.
        self.assertTrue(nor > 3)

    def test_untouched_pixels(self):

        nf = 2
        nra = 10
        ndec = 10
        map_size = 1.2
        scan_size = 2.0
        thermal_var = 0.001
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

    def make_map(self, Data, n_data=10, time=False):
        """This does all the map making for the input parameters.  Contains no
        tests."""
        
        nf = Data.nf
        nra = Data.nra
        ndec = Data.ndec
        map_size = Data.map_size
        scan_size = Data.scan_size
        thermal_var = Data.thermal
        over_f_pars = Data.correlated_noise

        params = {
            'dm_file_middles' : ['a'] * n_data,
            'dm_map_shape' : (nra, ndec),
            'dm_field_centre' : (21., 0.),
            'dm_pixel_spacing' : map_size/nra,
            'dm_time_block' : 'file',
            'dm_frequency_correlations' : 'mean',
            'dm_thermal_weight' : 'thermal',
            'dm_deweight_time_mean' : True,
            'dm_deweight_time_slope' : True
                  }
        Maker = dirty_map.DirtyMap(params, feedback=0)
        Maker.pol = 0
        Maker.band = 0
        # Replace the data reader with mock data generator.
        def iterate_data():
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
                #return (0.01, -1.2, 0.1)
            else :
                raise ValueError()
        Maker.get_noise_parameter = lambda p_name : get_noise_parameter(Maker,
                                                                        p_name)
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


class TestModuleIO(unittest.TestCase):
    
    def setUp(self):

        nra = 16
        ndec = 8
        self.params = {
            'dm_input_root' : './testdata/',
            'dm_file_middles' : ('testfile_guppi_combined',),
            'dm_input_end' : '.fits',
            'dm_output_root' : './testoutput_',
            'dm_scans' : (),
            'dm_bands' : (),
            'dm_polarizations' : ('I', 'U'),
            'dm_map_shape' : (nra, ndec),
            'dm_field_centre' : (218., 2.),
            'dm_pixel_spacing' : 0.15,
            'dm_time_block' : 'file',
            'dm_frequency_correlations' : 'mean',
            'number_frequency_modes' : 1,
            'dm_noise_parameter_file' : "",
            'dm_thermal_weight' : 'thermal',
            'dm_deweight_time_mean' : True,
            'dm_deweight_time_slope' : True
                  }
    
    def deactivate_test_normal(self) :
        params = self.params
        dirty_map.DirtyMap(params, feedback=0).execute(4)
        files = glob.glob('*testout*')
        # 17 files = 2 pols * 2 bands * 2 (map and noise) * 2 (.npy and
        # .npy.meta) + 1 (parameter file)
        self.assertTrue(len(files) == 17)

    def deactivate_test_separate_scans(self) :
        params = self.params
        params['dm_time_block'] = 'scan'
        dirty_map.DirtyMap(params, feedback=0).execute(4)
        files = glob.glob('*testout*')
        # 17 files = 2 pols * 2 bands * 2 (map and noise) * 2 (.npy and
        # .npy.meta) + 1 (parameter file)
        self.assertTrue(len(files) == 17)

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)



if __name__ == '__main__' :
    unittest.main()
