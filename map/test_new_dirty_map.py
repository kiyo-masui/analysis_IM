"""Unit tests for dirty map maker."""

import os
import glob
import time as time_module

import unittest
import scipy as sp
import numpy.ma as ma
import matplotlib.pyplot as plt

import core.algebra as al
import new_dirty_map as dirty_map
import tools

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
                 add_noise=False, add_ground=False):
        # Store parameters.
        self.nscans = nscans
        self.nt_scan = nt_scan
        self.nf = nf
        self.nra = nra
        self.ndec = ndec
        self.scan_size = float(scan_size)
        self.map_size = float(map_size)
        self.add_noise = add_noise
        self.add_ground = add_ground
        # Initialize a counter.
        self.data_set_number = 0
        self.nt = nt_scan * nscans

    def ground_az_el(self, az, el):
        time_stream = sp.zeros((self.nf, self.nt))
        time_stream += sp.sin(az*0.5)
        time_stream += sp.cos(el*0.3)
        time_stream *= -sp.arange(self.nf, dtype=float)[:,None]/self.nf + 2.0
        return time_stream

    def sky_ra_dec(self, ra, dec):
        time_stream = sp.zeros((self.nf, self.nt))
        time_stream += sp.sin(ra*2.3)*sp.cos(dec*2.7) # The sky map.
        time_stream *= sp.cos(sp.arange(self.nf, dtype=float))[:,None] + 1.2
        return time_stream

    def get_time(self, number):
        time = (sp.arange(self.nt_scan) + sp.arange(self.nscans)[:,None] 
                * 1.5 * self.nt_scan)
        time = time.flat[...] + 2.0 * self.nt*number
        time = time * 0.2 + 123456.0
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

    def get_all(self, number=None):
        if number == None:
            number = self.data_set_number
            self.data_set_number += 1
        time = self.get_time(number)
        az, el = self.get_az_el(number)
        ra, dec = self.get_ra_dec(number)
        time_stream = sp.zeros((self.nf, self.nt))
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
        for ii in xrange(nf_d):
            for jj in xrange(Noise.time_modes.shape[0]):
                tmp_mat[ii,:,ii,:] += (Noise.time_modes[jj,:]
                                       * Noise.time_modes[jj,:,None]
                                       * dirty_map.T_infinity)
        # Here I assume that the only frequency noise mode is the mean mode.
        # TODO: Make more like the above loop for time modes.
        tmp_mat += Noise.freq_mode_noise[0,None,:,None,:]/nf_d
        tmp_mat.shape = (nt*nf_d, nt*nf_d)
        # Check that the matrix I built for testing is indeed symetric.
        self.assertTrue(sp.allclose(tmp_mat, tmp_mat.transpose()))
        noise_inv = Noise.get_inverse()
        noise_inv.shape = (nt*nf_d, nt*nf_d)
        # Check that the production matrix is symetric.
        self.assertTrue(sp.allclose(noise_inv, noise_inv.transpose()))
        tmp_eye = sp.dot(tmp_mat, noise_inv)
        noise_inv.shape = (nf_d, nt, nf_d, nt)
        self.assertTrue(sp.allclose(tmp_eye, sp.identity(nt*nf_d)))
        #### Test the noise weighting of the data.
        noise_weighted_data = Noise.noise_weight_time_stream(time_stream)
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
        Noise.set_pointing(P)
        start = time_module.clock()
        for ii in xrange(nf_d):
            for jj in xrange(nra_d):
                for kk in xrange(ndec_d):
                    Noise.update_map_noise(ii, (jj, kk),
                                           map_noise_inv[ii,jj,kk,:,:,:])
        stop = time_module.clock()
        #print "Constructing map noise took %5.2f seconds." % (stop - start)
        self.assertTrue(sp.allclose(map_noise_inv, tmp_map_noise_inv))
        # Check the other way of doing it.
        map_noise_inv[...] = 0
        for ii in xrange(nf_d):
            Noise.update_map_noise(ii, None, map_noise_inv[ii,:,:,:,:,:])
        self.assertTrue(sp.allclose(map_noise_inv, tmp_map_noise_inv))

    def test_profile(self):
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
        Noise.set_pointing(P)
        stop = time_module.clock()
        print "Finalizing noise took %5.2f seconds." % (stop - start)
        # Do the profiling.
        map_noise_inv = sp.zeros((nf, nra, ndec, nf, nra, ndec),
                                 dtype=float)
        print "Frequency ind:",
        start = time_module.clock()
        for ii in xrange(1):
            print ii,
            Noise.update_map_noise(ii, None, map_noise_inv[ii,:,:,:,:,:])
        stop = time_module.clock()
        print
        print "Constructing map noise took %5.2f seconds." % (stop - start)

class TestMain(unittest.TestCase):

    def test_full(self):

        nf = 5
        nra = 20
        ndec = 20
        map_size = 5.
        scan_size = 6.

        Data = DataMaker(nscans=6, nt_scan=200, nf=nf, nra=nra,
                 ndec=ndec, scan_size=scan_size, map_size=map_size,
                 add_noise=False, add_ground=False)
        params = {
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
        def preprocess_data(MakerClass):
            return Data.get_all_trimmed(MakerClass.file_number)
        Maker.preprocess_data = lambda : preprocess_data(Maker)
        # Replace the noise parameter reader.
        def get_noise_parameter(MakerClass, parameter_name):
            if parameter_name == "thermal":
                thermal = sp.empty(nf)
                thermal[...] = 0.04
                return thermal
            elif parameter_name == "mean_over_f":
                return (0.01, -1.2, 0.1)
            else :
                raise ValueError()
        Maker.get_noise_parameter = lambda p_name : get_noise_parameter(Maker,
                                                                        p_name)
        Maker.n_chan = nf
        # Initialize the map maker's map and inverse covariance.
        Maker.map = Data.get_map()
        cov_inv = sp.zeros((nf_d, nra_d, ndec_d, nf_d, nra_d, ndec_d),
                           dtype=float)
        cov_inv = al.make_mat(cov_inv, axis_names=('freq', 'ra', 'dec',
            'freq', 'ra', 'dec'), row_axes=(0, 1, 2), col_axes=(3, 4, 5))
        Maker.cov_inv = cov_inv
        # Run the engine of the map maker.
        Maker.make_map()




if __name__ == '__main__' :
    unittest.main()
