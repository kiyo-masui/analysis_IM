"""Unit tests for dirty map maker."""

import os
import glob

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
nf_d = 10
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
        self.scan_size = scan_size
        self.map_size = map_size
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
        map.set_axis_info('ra', 21, 6.0/self.nra)
        map.set_axis_info('dec', 0, 6.0/self.ndec)
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


class TestClasses(unittest.TestCase):

    def setUp(self):
        self.DM = DataMaker()
        self.map = self.DM.get_map()
    
    def test_get_matrix(self):
        time_stream, ra, dec, az, el, time, mask_inds = self.DM.get_all()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), self.map, 'nearest')
        pointing_matrix = P.get_matrix()
        self.assertTrue(pointing_matrix.axes == ('time', 'ra', 'dec'))
        self.assertEqual(pointing_matrix.shape[0], ra.shape[0])
        self.assertEqual(self.map.shape[1], pointing_matrix.shape[1])
        self.assertEqual(self.map.shape[2], pointing_matrix.shape[2])

    def test_get_matrix_out_of_bounds(self):
        time_stream, ra, dec, az, el, time, mask_inds = self.DM.get_all()
        off_map_ind = 11
        ra[off_map_ind] = -22.0
        dec[off_map_ind] = 55.0
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), self.map, 'nearest')
        pointing_matrix = P.get_matrix()
        self.assertTrue(sp.allclose(pointing_matrix[off_map_ind,:], 0))
    
    def deactivated_test_apply_time_axis_nearest(self):
        # The function that this is testing hasn't been implemented. Test
        # disabled.
        time_stream, ra, dec, az, el, time, mask_inds = self.DM.get_all()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), self.map, 'nearest')
        pointing_matrix = P.get_matrix()
        gridded_data_mult = al.partial_dot(pointing_matrix.mat_transpose(),
                                           time_stream)
        gridded_data_fast = P.apply_to_time_axis(time_stream)
        self.assertTrue(sp.allclose(gridded_data_mult, gridded_data_fast))

    def test_build_noise(self):
        map = self.map
        time_stream, ra, dec, az, el, time, mask_inds = self.DM.get_all()
        Noise = dirty_map.Noise(time_stream, time)
        thermal_noise_levels = sp.zeros((nf_d)) + 0.04  # Kelvin**2
        Noise.add_thermal(thermal_noise_levels)
        Noise.add_mask(mask_inds)
        self.assertTrue(sp.alltrue(Noise.diagonal[mask_inds] > 10))
        Noise.deweight_time_mean()
        Noise.add_correlated_over_f(0.01, -1.2, 0.1)
        Noise.finalize()
        return
        #### Test that the first round of matrix inversion (using diagonal
        # noise bits) works.
        tmp_inv = Noise.get_diag_allfreq_inverse(2)
        tmp_mat = sp.copy(Noise.allfreq)
        tmp_mat.flat[::nt_d + 1] += Noise.diagonal[2,:]
        tmp_eye = sp.dot(tmp_mat, tmp_inv)
        self.assertTrue(sp.allclose(tmp_eye, sp.identity(nt_d)))
        #### Test the full inverse.
        tmp_mat = sp.zeros((nf_d, nt_d, nf_d, nt_d))
        tmp_mat.flat[::nt_d*nf_d + 1] += Noise.diagonal.flat
        for ii in xrange(nf_d):
            tmp_mat[ii,:,ii,:] += Noise.allfreq
        # Here I assume that the only frequency noise mode is the mean mode.
        tmp_mat += Noise.mode_noise[0,None,:,None,:]/nf_d
        tmp_mat.shape = (nt_d*nf_d, nt_d*nf_d)
        noise_inv = Noise.get_inverse()
        noise_inv.shape = (nt_d*nf_d, nt_d*nf_d)
        tmp_eye = sp.dot(tmp_mat, noise_inv)
        noise_inv.shape = (nf_d, nt_d, nf_d, nt_d)
        self.assertTrue(sp.allclose(tmp_eye, sp.identity(nt_d*nf_d)))
        #### Test the noise weighting of the data.
        noise_weighted_data = Noise.noise_weight_time_stream(time_stream)
        self.assertTrue(sp.allclose(noise_weighted_data, al.dot(noise_inv,
                                                                time_stream)))
        #### Test making noise in map space.
        # First make the noise matrix by brute force.
        P = dirty_map.Pointing(("ra", "dec"), (ra, dec), map, 'nearest')
        map_noise_inv = sp.zeros((nf_d, nra_d, ndec_d, nf_d, nra_d, ndec_d),
                                 dtype=float)
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
        for ii in xrange(nf_d):
            Noise.update_map_noise(P, ii, map_noise_inv[ii,:,:,:,:,:])
        self.assertTrue(sp.allclose(map_noise_inv, tmp_map_noise_inv))


if __name__ == '__main__' :
    unittest.main()
