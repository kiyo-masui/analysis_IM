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

nscans = 6
nt_scan = 50
nt = nt_scan * nscans
nf = 10

nra = 20
ndec = 20

# XXX: Prefilter that makes a time stream vect will deal with out of bounds
# pointings, adding infinite noise.  The pointing should not have to deal with
# this.


# First define a function of what the ground looks like in az-el.
def ground_az_el(az, el):
    time_stream = sp.zeros((nf, nt))
    time_stream += sp.sin(az*0.5)
    time_stream += sp.cos(el*0.3)
    time_stream *= -sp.arange(nf, dtype=float)[:,None]/nf + 2.0
    return time_stream

# Define a function that gives the sky in ra-dec.
def sky_ra_dec(ra, dec):
    time_stream = sp.zeros((nf, nt))
    time_stream += sp.sin(ra*2.3)*sp.cos(dec*2.7) # The sky map.
    time_stream *= sp.arange(nf, dtype=float)[:,None]/nf/2 + 1.0
    return time_stream


def make_test_data(ii=0):
    # Make some fake pointing data.  This isn't entirely consitant with
    # possible geometies on the sky, but that probably doesn't matter for
    # testing.
    
    # Time data.
    time = (sp.arange(nt_scan) + sp.arange(nscans)[:,None]*1.5*nt_scan 
            + 2.0*nt*ii).flat[...] * 0.2 + 123456.0
    # Elevation always constant with current scan strategy.
    el = sp.zeros(nt) + (13.13 * ii) % 60
    # Make az scan back and forth, 5 degrees.
    az = sp.zeros(nt)
    for ii in range(0, nscans, 2):
        az[ii*nt/nscans:(ii+1)*nt/nscans] = sp.arange(-nt/nscans/2.0,
                                                       nt/nscans/2.0,
                                                       dtype=float)
    for ii in range(1, nscans, 2):
        az[ii*nt/nscans:(ii+1)*nt/nscans] = -sp.arange(-nt/nscans/2.0, 
                                                        nt/nscans/2.0,
                                                        dtype=float)
    az /= sp.cos(el[0]*sp.pi/180.0)
    az += (26.64 * ii) % 290
    # Position angle of the telescope controls how much ra vs dec we scan.
    # Number between -70 and 70 degrees.
    # Using pi as an irrational number so we don't get repeats.
    pa = -70 + sp.pi*2*(ii+1) % 140
    # The angles in the drift direction.
    drifting = sp.arange(-nt/2.0, nt/2.0, dtype=float)/nt*5
    # The angles in the scan direction.
    scanning = sp.zeros(nt)
    for ii in range(0, nscans, 2):
        scanning[ii*nt//nscans:(ii+1)*nt//nscans] = sp.arange(-nt/nscans/2.0,
                                                              nt/nscans/2.0,
                                                              dtype=float)
    for ii in range(1, nscans, 2):
        scanning[ii*nt//nscans:(ii+1)*nt//nscans] = -sp.arange(-nt/nscans/2.0, 
                                                               nt/nscans/2.0,
                                                               dtype=float)
    scanning *= 5.0/(nt/nscans)
    # Use `pa` to rotate these to ra-dec.
    ra = sp.sin(pa*sp.pi/180)*drifting + sp.cos(pa*sp.pi/180)*scanning
    dec = sp.cos(pa*sp.pi/180)*drifting - sp.sin(pa*sp.pi/180)*scanning
    # Add maps centres and the some offset jitter.
    ra += 21 + (((sp.pi * (ii+3)) % 1.0) - 0.5) * 1.0
    dec += 0 + (((sp.pi * (ii+7)) % 1.0) - 0.5) * 1.0
    # Now make some data that has stuff that both correlates on the sky and on
    # the ground.
    time_stream = sp.zeros((nf, nt))
    time_stream += sky_ra_dec(ra, dec)
    time_stream += ground_az_el(az, el)
    time_stream = al.make_vect(time_stream, ('freq', 'time'))
    # Make of some mask_inds.
    mask_inds = (sp.empty(4, dtype=int), sp.empty(4, dtype=int))
    for jj in range(4):
        mask_inds[0][jj] = (13 * (ii*jj+1)) % nf
        mask_inds[1][jj] = (53 * (ii*jj+1)) % nf

    # Now make the map that we are going to accumulate this into.  If I did the
    # above correctly, the pointing should never be out of bounds.
    map = sp.zeros((nf, nra, ndec))
    map = al.make_vect(map, ('freq', 'ra', 'dec'))
    map.set_axis_info('freq', 800e6, 1e6)
    map.set_axis_info('ra', 21, 6.0/nra*sp.sqrt(2.0))
    map.set_axis_info('dec', 0, 6.0/nra*sp.sqrt(2.0))

    return map, time_stream, ra, dec, az, el, time, mask_inds


class TestPointing(unittest.TestCase):
    
    def test_get_matrix(self):
        map, time_stream, ra, dec, az, el, time, mask_inds = make_test_data()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), map, 'nearest')
        pointing_matrix = P.get_matrix()
        self.assertTrue(pointing_matrix.axes == ('time', 'ra', 'dec'))
        self.assertEqual(pointing_matrix.shape[0], ra.shape[0])
        self.assertEqual(map.shape[1], pointing_matrix.shape[1])
        self.assertEqual(map.shape[2], pointing_matrix.shape[2])
    
    def deactivated_test_apply_time_axis_nearest(self):
        # The function that this is testing hasn't been implemented. Test
        # disabled.
        map, time_stream, ra, dec, az, el, time, mask_inds = make_test_data()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), map, 'nearest')
        pointing_matrix = P.get_matrix()
        gridded_data_mult = al.partial_dot(pointing_matrix.mat_transpose(),
                                           time_stream)
        gridded_data_fast = P.apply_to_time_axis(time_stream)
        self.assertTrue(sp.allclose(gridded_data_mult, gridded_data_fast))


class TestNoiseClass(unittest.TestCase):

    def test_build_noise(self):
        map, time_stream, ra, dec, az, el, time, mask_inds = make_test_data()
        Noise = dirty_map.Noise(time_stream, time)
        thermal_noise_levels = sp.zeros((nf)) + 0.04  # Kelvin**2
        Noise.add_thermal(thermal_noise_levels)
        Noise.add_mask(mask_inds)
        self.assertTrue(sp.alltrue(Noise.diagonal[mask_inds] > 10))
        Noise.deweight_time_mean()
        Noise.add_correlated_over_f(0.01, -1.2, 0.1)
        Noise.finalize()
        #### Test that the first round of matrix inversion (using diagonal
        # noise bits) works.
        tmp_inv = Noise.get_diag_allfreq_inverse(2)
        tmp_mat = sp.copy(Noise.allfreq)
        tmp_mat.flat[::nt + 1] += Noise.diagonal[2,:]
        tmp_eye = sp.dot(tmp_mat, tmp_inv)
        self.assertTrue(sp.allclose(tmp_eye, sp.identity(nt)))
        #### Test the full inverse.
        tmp_mat = sp.zeros((nf, nt, nf, nt))
        tmp_mat.flat[::nt*nf+1] += Noise.diagonal.flat
        for ii in xrange(nf):
            tmp_mat[ii,:,ii,:] += Noise.allfreq
        # Here I assume that the only frequency noise mode is the mean mode.
        tmp_mat += Noise.mode_noise[0,None,:,None,:]/nf
        tmp_mat.shape = (nt*nf, nt*nf)
        noise_inv = Noise.get_inverse()
        noise_inv.shape = (nt*nf, nt*nf)
        tmp_eye = sp.dot(tmp_mat, noise_inv)
        noise_inv.shape = (nf, nt, nf, nt)
        self.assertTrue(sp.allclose(tmp_eye, sp.identity(nt*nf)))
        #### Test the noise weighting of the data.
        noise_weighted_data = Noise.noise_weight_time_stream(time_stream)
        self.assertTrue(sp.allclose(noise_weighted_data, al.dot(noise_inv,
                                                                time_stream)))
        #### Test making noise in map space.
        # First make the noise matrix by brute force.
        P = dirty_map.Pointing(("ra", "dec"), (ra, dec), map, 'nearest')
        map_noise_inv = sp.zeros((nf, nra, ndec, nf, nra, ndec), dtype=float)
        P_mat = P.get_matrix()
        tmp_map_noise_inv = al.partial_dot(noise_inv,
                                           P_mat)
        tmp_map_noise_inv = al.partial_dot(P_mat.mat_transpose(), 
                                           tmp_map_noise_inv)
        # I mess up the meta data by doing this, but rotate the axes so they
        # are in the desired order.
        tmp_map_noise_inv = sp.rollaxis(tmp_map_noise_inv, 2, 0)
        # Now use fast methods.
        map_noise_inv = sp.zeros((nf, nra, ndec, nf, nra, ndec), dtype=float)
        map_noise_inv = al.make_mat(map_noise_inv, axis_names=('freq', 'ra', 
            'dec', 'freq', 'ra', 'dec'), row_axes=(0, 1, 2), 
            col_axes=(3, 4, 5))
        for ii in xrange(nf):
            Noise.update_map_noise(P, ii, map_noise_inv[ii,:,:,:,:,:])
        self.assertTrue(sp.allclose(map_noise_inv, tmp_map_noise_inv))


if __name__ == '__main__' :
    unittest.main()
