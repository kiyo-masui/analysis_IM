"""Unit tests for dirty map maker."""

import os
import glob

import unittest
import scipy as sp
import numpy.ma as ma

import core.algebra as al
import new_dirty_map as dirty_map
import tools

nscans = 10
nt = 50 * nscans
nf = 10

nra = 20
ndec = 20

# XXX: Prefilter that makes a time stream vect will deal with out of bounds
# pointings, adding infinite noise.  The pointing should not have to deal with
# this.

def make_test_data():
    # Make soem fake pointing data.  This isn't entirely consitant with
    # possible geometies on the sky, but that probably doens't matter for
    # testing.
    # 5 degress, drifting, centred at 21.
    ra = sp.arange(-nt/2.0, nt/2.0, dtype=float)/nt*5 + 21  
    # Make dec scan back and forth, 5 degrees, centred at 0.
    dec = sp.zeros(nt)
    for ii in range(0, nscans, 2):
        dec[ii*nt/nscans:(ii+1)*nt/nscans] = sp.arange(-nt/nscans/2.0,
                                                       nt/nscans/2.0, dtype=float)
    for ii in range(1, nscans, 2):
        dec[ii*nt/nscans:(ii+1)*nt/nscans] = -sp.arange(-nt/nscans/2.0, 
                                                        nt/nscans/2.0,
                                                        dtype=float)
    dec *= 5.0/(nt/nscans)
    # Constant el.
    el = sp.zeros(nt) + 40.0
    # Scanning az.
    # XXX no cross linking.
    az = sp.copy(dec)
    # Now make some data that has stuff that both correlates on the sky and in
    # az-el.
    time_stream = sp.zeros((nf, nt))
    time_stream += sp.sin(ra*3)*sp.cos(dec*2) # The sky map.
    time_stream *= sp.arange(nf)[:,None]*0.1 + 2.0
    time_stream += sp.sin(az/2) # The ground spill.
    time_stream = al.make_vect(time_stream, ('freq', 'time'))

    # Now make the map that we are going to accumulate this into.
    map = sp.zeros((nf, nra, ndec))
    map = al.make_vect(map, ('freq', 'ra', 'dec'))
    map.set_axis_info('freq', 800e6, 1e6)
    map.set_axis_info('ra', 21, 6.0/nra)
    map.set_axis_info('dec', 0, 6.0/nra)

    return map, time_stream, ra, dec, az, el



class TestPointing(unittest.TestCase):
    
    def test_get_matrix(self):
        map, time_stream, ra, dec, az, el = make_test_data()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), map, 'nearest')
        pointing_matrix = P.get_matrix()
        self.assertTrue(pointing_matrix.axes == ('time', 'ra', 'dec'))
        self.assertEqual(pointing_matrix.shape[0], ra.shape[0])
        self.assertEqual(map.shape[1], pointing_matrix.shape[1])
        self.assertEqual(map.shape[2], pointing_matrix.shape[2])
    
    def test_apply_time_axis_nearest(self):
        map, time_stream, ra, dec, az, el = make_test_data()
        P = dirty_map.Pointing(('ra', 'dec'), (ra, dec), map, 'nearest')
        pointing_matrix = P.get_matrix()
        gridded_data_mult = al.dot(pointing_matrix, time_stream)
        gridded_data_fast = P.apply_to_time_axis(time_stream)
        self.assertTrue(sp.allclose(gridded_data_mult, gridded_data_fast))
        



if __name__ == '__main__' :
    unittest.main()
