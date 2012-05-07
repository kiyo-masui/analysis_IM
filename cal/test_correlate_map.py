import unittest
import os

import numpy as np
from scipy import interpolate
from numpy import random
import numpy.ma as ma

from core import fitsGBT, algebra
import correlate_map

test_file = 'testdata/testfile_guppi_rotated.fits'

class TestGetCorrelation(unittest.TestCase):

    def setUp(self):
        Reader = fitsGBT.Reader(test_file, feedback=0)
        Data = Reader.read(0, 0)
        Data.data[...] = 0
        Data.calc_freq()
        freq = Data.freq
        n_chan = len(freq)
        centre = freq[n_chan//2]
        delta = -abs(np.mean(np.diff(freq)))
        map = np.zeros((n_chan, 20, 10))
        map = algebra.make_vect(map, axis_names=('freq', 'ra', 'dec'))
        map.set_axis_info('freq', centre, delta)
        map.set_axis_info('ra', 0, -0.1)
        map.set_axis_info('dec', 0, 0.1)
        self.Data = Data
        self.map = map

    def test_perfect_correlation(self):
        map = self.map
        Data = self.Data
        ra_map = map.get_axis('ra')
        dec_map = map.get_axis('dec')
        freq = map.get_axis('freq')
        n_chan = len(freq)
        n_time = Data.dims[0]
        # Make the ra-dec map a quadratic so it's perfectly interpolable.
        map[...] = (ra_map[None,:,None] + 3*dec_map[None,None,:] + 6.)**2
        map[...] *= np.arange(n_chan)[:,None,None] + 10.
        # Make the time domain ra and dec.
        ra = ((random.rand(n_time) - 0.5)
              * len(ra_map) * np.diff(ra_map)[0] * 0.8)
        dec = ((random.rand(n_time) - 0.5)
               * len(dec_map) * np.diff(dec_map)[0] * 0.8)
        # Set the data equal to the map.
        for ii in range(n_time):
            Data.data[ii,0,0,:] = map.slice_interpolate((1, 2), 
                                        (ra[ii], dec[ii]), kind='cubic')
        # Mask things out all over the place.
        Data.data[:,0,0,12] = ma.masked
        Data.data[80,0,0,:] = ma.masked
        Data.data[8:20,0,0,3] = ma.masked
        # Rig the pointing of the Data object.
        def rigged_pointing() :
            Data.ra = ra
            Data.dec = dec
        Data.calc_pointing = rigged_pointing
        corr, norm = correlate_map.get_correlation(Data, (map,), 'cubic', 10)
        # Since there is no noise, the correlation should be exactly unity.
        self.assertTrue(np.allclose(corr, norm))




if __name__ == '__main__' :
    unittest.main()
