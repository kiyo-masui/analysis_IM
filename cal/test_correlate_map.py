import unittest
import os

import numpy as np

test_file = 'testdata/testfile_guppi_rotated.fits'

class TestGetCorrelation(unittest.TestCase):

    def setUp(self):
        Reader = fitsGBT.Reader(this_test_file, feedback=0)
        Data = Reader.read(0, 0)
        Data.data[...] = 0
        Data.calc_freq()
        freq = Data.freq
        n_chan = len(freq)
        centre = freq[n_chan//2]
        delta = -abs(np.mean(np.diff(freq)))
        map = sp.zeros((n_chan, 10, 11))
        map = algebra.make_vect(map, axis_names=('freq', 'ra', 'dec'))
        map.set_axis_info('freq', centre, delta)
        map.set_axis_info('ra', 325.6, -0.2)
        map.set_axis_info('dec', 0, 0.2)
        self.Data = Data
        self.map = map



if __name__ == '__main__' :
    unittest.main()
