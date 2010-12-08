"""Unit tests for subtract_map_data.py"""

import unittest

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
from core import fits_map, fitsGBT
import subtract_map_data as smd
import rotate_pol


test_file = 'testfile_GBTfits.fits'
test_map = 'testfile_map.fits'
n_pointings = 10 # Know property of test_file.  Per scan.

class TestSubMap(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.blocks = Reader.read((),())
        self.Map = fits_map.read(test_map, feedback=0)
        for Data in self.blocks :
            rotate_pol.rotate(Data, (1,))

    def test_runs(self) :
        for Data in self.blocks :
            smd.sub_map(Data, self.Map)

    def test_checks_pols(self) :
        Data = self.blocks[2]
        # This should work (for future compatibility with multiple
        # polarizations).
        smd.sub_map(Data, (self.Map,))
        # This should fail since the polarizations have to match.
        Data.field['CRVAL4'][0] = -5
        self.assertRaises(ce.DataError, smd.sub_map, Data, self.Map)

    def test_rigged_pointing(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        Map = self.Map
        # Set all data = (f + cal_ind)*time_ind
        Data.data[:,:,:,:] = (sp.arange(1,11)[:,sp.newaxis,sp.newaxis,sp.newaxis]
                              *(Data.freq/1e6+sp.arange(2).reshape((1,1,2,1))))
        Map.calc_axes()
        Map.data[:,:,:] = 0.0
        # Set 10 pixels to match data (except for cal_ind part).
        Map.data[range(10), range(10), :] = (sp.arange(1,11)[:,sp.newaxis]
                                             * Map.freq/1e6)
        # Rig the pointing to point to those 10 pixels.
        def rigged_pointing() :
            Data.ra = Map.long[range(10)]
            Data.dec = Map.lat[range(10)]
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, Map)
        # Now data should be just time_ind*cal_ind, within 2.0 MHz.
        Data.data /= sp.arange(1,11)[:,sp.newaxis,sp.newaxis,sp.newaxis]
        self.assertTrue(sp.allclose(Data.data[:,:,0,:], 0.0, atol=1.05))
        self.assertTrue(sp.allclose(Data.data[:,:,1,:], 1.0, atol=1.05))

    def test_off_map(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        Map = self.Map
        Map.data[:,:,:] = 0.0
        Data.data[:,:,:,:] = 0.0
        Map.calc_axes()
        # Rig the pointing but put one off the map.
        def rigged_pointing() :
            Data.ra = Map.long[range(10)]
            Data.dec = Map.lat[range(10)]
            Data.ra[3] = Data.ra[3] - 8.0
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, Map)
        self.assertTrue(sp.alltrue(ma.getmaskarray(Data.data[3,:,:,:])))
        self.assertTrue(sp.alltrue(sp.logical_not(
                    ma.getmaskarray((Data.data[[0,1,2,4,5,6,7,8,9],:,:,:])))))

    def tearDown(self) :
        del self.blocks
        del self.Map



if __name__ == '__main__' :
    unittest.main()
