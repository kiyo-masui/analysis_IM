"""Unit tests for subtract_map_data.py"""

import unittest
import os
import cPickle

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
        Data.data[:,:,:,:] = (sp.arange(-4.5, 5)
                              [:,sp.newaxis,sp.newaxis,sp.newaxis]
                              *(Data.freq/1e6+sp.arange(2).reshape((1,1,2,1))))
        Map.calc_axes()
        Map.data[:,:,:] = 0.0
        # Set 10 pixels to match data (except for cal_ind part).
        Map.data[range(10), range(10), :] = (sp.arange(-4.5, 5)[:,sp.newaxis]
                                             * Map.freq/1e6)
        # Rig the pointing to point to those 10 pixels.
        def rigged_pointing() :
            Data.ra = Map.long[range(10)]
            Data.dec = Map.lat[range(10)]
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, Map)
        # Now data should be just time_ind*cal_ind, within 2.0 MHz/2.
        Data.data /= sp.arange(-4.5,5)[:,sp.newaxis,sp.newaxis,sp.newaxis]
        self.assertTrue(sp.allclose(Data.data[:,:,0,:], 0.0, atol=1.05))
        self.assertTrue(sp.allclose(Data.data[:,:,1,:], 1.0, atol=1.05))

    def test_correlate(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        Map = self.Map
        gain = 3.45
        const = 2.14
        # Set all data = gain*(cos(time_ind)).
        Data.data[:,:,:,:] = gain*sp.cos(sp.arange(1,11)
                                    [:,sp.newaxis,sp.newaxis,sp.newaxis])
        # Explicitly set time mean to something known.
        Data.data -= ma.mean(Data.data, 0)
        Data.data += gain*const*Data.freq/1e6
        # Now the Map.
        Map.calc_axes()
        Map.data[:,:,:] = 0.0
        # Set 10 pixels to match cos part of data.
        Map.data[range(10), range(10), :] = (
                    sp.cos(sp.arange(1,11)[:,sp.newaxis]))
        Map.data[range(10), range(10), :] -= ma.mean(Map.data[range(10), 
                    range(10), :], 0)
        # Give Map a mean o test things out.
        Map.data += 2.5
        # Rig the pointing to point to those 10 pixels.
        def rigged_pointing() :
            Data.ra = Map.long[range(10)]
            Data.dec = Map.lat[range(10)]
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, Map, correlate=True)
        # Subtract the Map mean part.
        Data.data += gain*2.5
        # Now data should be just be gain*const*f, within machine precision.
        Data.data /= gain*Data.freq/1e6
        self.assertTrue(sp.allclose(Data.data[:,:,:,:], const))
    
    def test_masked_map(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        Map = self.Map
        # Set all data = f*time_ind
        Data.data[:,:,:,:] = (Data.freq/1e6*sp.arange(-4.5, 5)
                              [:,sp.newaxis,sp.newaxis,sp.newaxis])
        Data.data[3,:,:,:] = ma.masked
        Map.calc_axes()
        Map.data[:,:,:] = 0.0
        # Set 10 pixels to match data.
        Map.data[range(10), range(10), :] = (sp.arange(-4.5, 5)[:,sp.newaxis]
                                             * Map.freq/1e6)
        # Mask one of these pixels. Masking 4 and 5 preserves the 0 mean.
        Map.data[5,5,:] = 99999999.
        Map.data[5,5,:] = ma.masked
        Map.data[4,4,:] = ma.masked
        # Rig the pointing to point to those 10 pixels.
        def rigged_pointing() :
            Data.ra = Map.long[range(10)]
            Data.dec = Map.lat[range(10)]
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, Map)
        self.assertTrue(sp.alltrue(Data.data[5,:,:,:].mask))
        self.assertTrue(sp.alltrue(Data.data[4,:,:,:].mask))
        self.assertTrue(sp.alltrue(Data.data[3,:,:,:].mask))
        Data.data[5,:,:,:] = 0
        Data.data[4,:,:,:] = 0
        Data.data[3,:,:,:] = 0
        # Now data should be just time_ind*cal_ind, within 2.0 MHz/2.
        Data.data /= sp.arange(-4.5,5)[:,sp.newaxis,sp.newaxis,sp.newaxis]
        self.assertTrue(sp.allclose(Data.data, 0.0, atol=1.05))
   
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

class TestModule(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        blocks = Reader.read((),())
        for Data in blocks :
            rotate_pol.rotate(Data, (1,))
        Writer = fitsGBT.Writer(blocks, feedback=0)
        Writer.write('test_rot.testout.fits')

        self.params = {'sm_file_middles' : ("test",),
                       'sm_input_end' : "_rot.testout.fits",
                       'sm_output_end' : "_sub.testout.fits",
                       'sm_solve_for_gain' : True,
                       'sm_gain_output_end' : '_gain.pickle'
                       }

    def test_history(self) :
        smd.Subtract(self.params, feedback=0).execute()
        Data = fitsGBT.Reader('test_sub.testout.fits', feedback=0).read(0,0)
        self.assertTrue(Data.history.has_key('003: Subtracted map from data.'))

    def tearDown(self) :
        os.remove('test_rot.testout.fits')
        os.remove('test_sub.testout.fits')
        os.remove('params.ini')
        os.remove('test_gain.pickle')
        



if __name__ == '__main__' :
    unittest.main()
