"""Unit tests for subtract_map_data.py"""

import unittest
import os
import cPickle

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
from core import algebra, fitsGBT
import subtract_map_data as smd
import rotate_pol


test_file = 'testfile_GBTfits.fits'
n_pointings = 10 # Known property of test_file.  Per scan.

class TestSubMap(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.blocks = Reader.read((),())
        b = self.blocks[0]
        b.calc_freq()
        map = sp.zeros((150, 10, 11))
        map = algebra.make_vect(map, axis_names=('freq', 'ra', 'dec'))
        map.set_axis_info('freq', 700.0e6, -2.0e6)
        map.set_axis_info('ra', 325.6, -0.2)
        map.set_axis_info('dec', 0, 0.2)
        self.map = map

        for Data in self.blocks :
            rotate_pol.rotate(Data, (1, 2, 3, 4))

    def test_runs(self) :
        for Data in self.blocks :
            smd.sub_map(Data, self.map)

    def test_checks_pols(self) :
        Data = self.blocks[2]
        # This should work (for future compatibility with multiple
        # polarizations).
        smd.sub_map(Data, (self.map,), pols=(0,))
        # But if I supply, multiple maps, it had better be one per
        # polarization.
        self.assertRaises(ValueError, smd.sub_map, Data, (self.map, self.map),
                          pols=(0,2,3))

    def test_rigged_pointing(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        map = self.map
        # Set all data = (f + cal_ind)*time_ind
        Data.data[:,:,:,:] = (sp.arange(-4.5, 5)
                              [:,sp.newaxis,sp.newaxis,sp.newaxis]
                              *(Data.freq/100e6))
        Data.data[...] -= sp.mean(Data.data, 0)
        Data.data[...] += (sp.arange(6,8).reshape((1,1,2,1)) * (Data.freq/100e6) 
                           * sp.arange(-4.5, 5).reshape((10, 1, 1, 1)))
        map[:,:,:] = 0.0
        # Set 10 pixels to match data (except for cal_ind part).
        map[:, range(10), range(10)] = (sp.arange(-4.5, 5)[None,:]
                                        * map.get_axis('freq')[:,None]/100e6)
        # We should be completely insensitive to the map mean.  Th following
        # should have no effect.
        map[...] += 0.352*map.get_axis('freq')[:, None, None]/800.0e7
        # Rig the pointing to point to those 10 pixels.
        def rigged_pointing() :
            Data.ra = map.get_axis('ra')[range(10)]
            Data.dec = map.get_axis('dec')[range(10)]
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, map)
        # Now data should be just f*time_ind*(cal_ind+6), within 2.0 MHz/2.
        Data.data /= sp.arange(-4.5, 5)[:,sp.newaxis,sp.newaxis,sp.newaxis]
        Data.data /= Data.freq/100e6
        # Relative tol of 1/700, is the frequency bin width.
        self.assertTrue(sp.allclose(Data.data[:,:,0,:], 6.0, rtol=1.0/700))
        self.assertTrue(sp.allclose(Data.data[:,:,1,:], 7.0, rtol=1.0/700))

    def test_correlate(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        map = self.map
        gain = 3.45
        const = 2.14
        # Set all data = gain*(cos(time_ind)).
        Data.data[:,:,:,:] = gain*sp.cos(sp.arange(1,11)
                                    [:,sp.newaxis,sp.newaxis,sp.newaxis])
        # Explicitly set time mean to something known.
        Data.data -= ma.mean(Data.data, 0)
        Data.data += gain*const*Data.freq/800.0e6
        # Now the Map.
        map[:,:,:] = 0.0
        # Set 10 pixels to match cos part of data.
        map[:, range(10), range(10)] = (
                    sp.cos(sp.arange(1,11)[None, :]))
        map[:, range(10), range(10)] -= ma.mean(
            map[:, range(10), range(10)], 1)[:, None]
        # Give Map a mean to test things out. Should really have no effect.
        map[...] += 0.352*map.get_axis('freq')[:, None, None]/800.0e6
        # Rig the pointing to point to those 10 pixels.
        def rigged_pointing() :
            Data.ra = map.get_axis('ra')[range(10)]
            Data.dec = map.get_axis('dec')[range(10)]
        Data.calc_pointing = rigged_pointing
        solved_gains = smd.sub_map(Data, map, correlate=True)
        # Now data should be just be gain*const*f, within machine precision.
        Data.data /= gain*Data.freq/800.0e6
        self.assertTrue(sp.allclose(Data.data[:,:,:,:], const))
        self.assertTrue(sp.allclose(solved_gains, gain))
    
    def test_off_map(self) :
        Data = self.blocks[0]
        Data.calc_freq()
        map = self.map
        map[:,:,:] = 0.0
        Data.data[:,:,:,:] = 0.0
        # Rig the pointing but put one off the map.
        def rigged_pointing() :
            Data.ra = map.get_axis('ra')[range(10)]
            Data.dec = map.get_axis('dec')[range(10)]
            Data.ra[3] = Data.ra[3] - 8.0
        Data.calc_pointing = rigged_pointing
        smd.sub_map(Data, map)
        self.assertTrue(sp.alltrue(ma.getmaskarray(Data.data[3,:,:,:])))
        self.assertTrue(sp.alltrue(sp.logical_not(
                    ma.getmaskarray((Data.data[[0,1,2,4,5,6,7,8,9],:,:,:])))))

    def tearDown(self) :
        del self.blocks
        del self.map

class TestModule(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        blocks = Reader.read((),())
        for Data in blocks :
            rotate_pol.rotate(Data, (1,2,3,4))
        Writer = fitsGBT.Writer(blocks, feedback=0)
        Writer.write('test_rot.testout.fits')
        
        map = sp.zeros((150, 10, 11))
        map = algebra.make_vect(map, axis_names=('freq', 'ra', 'dec'))
        map.set_axis_info('freq', 700.0e6, -2.0e6)
        map.set_axis_info('ra', 325.6, -0.2)
        map.set_axis_info('dec', 0, 0.2)
        algebra.save('test_inmap.npy', map)

        self.params = {'sm_file_middles' : ("test",),
                       'sm_input_end' : "_rot.testout.fits",
                       'sm_output_end' : "_sub.testout.fits",
                       'sm_solve_for_gain' : True,
                       'sm_gain_output_end' : '_gain.pickle',
                       'sm_map_file' : 'test_inmap.npy'
                       }

    def test_history(self) :
        smd.Subtract(self.params, feedback=0).execute()
        Data = fitsGBT.Reader('test_sub.testout.fits', feedback=0).read(0,0)
        self.assertTrue(Data.history.has_key('003: Subtracted map from data.'))

    def tearDown(self) :
        os.remove('test_inmap.npy')
        os.remove('test_rot.testout.fits')
        os.remove('test_sub.testout.fits')
        os.remove('params.ini')
        os.remove('test_gain.pickle')
        
if __name__ == '__main__' :
    unittest.main()
