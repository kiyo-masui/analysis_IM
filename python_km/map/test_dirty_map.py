"""Unit tests for dirty map maker."""

import os
import glob

import unittest
import scipy as sp
import numpy.ma as ma

import dirty_map
import tools

class TestAddingData(unittest.TestCase) :
    
    def setUp(self) : 
        self.shape = (5, 5, 20) 
        self.ntime = 15
        self.map = ma.zeros(self.shape, dtype=float)
        self.noise_inv = sp.zeros(self.shape, dtype=float)
        self.data = ma.ones((self.ntime, self.shape[-1]))
        self.centre = (90.0,0.0,100.0)
        self.spacing = 1.0
        self.ra = 0.2*(sp.arange(self.ntime)-self.ntime/2.0) + self.centre[0]
        self.dec = 0.2*(sp.arange(self.ntime)-self.ntime/2.0) + self.centre[1]
        self.ra_inds = tools.calc_inds(self.ra, self.centre[0], self.shape[0],
                                    self.spacing)
        self.dec_inds = tools.calc_inds(self.dec, self.centre[1],
                                        self.shape[1], self.spacing)
        
    def test_adds_data(self) :
        dirty_map.add_data_2_map(self.data, self.ra_inds,  self.dec_inds, 
                                 self.map, self.noise_inv)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1])

    def test_mask_data(self) :
        self.data[5,8] = ma.masked 
        dirty_map.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.noise_inv)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1]-1)

    def test_off_map(self) :
        self.ra_inds[3] = 10 
        self.dec_inds[6] = -1
        dirty_map.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.noise_inv)
        self.assertAlmostEqual(ma.sum(self.map),
                               self.ntime*self.shape[-1]-2*self.shape[2])

    def test_checks_shapes(self) :
        self.assertRaises(ValueError, dirty_map.add_data_2_map, 
                          self.data, self.ra, sp.arange(self.ntime+1), 
                          self.map, self.noise_inv)
        self.assertRaises(ValueError, dirty_map.add_data_2_map, 
                          self.data, self.ra, self.dec, self.map[:,:,0:-1], 
                          self.noise_inv)

    def test_no_noise(self) :
        dirty_map.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, None)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1])
        
    def tearDown(self) :
        del self.data 
        del self.map

class TestDisjointScans(unittest.TestCase) :
    
    def setUp(self) :
        self.data = ma.ones((10, 20))*5 
        self.ra_inds = sp.array([5, 5, 5, 9, 7, 8, 4, 5, 8, 8], dtype=int)
        self.dec_inds = sp.array([2, 2, 2, 5, 3, 1, 3, 3, 1, 4], dtype=int)
        self.pix_counts = sp.zeros((10, 20), dtype=int)
        self.noise_inv = sp.zeros((10,6,10,6,20))

    def test_pix_counts_basic_cases(self) :
        pixels = dirty_map.pixel_counts(self.data, self.ra_inds, self.dec_inds,
                                      self.pix_counts)
        self.assertTrue(sp.alltrue(self.pix_counts[pixels.index((8,1))] == 2))
        self.assertTrue(sp.alltrue(self.pix_counts[pixels.index((5,2))] == 3))
        self.assertTrue(sp.alltrue(self.pix_counts[pixels.index((9,5))] == 1))
        self.assertTrue(sp.alltrue(self.pix_counts[pixels.index((7,3))] == 1))
        self.assertTrue(sp.alltrue(self.pix_counts[pixels.index((8,4))] == 1))

    def test_pix_counts_masked(self) :
        self.data[2, 15] = ma.masked
        self.data[4, 8] = ma.masked
        pixels = dirty_map.pixel_counts(self.data, self.ra_inds, self.dec_inds,
                                      self.pix_counts)
        self.assertEqual(self.pix_counts[pixels.index((5,2)),15], 2)
        self.assertEqual(self.pix_counts[pixels.index((5,2)),8], 3)
        self.assertEqual(self.pix_counts[pixels.index((7,3)),8], 0)

    def test_off_map(self) :
        pixels = dirty_map.pixel_counts(self.data, self.ra_inds, self.dec_inds,
                                      self.pix_counts, map_shape=(7,7))
        self.assertTrue(not (8,1) in pixels)
        self.assertTrue(not (9,5) in pixels)
        self.assertTrue(not (8,4) in pixels)
        self.assertTrue(not (7,3) in pixels)
        self.assertTrue(sp.alltrue(self.pix_counts[pixels.index((5,2))] == 3))

    def test_pix_counts_raises(self) :
        self.ra_inds = self.ra_inds[:-1]
        self.assertRaises(ValueError, dirty_map.pixel_counts, self.data,
                          self.ra_inds, self.dec_inds, self.pix_counts)
    
    def test_pix_counts_raises2(self) :
        self.pix_counts = self.pix_counts[:,:-1]
        self.assertRaises(ValueError, dirty_map.pixel_counts, self.data, 
                          self.ra_inds, self.dec_inds, self.pix_counts)
    
    def test_add_scan_noise(self) :
        pixels = dirty_map.pixel_counts(self.data, self.ra_inds,
                                        self.dec_inds, self.pix_counts)
        dirty_map.add_scan_noise(pixels, self.pix_counts, 2.0, self.noise_inv)
        self.assertTrue(sp.allclose(self.noise_inv[5,2,5,2,:], 
                                    (3.0-9.0/10)/2.0))
        self.assertTrue(sp.allclose(self.noise_inv[5,3,5,3,:], 
                                    (1.0-1.0/10)/2.0))
        self.assertTrue(sp.allclose(self.noise_inv[8,1,5,2,:], 
                                    (-6.0/10)/2.0))
        self.assertTrue(sp.allclose(self.noise_inv[5,2,8,1,:], 
                                    (-6.0/10)/2.0))

    def test_add_scan_noise_masked(self) :
        self.data[2, 15] = ma.masked
        self.data[4, 8] = ma.masked
        pixels = dirty_map.pixel_counts(self.data, self.ra_inds, self.dec_inds,
                                      self.pix_counts)
        dirty_map.add_scan_noise(pixels, self.pix_counts, 2.0, self.noise_inv)
        self.assertEqual(self.noise_inv[5,2,5,2,15], (2.0-4.0/9.0)/2.0)
        self.assertEqual(self.noise_inv[5,2,5,2,4], (3.0-9.0/10.0)/2.0)
        self.assertTrue(sp.allclose(self.noise_inv[:,:,7,3,8], 0.0))
        self.assertEqual(self.noise_inv[5,2,7,3,15], (-2.0/9.0)/2.0)

    def test_add_raises(self) :
        pixels = dirty_map.pixel_counts(self.data, self.ra_inds, self.dec_inds,
                                      self.pix_counts)
        self.noise_inv = self.noise_inv[1,...]
        self.assertRaises(ValueError, dirty_map.add_scan_noise, pixels, 
                          self.pix_counts, 2.0, self.noise_inv)

    def tearDown(self) :
        del self.data
        del self.ra_inds
        del self.dec_inds
        del self.pix_counts
        del self.noise_inv

class TestRuns(unittest.TestCase) :

    def setUp(self) :
        self.pars = { "dm_output_root" : "./testoutput_",
                      'dm_file_middles' : ("testfile_GBTfits",),
                      'dm_input_end' : ".fits",
                      'dm_scans' : (),
                      'dm_IFs' : (0,),
                      'dm_field_centre' : (325.0, 0.0),
                      'dm_map_shape' : (5, 5),
                      'dm_pixel_spacing' : 0.5
                    }
    
    def test_gridder(self) :
        pars = { 'dm_noise_model' : 'grid'}
        pars.update(self.pars)
        Maker = dirty_map.DirtyMapMaker(pars, 0)
        Maker.execute()

    def test_diag(self) :
        pars = { 'dm_noise_model' : 'diag_file'}
        pars.update(self.pars)
        Maker = dirty_map.DirtyMapMaker(pars, 0)
        Maker.execute()

    def test_full(self) :
        pars = { 'dm_noise_model' : 'disjoint_scans'}
        pars.update(self.pars)
        Maker = dirty_map.DirtyMapMaker(pars, 0)
        Maker.execute()

    def tearDown(self) :
        files = glob.glob('testoutput_*')
        for f in files :
            os.remove(f)

if __name__ == '__main__' :
    unittest.main()
