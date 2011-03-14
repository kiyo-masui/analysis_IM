"""Unit tests for dirty map maker."""

import os
import glob

import unittest
import scipy as sp
import numpy.ma as ma

import clean_map
import dirty_map

class TestRuns(unittest.TestCase) :

    def setUp(self) :
        self.pars = { "dm_output_root" : "./testoutput_",
                      'dm_file_middles' : ("testfile_GBTfits",),
                      'dm_input_end' : ".fits",
                      'dm_scans' : (),
                      'dm_IFs' : (0,),
                      'dm_field_centre' : (325.0, 0.0),
                      'dm_map_shape' : (5, 5),
                      'dm_pixel_spacing' : 0.5,
                      # Clean map params:
                      'cm_input_root' :  "./testoutput_",
                      'cm_output_root' : "./testoutput_",
                      'cm_dirty_map_files' : ('dirty_map_XX.npy',),
                      'cm_noise_inverse_files' : ('noise_inv_XX.npy',)
                      }
    
    def test_gridder(self) :
        pars = { 'dm_noise_model' : 'grid'}
        pars.update(self.pars)
        Maker = dirty_map.DirtyMapMaker(pars, 0)
        Maker.execute()
        Maker = clean_map.CleanMapMaker(pars, 0)
        Maker.execute()

    def test_full(self) :
        pars = { 'dm_noise_model' : 'disjoint_scans'}
        pars.update(self.pars)
        Maker = dirty_map.DirtyMapMaker(pars, 0)
        Maker.execute()
        Maker = clean_map.CleanMapMaker(pars, 0)
        Maker.execute()

    def tearDown(self) :
        files = glob.glob('testoutput_*')
        for f in files :
            os.remove(f)

if __name__ == '__main__' :
    unittest.main()
