"""unit tests for covariance.py."""

import unittest
import copy
import os
import glob

import scipy as sp

import covariance
from core import algebra
import covariance
import beam
import kiyopy.custom_exceptions as ce

class TestModule(unittest.TestCase) :
    """Tests that the whole module runs."""

    def setUp(self) :
        # Make a map that the Covariance gets it's axis info from
        map = sp.zeros((5, 10, 8), dtype=float)
        map = algebra.make_vect(map, axis_names=('freq', 'ra', 'dec'))
        map.set_axis_info("freq", 800.0e6, 2.0e6)
        dec = 10
        map.set_axis_info("ra", 215., 0.1/sp.cos(dec*sp.pi/180))
        map.set_axis_info("dec", dec, 0.1)
        algebra.save("tmp_mapfile.npy", map)
        # Set up some parameters to pass.
        self.params = {
            "cv_map_file" : "tmp_mapfile.npy",
            "cv_unit_system" : "deg-freq", 
            "cv_out_file_name" : "testout_covariance.npy"
            }
    
    def tearDown(self) :
        os.remove("tmp_mapfile.npy")    
        os.remove("tmp_mapfile.npy.meta")
        files = glob.glob('testout_*')
        for f in files :
            os.remove(f)

    def  test_runs_deg_freq(self) :
        #self.params["cv_unit_system"] = "deg-freq"
        covariance.Covariance(self.params, feedback=0).execute()
    
    def  test_runs_Mpc(self) :
        self.params["cv_unit_system"] = "Mpc/h-Mpc/h"
        covariance.Covariance(self.params, feedback=0).execute()

class TestIntegration(unittest.TestCase) :

    def setUp(self) :
        self.width = (5 + sp.arange(6, dtype=float))/5 # between, 1 and 2.
        self.frequencies = sp.arange(6, dtype=float)
        self.Beam = beam.GaussianBeam(self.width, self.frequencies)
        # 1/r^3 correlation function.
        n = 50
        r = sp.arange(n)
        r[0] = 0.1
        self.corr = 1.0/r[:, None]**2/r
        self.corr *= 1.0/n**3

if __name__ == '__main__' :
    unittest.main()
