"""Unit tests for pol_cal.py"""

import unittest
import os

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
from core import fits_map, fitsGBT
import pol_cal
import rotate_pol

test_file = 'testfile_GBTfits.fits'
test_mueler_file = 'fname'
n_pointings = 10 # Known property of test_file.  Per scan.


# Tests for mueller.
class TestMueller(unittest.TestCase) :
    
    def test_runs(self) :
        mueler = pol_cal.mueller()


# Tests for the calibrate_pol function
class TestCalPol(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.blocks = Reader.read((),())
        self.mueler = pol_cal.mueller()
 	print self.mueler   

    def test_runs(self) :
        for Data in self.blocks :
            rotate_pol.rotate(Data, (1,2,3,4))
            pol_cal.calibrate_pol(Data, self.mueler)
        
    def test_raises(self) :
        Data = self.blocks[0]
        rotate_pol.rotate(Data, (1,))
        self.assertRaises(ce.DataError, pol_cal.calibrate_pol, Data,self.mueler)
        Data = self.blocks[1]
        rotate_pol.rotate(Data, (1,2,3,4))
        Data.field['CRVAL4'][0] = 3
        Data.field['CRVAL4'][2] = 1
        self.assertRaises(ce.DataError, pol_cal.calibrate_pol, Data,
                          self.mueler)

    def test_gets_answer_right(self) :
        pass
   
    def tearDown(self) :
        del self.blocks
        del self.mueler


# Test that the whole thing runs as a module.
class TestModule(unittest.TestCase) :
    
    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        blocks = Reader.read((),())
        for Data in blocks :
            rotate_pol.rotate(Data, (1,2,3,4))
        Writer = fitsGBT.Writer(blocks, feedback=0)
        Writer.write('test_rot.testout.fits')

        self.params = {'pc_file_middles' : ("test",),
                       'pc_input_end' : "_rot.testout.fits",
                       'pc_output_end' : "_polcal.testout.fits",
                       'pc_mueler_file' : test_mueler_file
                       }

    def test_runs(self) :
        pol_cal.Calibrate(self.params, feedback=0).execute()

    def tearDown(self) :
        os.remove('test_rot.testout.fits')
        os.remove('test_polcal.testout.fits')
        os.remove('params.ini')
        



if __name__ == '__main__' :
    unittest.main()
