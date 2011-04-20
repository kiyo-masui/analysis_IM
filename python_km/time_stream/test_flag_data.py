"""Unit tests for flag_data.py"""

import unittest

import scipy as sp
import numpy.ma as ma
import numpy.random as rand
import matplotlib.pyplot as plt

import kiyopy.custom_exceptions as ce
import flag_data
import hanning
from core import data_block, fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestFlagData(unittest.TestCase) :

    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,0)

    def test_different_lenght(self) :
        self.Data.data = self.Data.data[:,:,:,:1024]
        self.Data.verify()
        ind = (6,1,1,125)
        ind2 = (3,2,0,987)
        self.Data.data[ind] = 0.2
        self.Data.data[ind2] = 0.2
        flag_data.apply_cuts(self.Data, -1, 5.0, 0, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)

    def test_different_lenght2(self) :
        self.Data.data = ma.copy(ma.concatenate((ma.copy(self.Data.data),
                                                 self.Data.data), -1))
        self.Data.verify()
        ind = (6,1,1,362)
        ind2 = (3,2,0,1942)
        self.Data.data[ind] = 0.2
        self.Data.data[ind2] = 0.2
        flag_data.apply_cuts(self.Data, -1, 5.0, 0, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)

    def test_pol_cut(self) :
        ind = (6,1,1,676)
        ind2 = (3,2,0,245)
        self.Data.data[ind] = 0.2
        self.Data.data[ind2] = 0.2
        flag_data.apply_cuts(self.Data, -1, 5.0, 0, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)
        self.assertTrue(float(ma.count_masked(self.Data.data)) / 
                        float(self.Data.data.size) < 0.1)

    def test_pol_cut_IQUV(self) :
        self.Data.field['CRVAL4'] = sp.array([1, 2, 3, 4], dtype=int)
        self.Data.data[0::2,:,:,:] = 1
        self.Data.data[1::2,:,:,:] = 2
        ind = (6,3,1,676)
        ind2 = (3,2,0,245)
        self.Data.data[ind] = 8
        self.Data.data[ind2] = 10
        flag_data.apply_cuts(self.Data, -1, 5.0, 0, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)
        self.assertTrue(float(ma.count_masked(self.Data.data)) / 
                        float(self.Data.data.size) < 0.1)

    def test_flags_crazy_data(self) :
        ind = (5,0,0,1345)
        ind2 = (2,3,1,425)
        self.Data.data[ind] = 50
        self.Data.data[ind2] = 67
        flag_data.apply_cuts(self.Data, 5, -1)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)

    def test_flags_controled_data(self) :
        self.Data.data[0::2,:,:,:] = 1
        self.Data.data[1::2,:,:,:] = 2
        ind = (5,0,0,1345)
        ind2 = (2,3,1,425)
        ind3 = (8,1,1,754)
        self.Data.data[ind] = 5
        self.Data.data[ind2] = 6
        self.Data.data[ind3] = 6
        flag_data.apply_cuts(self.Data, 5, -1)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)
        self.assertTrue(self.Data.data[ind3] is ma.masked)
        # 8 = 2 bad data * 4 polarizations
        self.assertTrue(ma.count_masked(self.Data.data) == 12)
    
    def test_pol_cal_off_controled(self) :
        self.Data.data[:,[0,3],1,:] = rand.normal(5, 0.1, (10, 2, 2048))
        self.Data.data[:,[0,3],0,:] = rand.normal(5.5, 0.1, (10, 2, 2048))
        self.Data.data[:,[1,2],1,:] = rand.normal(0, 0.1, (10, 2, 2048))
        self.Data.data[:,[1,2],0,:] = rand.normal(0.5, 0.1, (10, 2, 2048))
        ind = (5,1,1,1345)
        ind2 = (2,2,0,425)
        self.Data.data[ind] += 3.0
        self.Data.data[ind2] += 3.0
        # Threshold rediculously high because distribution of cross has a long
        # tail.
        flag_data.apply_cuts(self.Data, -1, 30.0, 0, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)
        # 8 = 2 bad data * 4 polarizations
        self.assertTrue(ma.count_masked(self.Data.data) == 8)

    def test_flattens_and_width(self) :
        self.Data.data[:,[0,3],1,:] = rand.normal(5, 0.1, (10, 2, 2048))
        self.Data.data[:,[0,3],0,:] = rand.normal(5.5, 0.1, (10, 2, 2048))
        self.Data.data[:,[1,2],1,:] = rand.normal(0, 0.1, (10, 2, 2048))
        self.Data.data[:,[1,2],0,:] = rand.normal(0.5, 0.1, (10, 2, 2048))
        # Give the cross pol some structure.
        self.Data.data[:,[1,2],0,:] += sp.sin(sp.arange(2048)/100.)
        ind = (3,1,1,543)
        ind2 = (7,2,0,765)
        self.Data.data[ind] += 3.0
        self.Data.data[ind2] += 3.0
        width = 2
        # Threshold rediculously high because distribution of cross has a long
        # tail.
        flag_data.apply_cuts(self.Data, -1, 30.0, width, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)
        # 8 = 2 bad data * 4 polarizations
        self.assertTrue(ma.count_masked(self.Data.data) == 8*(2*width+1))

    def test_does_nothing_negitive(self) :
        self.Data.data[0::2,:,:,:] = 1
        self.Data.data[1::2,:,:,:] = 2
        ind = (5,0,0,1345)
        ind2 = (2,3,1,425)
        ind3 = (3,1,1,634)
        self.Data.data[ind] = 20
        self.Data.data[ind2] = 60
        self.Data.data[ind3] = 43
        flag_data.apply_cuts(self.Data, -1, -1)
        self.assertTrue(ma.count_masked(self.Data.data) == 0)

    def tearDown(self) :
        del self.Data

class TestHanning(unittest.TestCase) :

    def setUp(self) :
        Reader = core.fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,1)
        self.Data.verify()
        self.Data_copy = copy.deepcopy(self.Data)

    def test_hanning_data_changed(self) :
        """Copy the data, see that we did something."""
        hanning.hanning_smooth(self.Data)
        # For lack of anything better to test:
        self.Data.verify()
        # Make sure we actually did something.
        self.assertTrue(not ma.allclose(self.Data.data, self.Data_copy.data))
        # Make sure we didn't change other fields, like LST.
        self.assertTrue(sp.allclose(self.Data.field['LST'],
                        self.Data_copy.field['LST']))

    def test_hanning_cases(self) :
        data = self.Data.data
        data[:,:,:,:] = 1.
        data[5,2,1,631] = 4.
        data[7,3,1,853] = ma.masked
        hanning.hanning_smooth(self.Data)
        self.assertTrue(data[7,3,1,853] is ma.masked)
        self.assertTrue(data[7,3,1,852] is ma.masked)
        self.assertTrue(data[7,3,1,854] is ma.masked)
        self.assertAlmostEqual(data[5,2,1,631], 2.5)
        self.assertAlmostEqual(data[5,2,1,630], 1.75)
        self.assertAlmostEqual(data[5,2,1,632], 1.75)


if __name__ == '__main__' :
    unittest.main()

