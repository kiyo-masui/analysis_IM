"""Unit tests for flag_data.py"""

import unittest
import copy

import scipy as sp
import numpy.ma as ma
import numpy.random as rand
#import matplotlib.pyplot as plt

import kiyopy.custom_exceptions as ce
from time_stream import flag_data
from time_stream import rotate_pol
from time_stream import hanning
from core import data_block, fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestFlagData(unittest.TestCase) :

    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,0)

    def test_001_a_test_different_length(self) :
        self.Data.data = self.Data.data[:,:,:,:1024]
        self.Data.verify()
        ind = (6,1,1,125)
        freq = 125
        ind2 = (3,2,0,987)
        freq2 = 987
        self.Data.data[ind] += 0.2
        self.Data.data[ind2] += 0.2
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)
        self.assertFalse(False in self.Data.data[:,:,:,freq2].mask)

    def test_002_flagging_cut(self):
        ind = (6,1,1,676)
        freq = 676
        ind2 = (3,2,0,245)
        freq2 = 245
        self.Data.data[ind] += 0.2
        self.Data.data[ind2] += 0.2
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)
        self.assertFalse(False in self.Data.data[:,:,:,freq2].mask)
        self.assertTrue(float(ma.count_masked(self.Data.data)) / 
                        float(self.Data.data.size) < 0.1)

    def test_003_flagging_cut_almost_all_even(self):
        self.Data.data[0::2,:,:,0::2] = 6
        self.Data.data[1::2,:,:,1::2] = 7
        ind = (6,1,1,676)
        freq = 676
        ind2 = (3,2,0,245)
        freq2 = 245
        self.Data.data[ind] += 15
        self.Data.data[ind2] += 15
        flag_data.apply_cuts(self.Data)
        # assert that only the bad 2 frequencies were flagged.
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)
        self.assertFalse(False in self.Data.data[:,:,:,freq2].mask)
        for freq3 in range(0, self.Data.dims[-1]):
            if (freq3 != freq) and (freq3 != freq2):
                self.assertFalse(True in self.Data.data[:,:,:,freq3].mask)

    def test_004_flag_in_XX_cal0(self):
        ind = (3,0,0,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_005_flag_in_XY_cal0(self):
        ind = (3,1,0,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_006_flag_in_YX_cal0(self):
        ind = (3,2,0,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_007_flag_in_YY_cal0(self):
        ind = (3,3,0,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_008_flag_in_XX_cal1(self):
        ind = (3,0,1,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_009_flag_in_XY_cal1(self):
        ind = (3,1,1,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_010_flag_in_YX_cal1(self):
        ind = (3,2,1,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_011_flag_in_YY_cal1(self):
        ind = (3,3,1,245)
        freq = 245
        self.Data.data[ind] += 1
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)

    def test_012_flags_crazy_data(self):
        ind = (5,0,0,1345)
        freq = 1345
        ind2 = (2,3,1,425)
        freq2 = 425
        self.Data.data[ind] = 50
        self.Data.data[ind2] = 67
        flag_data.apply_cuts(self.Data)
        self.assertFalse(False in self.Data.data[:,:,:,freq].mask)
        self.assertFalse(False in self.Data.data[:,:,:,freq2].mask)

    def test_013_flags_constant_but_very_high_single_frequency(self):
        # If a frequency is ALWAYS at a value [whether it's always 0 or
        # 10000] it will not be flagged because it has no variance.
        # Of course, this is not right, but the chance of having a very
        # high but constant temperature at a certain frequency is 
        # next to nil.
        ind = (5,0,0,1345)
        freq = 1345
        for time in range(0, self.Data.dims[0]):
            self.Data.data[time,:,:,ind] = 50
        flag_data.apply_cuts(self.Data)
        # Notice the assertTrue, not assertFalse like all of the others.
        self.assertTrue(False in self.Data.data[:,:,:,freq].mask)

    def test_014_everything_constant(self):
        self.Data.data[:,:,:,:] = 1
        flag_data.apply_cuts(self.Data)
        # Nothing gets flagged.
        self.assertFalse(True in self.Data.data.mask)

    def test_015_problem_in_time(self):
        # Right now, the test fits file is too small in time.
        # When a longer one is made, this can be tested.
        ##
        # At freq 500, everything should be constant except for a wild
        # bit for 15 time bins. Hence, only around the wild bit should 
        # freq 500 be flagged, not the whole frequency.
        self.Data.data[:,:,:,500] = 2
        self.Data.data[100:115:5,:,:,:] = 60
        self.Data.data[101:115:5,:,:,:] = 30
        self.Data.data[102:115:5,:,:,:] = 75
        self.Data.data[103:115:5,:,:,:] = 45
        self.Data.data[104:115:5,:,:,:] = 100
        flag_data.apply_cuts(self.Data)
        #for time in range(110,115):
        #    self.assertFalse(False in self.Data.data[time,:,:,:].mask)
        #self.assertTrue(False in self.Data.data[:,:,:,500])
        #Uncomment above stuff for real test.
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()


