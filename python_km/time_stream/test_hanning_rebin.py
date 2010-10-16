"""Unit tests for preprocess module."""

import unittest
import copy

import scipy as sp
import numpy.ma as ma

import preprocess
import core.data_block
import core.fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestFunction(unittest.TestCase) :
    """Since these operations actually changes the data, these are only sanity
    tests, far from thorough."""
    
    def setUp(self) :
        Reader = core.fitsGBT.Reader(test_file)
        self.Data = Reader.read(1,1)
        self.Data.verify()
        self.Data_copy = copy.deepcopy(self.Data)

    def sanity_test(self) :
        # For lack of anything better to test:
        self.Data.verify()
        # Make sure we actually did something.
        self.assertTrue(not ma.allclose(self.Data.data, self.Data_copy.data))
        # Make sure we didn't change other fields, like LST.
        self.assertTrue(sp.allclose(self.Data.field['LST'],
                        self.Data_copy.field['LST']))

    def test_hanning_data_changed(self) :
        """Copy the data, see that we did something."""
        preprocess.hanning(self.Data)
        self.sanity_test()

    # Could test that end pointes end up masked.  Could test that points
    # adjacent to masked data end up masked.
        
    def tearDown(self) :
        del self.Data
        del self.Data_copy


if __name__ == '__main__' :
    unittest.main()

