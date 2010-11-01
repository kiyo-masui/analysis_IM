"""Unit tests for window stitching modules."""

import unittest

import scipy as sp
import numpy.ma as ma
import matplotlib.pyplot as plt

import kiyopy.custom_exceptions as ce
import stitch_windows_crude as swc
from core import data_block, fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestStitchWindowsCrude(unittest.TestCase) :

    def setUp(self) :
        Reader = fitsGBT.Reader('testfile_GBTfits.fits', feedback=0)
        self.data_blocks = Reader.read(0,())

    def test_runs(self) :
        swc.stitch(self.data_blocks)

    def test_has_required_keys(self) :
        del self.data_blocks[0].field['CRVAL1']
        self.assertRaises(ce.DataError, swc.stitch, self.data_blocks)

    def test_all_keys_same(self) :
        self.data_blocks[0].field['SCAN'] = 3
        self.assertRaises(ce.DataError, swc.stitch, self.data_blocks)

    def tearDown(self) :
        del self.data_blocks

if __name__ == '__main__' :
    unittest.main()

