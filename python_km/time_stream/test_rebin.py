"""Unit tests for hanning a rebin_freq modules."""

import unittest
import copy

import scipy as sp
import numpy.ma as ma

import sys

import rebin_freq
import core.data_block
import core.fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestFunctions(unittest.TestCase) :
    """Since these operations actually changes the data, these are only sanity
    tests, far from thorough."""
    
    def setUp(self) :
        Reader = core.fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,1)
        self.Data.verify()
        self.Data_copy = copy.deepcopy(self.Data)

    def test_rebin_runs(self) :
        rebin_freq.rebin(self.Data, 1.0)
        self.Data.verify()
        self.assertAlmostEqual(-1.0, self.Data.field['CDELT1']/1.0e6)

    def test_rebin_get_freqs_right(self) :
        """Make sure that the frequencies calculated by the rebinner are
        consistant with the ones calculated by calc_freq()."""
        rebin_freq.rebin(self.Data, 2.0)
        old_freq = sp.array(self.Data.freq)
        self.Data.calc_freq()
        self.assertTrue(sp.allclose(old_freq, self.Data.freq))
        self.Data.verify()
        self.assertEqual(self.Data.dims[-1], len(self.Data.freq))

    def test_end_freqs(self) :
        self.Data.calc_freq()
        old_top = self.Data.freq[-1]
        old_bot = self.Data.freq[0]
        rebin_freq.rebin(self.Data, 0.9)
        self.assertAlmostEqual(old_top, self.Data.freq[-1], places=-7)
        self.assertAlmostEqual(old_bot, self.Data.freq[0], places=-7)

    def test_rebins_data_right(self) :
        """set all data to 0 except 1 entry, figure out the f, rebin and
        check that only the closest new f is non-zero."""
        self.Data.data[:,:,:,:] = 0
        self.Data.data[5,2,0,324] = 1
        self.Data.calc_freq()
        oldf = self.Data.freq[324]
        rebin_freq.rebin(self.Data, 2.0, mean=True)
        new_ind = (abs(self.Data.freq - oldf)).argmin()
        self.assertNotAlmostEqual(0, self.Data.data[5,2,0,new_ind])
        self.Data.data[5,2,0,new_ind] = 0.0
        self.assertTrue(sp.allclose(self.Data.data, 0.0))

    def test_by_number_axis(self) :
        # For now rebin to a power of 2 for no residual
        wbins = 2**4
        # Get what we'd expect the new frequency axis to be
        self.Data.calc_freq()
        old_freq = self.Data.freq
        new_freq = old_freq[wbins//2::wbins]
        delta =  self.Data.field['CDELT1']
        # Rebin.
        rebin_freq.rebin(self.Data, wbins, by_nbins=True)
        self.Data.verify()
        # Check that the spaceing is right.
        self.assertAlmostEqual(wbins*delta, self.Data.field['CDELT1'])
        # Check that the axis is what we would think.
        self.Data.calc_freq()
        self.assertTrue(sp.allclose(new_freq, self.Data.freq, atol=abs(delta)))

    def test_by_number_data(self) :
        wbins = 9
        # Set data = f.
        self.Data.calc_freq()
        self.Data.data[...] = self.Data.freq
        delta =  self.Data.field['CDELT1']        
        # Rebin.
        rebin_freq.rebin(self.Data, wbins, by_nbins=True)
        self.Data.verify()
        # Except for the last bin, Data should still be freq.
        self.Data.calc_freq()
        self.assertTrue(sp.allclose(self.Data.data[:,:,:,:-1], 
                                    self.Data.freq[:-1], abs(delta)))
        # In last bin, same, but tolerance is worse.
        self.assertTrue(sp.allclose(self.Data.data[:,:,:,-1], 
            self.Data.freq[-1], abs(self.Data.field['CDELT1'])))




    def tearDown(self) :
        del self.Data
        del self.Data_copy


if __name__ == '__main__' :
    unittest.main()

