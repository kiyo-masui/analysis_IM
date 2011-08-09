"""Unit tests for data_block module."""

import unittest

import scipy as sp
import numpy.ma as ma

import data_block as db
import fitsGBT
import kiyopy.custom_exceptions as ce

ntimes = 5
npol = 4
ncal = 2
nfreq = 10
dims = (ntimes, npol, ncal, nfreq)

test_fits_file_name = 'testdata/testfile_GBTfits.fits'

def standard_data_test(utest_case) :
    """Encapsulates a few tests on TestDB.data.

    These test should pass any time we expect the a data block to have valid
    data."""

    for ii in range(4) :
        utest_case.assertEqual(dims[ii], utest_case.TestDB.dims[ii])
    utest_case.assertTrue(utest_case.TestDB.data_set)
    dir_data = dir(utest_case.TestDB.data)
    utest_case.assertTrue(dir_data.count('mask'))


class TestDataSetup(unittest.TestCase) :
    
    def setUp(self) :
        self.testdata = sp.reshape(sp.arange(ntimes*npol*ncal*nfreq,
                                   dtype=float), (ntimes,npol,ncal,nfreq))
        
    def test_init_sets_data(self) :
        self.TestDB = db.DataBlock(self.testdata)
        standard_data_test(self)

    def test_init_no_data(self) :
        self.TestDB = db.DataBlock()
        self.assertTrue(not self.TestDB.data_set)
        self.TestDB.set_data(self.testdata)
        standard_data_test(self)

    def tearDown(self) :
        del self.TestDB
        del self.testdata

class TestFields(unittest.TestCase) :
    
    def setUp(self) :
        # make valid data
        self.testdata = sp.reshape(sp.arange(ntimes*npol*ncal*nfreq,
                                   dtype=float), (ntimes,npol,ncal,nfreq))
        self.TestDB = db.DataBlock(self.testdata)
        self.LST_test = sp.arange(ntimes)
        self.TestDB.set_field('LST', self.LST_test, ('time',), '1E')
        # Make sure that 'verify' passes on good data:
        self.TestDB.verify()

    def test_set_get_field(self) :
        stored_LST = sp.array(self.TestDB.field['LST'])
        # The obviouse tests.
        for ii in range(ntimes) :
            self.assertAlmostEqual(self.TestDB.field['LST'][ii], 
                                   stored_LST[ii])
        self.assertEqual(self.TestDB.field_axes['LST'][0], 'time')
        self.assertEqual(self.TestDB.field_formats['LST'], '1E')
        # Test that it's checking for a valid axis.
        self.assertRaises(ValueError, self.TestDB.set_field, 'abc', 5, 
                          ('not_a_axis',), '1E')
        # Test that we've copied dereferenced data into the DataBlock
        self.LST_test[0] = -100.
        self.assertEqual(self.TestDB.field['LST'][0], stored_LST[0])

    def test_verify_keys_match1(self) :
        self.TestDB.field['stuff'] = range(5)
        self.assertRaises(ce.DataError, self.TestDB.verify)
        
    def test_verify_keys_match2(self) :
        """Oppossite case of the last test."""
        self.TestDB.field_axes['stuff'] = ('time',)
        self.assertRaises(ce.DataError, self.TestDB.verify)

    def test_verify_valid_axes(self) :
        self.TestDB.field_axes['LST'] = ('not_an_axis',)
        self.assertRaises(ValueError, self.TestDB.verify)

    def test_multiD_field(self) :
        self.TestDB.set_field('afield', sp.ones((npol,ncal)), ('pol', 'cal'), '1I')
        self.TestDB.verify()
        self.assertTrue(sp.allclose(self.TestDB.field['afield'],
                                    sp.ones((npol,ncal))))
        self.assertEqual(self.TestDB.field_axes['afield'], ('pol', 'cal'))
        # fitsGBT depends on these tests.  If you generalize this, add tests to
        # fitsGBT.
        self.assertRaises(ValueError, self.TestDB.set_field, 'afield',
                          sp.ones((npol,npol)), ('pol', 'pol'), '1I')
        self.assertRaises(ValueError, self.TestDB.set_field, 'afield',
                          sp.ones((ncal,npol)), ('cal', 'pol'), '1I')

    def test_zeroD_fields(self) :
        self.TestDB.set_field('SCAN', 113, format='1I')
        self.TestDB.verify()
        self.assertEqual(self.TestDB.field['SCAN'], 113)
        self.assertEqual(len(self.TestDB.field_axes['SCAN']), 0)

    def test_verify_field_shape(self) :
        self.TestDB.set_field('CRVAL4', [-5,-6,-7,-8,-9], 'pol', format='1I')
        self.assertRaises(ce.DataError, self.TestDB.verify)
        
    def test_verify_field_format(self) :
        """The class should really do something more sophisticated and
        acctually check that the fits format string is a valid one for the
        data.  Someone should code this."""
        self.TestDB.set_field('CRVAL4', [-5,-6,-7,-8], 'pol', format=10)
        self.assertRaises(ce.DataError, self.TestDB.verify)
        
    def tearDown(self) :
        del self.TestDB
        del self.testdata
        del self.LST_test

class TestHistory(unittest.TestCase) :
    
    def setUp(self) :
        self.TestDB = db.DataBlock()
        self.hist_str = 'For example, saved to file.'
        self.hist_detail = ('the file name', 'the date perhaps')

    def test_add_history(self) :
        # The basics:
        self.TestDB.add_history(self.hist_str, self.hist_detail)
        self.assertTrue(self.TestDB.history.has_key('000: '+self.hist_str))
        self.assertTrue(self.TestDB.history['000: '+self.hist_str][0] == 
                        self.hist_detail[0])
        # Limit the length
        self.assertRaises(ValueError, self.TestDB.add_history, str(10**70), 
                          self.hist_detail)
        self.assertRaises(ValueError, self.TestDB.add_history, self.hist_str, 
                          (str(10**70),))
        self.assertRaises(TypeError, self.TestDB.add_history, 5, 
                          self.hist_detail)
        self.assertRaises(TypeError, self.TestDB.add_history, self.hist_str, 
                          (1,2))

    def test_add_hist_no_detail(self) :
        self.TestDB.add_history(self.hist_str)
        self.assertTrue(self.TestDB.history.has_key('000: '+self.hist_str))
        self.assertEqual(len(self.TestDB.history['000: '+self.hist_str]), 0)

    def test_add_string_detail(self) :
        self.TestDB.add_history(self.hist_str, self.hist_detail[0])
        self.assertTrue(self.TestDB.history['000: '+self.hist_str][0] == 
                        self.hist_detail[0])
    
    # No unit tests for print history.  It is a function for user reading and
    # is thus not prone to propagable errors.

    def test_merge_histories(self) :
        # Basic tests
        self.TestDB.add_history(self.hist_str, self.hist_detail)
        SecondDB = db.DataBlock()
        second_details = ('second file name', )
        SecondDB.add_history(self.hist_str, second_details)
        merged = db.merge_histories(self.TestDB, SecondDB)
        self.assertEqual(len(merged.keys()), 1)
        self.assertTrue(second_details[0] in merged['000: '+self.hist_str])
        self.assertTrue(self.hist_detail[0] in merged['000: '+self.hist_str])
        self.assertEqual(len(merged['000: '+self.hist_str]), 3)
        # Enforce matching.
        ThirdDB = db.DataBlock()
        ThirdDB.add_history(self.hist_str, self.hist_detail)
        ThirdDB.add_history('Read from file.', self.hist_detail)
        self.assertRaises(ce.DataError, db.merge_histories, SecondDB, ThirdDB)

    def test_merge_multiple_histories(self) :
        entry1 = 'Read from file.'
        entry2 = 'Saved to file.'
        DBs = ()
        n = 10
        for ii in range(n) :
            tempDB = db.DataBlock()
            tempDB.add_history(entry1, 'filename: ' + str(ii))
            tempDB.add_history(entry2, 'saved filename not iterated')
            DBs = DBs + (tempDB, )
        merged = db.merge_histories(*DBs)
        self.assertEqual(len(merged['000: '+entry1]), n)
        self.assertEqual(len(merged['001: '+entry2]), 1)
    
    def tearDown(self) :
        del self.TestDB

class TestGBTmethods(unittest.TestCase) :
    """Unit tests for methods that should only work if the Data Block is from a
    GBT fits file.  Assumes specific fields have been set."""
    
    # Parameters of the test file data. Accessed by self.ntime (don't use the
    # globals)
    ntime = 10
    nfreq = 2048

    def setUp(self) :
        self.Reader = fitsGBT.Reader(test_fits_file_name, 0)
        self.Data = self.Reader.read(0, 0)

    def test_calculates_pointing(self) :
        self.Data.calc_pointing()
        self.assertTrue(hasattr(self.Data, 'ra'))
        self.assertTrue(hasattr(self.Data, 'dec'))
        self.assertEqual(len(self.Data.ra), self.ntime)
        self.assertEqual(len(self.Data.dec), self.ntime)

    def test_calculates_frequencies(self) :
        self.Data.calc_freq()
        self.assertTrue(hasattr(self.Data, 'freq'))
        self.assertEqual(len(self.Data.freq), self.nfreq)
        self.assertAlmostEqual(self.Data.field['BANDWID'], 
                               sp.amax(self.Data.freq) -
                               sp.amin(self.Data.freq), -5)

    def test_calculates_time(self) :
        self.Data.calc_time()
        self.assertTrue(hasattr(self.Data, 'time'))
        t_copy = sp.copy(self.Data.time)
        t_copy.sort()
        self.assertTrue(sp.allclose(t_copy, self.Data.time))


if __name__ == '__main__' :
    unittest.main()
