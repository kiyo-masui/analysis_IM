"""Unit tests for utils."""

import unittest
import os

import utils


class TestMKDIR(unittest.TestCase) :
    """Unit tests for mkdir_p."""
    
    def test_works(self) :
        utils.mkdir_p('tmp/tmp2')
        self.assertTrue(os.path.exists('tmp/tmp2'))
        # This should not raise an error.
        utils.mkdir_p('tmp/tmp2')
    
    def test_parents(self) :
        utils.mkparents('tmp/tmp2/stuff')
        self.assertTrue(os.path.exists('tmp/tmp2'))
        self.assertTrue(not os.path.exists('tmp/tmp2/stuff'))
        # Test that it does nothing for a simple file name but doesn't fail.
        utils.mkparents('stuff')
        self.assertTrue(not os.path.exists('stuff'))

    def tearDown(self) :
        if os.path.exists ('tmp/tmp2') :
            os.rmdir('tmp/tmp2')
            os.rmdir('tmp')

class TestFileNameAbbr(unittest.TestCase) :
    
    def test_works(self) :
        fname ='/a_place/in_a_very/deep/directory/with_a_very/long/name.dat'
        abbr = utils.abbreviate_file_path(fname)
        self.assertEqual(abbr, 'long/name.dat')
        abbr = utils.abbreviate_file_path(abbr)
        self.assertEqual(abbr, 'long/name.dat')
        fname = 'name'
        abbr = utils.abbreviate_file_path(fname)
        self.assertEqual(abbr, fname)

if __name__ == '__main__' :
    unittest.main()

