"""Unit tests for pipeline manager.

These aren't really tests, this just makes sure the whole thing runs and
provides examples of how it would work.
"""

import unittest
import os

import manager

# Companion test input file.
input_test_file_name = 'pipeline/test.pipe'

class TestManager(unittest.TestCase) :
    
    def test_executes(self) :
        manager.execute(input_test_file_name)

    def tearDown(self) :
        os.remove('testoutput_params.ini')
        os.remove('testoutput_testfile_GBTfits.testhanning.fits')
        os.remove('testoutput_testfile_GBTfits.testflag.fits')
        os.remove('testoutput_testfile_GBTfits.testrebin.fits')
        os.remove('testoutput_map.testmap.fits')
        os.remove('testoutput_noise.testmap.fits')

if __name__ == '__main__' :
    unittest.main()
