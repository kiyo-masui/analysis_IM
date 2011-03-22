"""Unit tests for pipeline manager.

These aren't really tests, this just makes sure the whole thing runs and
provides examples of how it would work.
"""

import unittest
import os
import glob

import manager

# Companion test input file.
input_test_file_name = 'pipeline/test.pipe'

class TestManager(unittest.TestCase) :
    
    def test_executes(self) :
        manager.execute(input_test_file_name, 0)

    def tearDown(self) :
        out_files = glob.glob("testoutput*")
        for f in out_files:
            os.remove(f)

if __name__ == '__main__' :
    unittest.main()
