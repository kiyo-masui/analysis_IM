"""Unit tests for freq_slices."""

import os
import unittest

import freq_slices

class Test(unittest.TestCase) :
    
    def test_runs(self):
        freq_slices.FreqSlices().execute()

    def tearDown(self):
        os.remove('testoutputparams.ini')

if __name__ == '__main__' :
    unittest.main()
