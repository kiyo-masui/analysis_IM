"""Unit tests for reflag.py"""

import unittest
import os
import glob

import measure_noise

class TestModule(unittest.TestCase) :

    def test_module(self) :
        params = {"mn_output_root" : "testout_"}
        measure_noise.Measure(params, feedback=0).execute()
        files = glob.glob('*testout*')
        self.assertTrue(len(files) > 1)

    def tearDown(self) :
        files = glob.glob('*testout*')
        for f in files :
            os.remove(f)

# Should probably do a real test here.


if __name__ == '__main__' :
    unittest.main()

