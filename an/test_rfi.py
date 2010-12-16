"""Check that the SWIG build rfi module works."""

import unittest

import scipy as sp

import rfi

class TestWorks(unittest.TestCase) :
    
    def test_get_fit(self) :
        f = sp.arange(2048, dtype=float)
        a = sp.arange(2048, dtype=float)
        b = sp.empty(2048, dtype=float)
        rfi.get_fit(a, f, b)
        self.assertTrue(sp.allclose(a, b))

if __name__ == '__main__' :
    unittest.main()

