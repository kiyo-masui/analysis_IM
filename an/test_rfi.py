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

    def test_get_fit_size(self) :
        f = sp.arange(32*20, dtype=float)
        a = sp.arange(32*20, dtype=float)
        b = sp.empty(32*20, dtype=float)
        status = rfi.get_fit(a, f, b)
        self.assertTrue(sp.allclose(a, b))

    def test_get_fit_fails_32(self) :
        f = sp.arange(31*20, dtype=float)
        a = sp.arange(31*20, dtype=float)
        b = sp.empty(31*20, dtype=float)
        status = rfi.get_fit(a, f, b)
        self.assertTrue(status)
    
    def test_fails_diff_sizes(self) :
        f = sp.arange(5432, dtype=float)
        a = sp.arange(5432, dtype=float)
        b = sp.empty(5153, dtype=float)
        status  = rfi.get_fit(a, f, b)
        self.assertTrue(status)

if __name__ == '__main__' :
    unittest.main()

