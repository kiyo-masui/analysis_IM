"""Unit tests for exceptions.py"""

import unittest
import custom_exceptions as ce

class TestDataError(unittest.TestCase) :
    
    def test_simple_test(self) :
        self.assertRaises(ce.DataError, self.raise_DE) 

    def raise_DE(self) :
        raise ce.DataError()

if __name__ == '__main__' :
    unittest.main()
