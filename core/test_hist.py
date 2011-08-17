"""Unit tests for hist.py."""

import unittest
import os
import glob

import hist

class TestReadWrite(unittest.TestCase) :

    def setUp(self) :
        self.history = hist.History()
        self.history.add('Did 1 thing', 'the details of what I did')
        self.history.add('Did something else', ('Details about the other'
            + ' thing', 'another detail', 'yet another one'))

    def test_write_runs(self) :
        self.history.write('testoutput_history.hist')

    def test_write_read(self) :
        self.history.write('testoutput_history.hist')
        history2 = hist.read('testoutput_history.hist')
        self.assertEqual(self.history, history2)

    def tearDown(self) :
        files = glob.glob('testoutput_*')
        for f in files :
            os.remove(f)



if __name__ == '__main__' :
    unittest.main()

