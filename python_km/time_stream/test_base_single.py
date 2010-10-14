"""Unit tests for base_single.py."""

import unittest
import os

import base_single
from core import fitsGBT, data_block

class ExampleProcessor(base_single.BaseSingle) :
    
    def nothing_to_do(self, Data) :
        Data.add_history("Did Nothing.", self.params['a_parameter'])
        return Data

    action = nothing_to_do
    params_init = {'a_parameter': 'not much'}
    feedback = 0 # suppresses feedback for testing
    prefix = 'nt_'


input_parameters = {
          'nt_input_root' : './testfile_GBTfits',
          'nt_file_middles' : (""),
          'nt_input_end' : ".fits",
          'nt_output_root' : "./temp_test_",
          'nt_output_end' : ".fits",
          'nt_scans' : (),
          'nt_IFs' : (),
          'nt_a_parameter' : 'nothing',
          'a_parameter' : 'never read'
          }

class TestBasicUsage(unittest.TestCase) :
    
    def test_runs(self) :
        ExampleProcessor(input_parameters).execute()
        Reader = fitsGBT.Reader('temp_test_.fits', feedback=0)
        Data = Reader.read(0, 0)
        self.assertTrue(Data.history.has_key('001: Did Nothing.'))
        self.assertEqual(Data.history['001: Did Nothing.'][0], 'nothing')

    def tearDown(self) :
        os.remove('temp_test_.fits')
        #os.remove('temp_test_params.ini')
    
    

if __name__ == '__main__' :
    unittest.main()
