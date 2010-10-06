"""Unit tests for base_single.py."""

import unittest
import os

import base_single
from core import fitsGBT, data_block

class ExampleProcessor(base_single.BaseSingle) :
    
    def nothing_to_do(self, Data) :
        Data.add_history("Did Nothing.", self.params['a_parameter'])

    actions = (nothing_to_do, )
    params_init = {'a_parameter': 'not much'}
    feedback = 0 # suppresses feedback for testing 

class AnotherExample(base_single.BaseSingle) :
    
    def something_to_do(self, Data) :
        Data.add_history("Did Something.", self.params['another_parameter'])
        Data.data[:,:,:,:] = 1.

    actions = (something_to_do, )
    params_init = {'another_parameter': 'much'}
    feedback = 0


input_parameters = {
          'input_root' : './testfile_GBTfits',
          'file_middles' : [""],
          'input_end' : ".fits",
          'output_root' : "./temp_test_",
          'output_end' : ".fits",
          'scans' : [],
          'IFs' : [],
          'a_parameter': 'nothing',
          'another_parameter': 'something'
          }

class TestBasicUsage(unittest.TestCase) :
    
    def test_runs(self) :
        ExampleProcessor(input_parameters)
        Reader = fitsGBT.Reader('temp_test_.fits', feedback=0)
        Data = Reader.read(0, 0)
        self.assertTrue(Data.history.has_key('001: Did Nothing.'))
        self.assertEqual(Data.history['001: Did Nothing.'][0], 'nothing')

    #def test_merged(self) :
    #    pass
        

    def tearDown(self) :
        os.remove('temp_test_.fits')
        os.remove('temp_test_params.ini')
    
    

if __name__ == '__main__' :
    unittest.main()
