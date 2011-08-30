"""Unit tests for data_map and fits_map module."""

import unittest
import copy
import os

import scipy as sp

import data_map
import fits_map
import kiyopy.custom_exceptions as ce

nlat = 10
nlong = 20
nfreq = 50

class Tests(unittest.TestCase) :
    
    def setUp(self) :
        self.test_map = data_map.DataMap()
        map = (sp.arange(nlong)[:,sp.newaxis,sp.newaxis] * 
               sp.arange(nlat)[sp.newaxis,:,sp.newaxis] *
               sp.arange(nfreq)[sp.newaxis,sp.newaxis,:])
        self.test_map.set_data(map)
        for ii in range(len(fits_map.fields)) :
            self.test_map.set_field(fits_map.fields[ii], ii, format='1E')
        self.test_map.add_history('Created from scratch.', ('Detail1',
                                                            'Detail2'))
        self.test_map.verify()

    def test_verify_fails_multiD_field(self) :
        self.test_map.set_field('BLAH', sp.arange(nlat), axis_names=('lat',), 
                                format='1E')
        self.assertRaises(ce.DataError, self.test_map.verify)

    def test_circle_write_read(self) :
        map_copy = copy.deepcopy(self.test_map)
        fits_map.write(map_copy, 'temp_testmap.fits', feedback=0)
        read_map = fits_map.read('temp_testmap.fits', feedback=0)
        os.remove('temp_testmap.fits')

        self.assertTrue(sp.allclose(self.test_map.data, read_map.data))
        for field_name in fits_map.fields :
            self.assertTrue(read_map.field.has_key(field_name))
            self.assertAlmostEqual(self.test_map.field[field_name], 
                                   read_map.field[field_name])
        # Finally, check the history.
        hist = read_map.history
        self.assertTrue(hist.has_key('000: Created from scratch.'))

    def test_read_works_missing_field(self) :
        del self.test_map.field['POL']
        fits_map.write(self.test_map, 'temp_testmap.fits', feedback=0)
        read_map = fits_map.read('temp_testmap.fits', feedback=0)
        os.remove('temp_testmap.fits')


        
    def tearDown(self) :
        del self.test_map



if __name__ == '__main__' :
    unittest.main()
