"""Unit tests for map maker."""

import unittest
import scipy as sp
import numpy.ma as ma

import map_maker_simple as mm
import tools

class TestAddingData(unittest.TestCase) :
    
    def setUp(self) :
        self.shape = (5, 5, 20)
        self.ntime = 15
        self.map = ma.zeros(self.shape, dtype=float)
        self.counts = sp.zeros(self.shape, dtype=int)
        self.data = ma.ones((self.ntime, self.shape[-1]))
        self.centre = (90.0,0.0,100.0)
        self.spacing = 1.0
        self.ra = 0.2*(sp.arange(self.ntime)-self.ntime/2.0) + self.centre[0]
        self.dec = 0.2*(sp.arange(self.ntime)-self.ntime/2.0) + self.centre[1]
        self.ra_inds = tools.calc_inds(self.ra, self.centre[0], self.shape[0],
                                    self.spacing)
        self.dec_inds = tools.calc_inds(self.dec, self.centre[1], self.shape[1],
                                    self.spacing)
        
    def test_adds_data(self) :
        mm.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.counts)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1])
        self.assertAlmostEqual(ma.sum(self.counts), self.ntime*self.shape[-1])

    def test_mask_data(self) :
        self.data[5,8] = ma.masked
        mm.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.counts)
        self.assertAlmostEqual(ma.sum(self.map), self.ntime*self.shape[-1]-1)
        self.assertAlmostEqual(ma.sum(self.counts), self.ntime*self.shape[-1]-1)

    def test_off_map(self) :
        self.ra_inds[3] = 10
        self.dec_inds[6] = -1
        mm.add_data_2_map(self.data, self.ra_inds, self.dec_inds, 
                          self.map, self.counts)
        self.assertAlmostEqual(ma.sum(self.map),
                               self.ntime*self.shape[-1]-2*self.shape[2])
        self.assertAlmostEqual(ma.sum(self.counts), 
                               self.ntime*self.shape[-1]-2*self.shape[2])

    # test a pointing we know is right

    def test_checks_shapes(self) :
        self.assertRaises(ValueError, mm.add_data_2_map, self.data, self.ra,
                          sp.arange(self.ntime+1), self.map, self.counts)
        self.assertRaises(ValueError, mm.add_data_2_map, self.data, self.ra,
                          self.dec, self.map[:,:,0:-1], self.counts)

    def tearDown(self) :
        del self.data
        del self.map
        del self.counts
        


if __name__ == '__main__' :
    unittest.main()
