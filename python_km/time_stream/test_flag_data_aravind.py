"""Unit tests for flag_data.py"""

import unittest

import scipy as sp
import numpy.ma as ma
import numpy.random as rand
import matplotlib.pyplot as plt

import kiyopy.custom_exceptions as ce
import time_stream.flag_data_aravind as flag_data
from core import data_block, fitsGBT

test_file = 'testfile_GBTfits.fits'

class TestFlagData(unittest.TestCase) :

    def setUp(self) :
        Reader = fitsGBT.Reader(test_file, feedback=0)
        self.Data = Reader.read(1,0)

    def test_pol_cut_cal_off(self) :
        ind = (6,1,1,676)
        self.Data.data[ind] = 10
        flag_data.apply_cuts(self.Data, -1, 3)
        self.assertTrue(self.Data.data[ind] is ma.masked)
    
    def test_pol_cal_off_controled(self) :
        self.Data.data[0::2,[0,3],:,:] = rand.normal(0,1,2048)
        self.Data.data[1::2,[0,3],:,:] = rand.normal(0,1.1,2048)
        self.Data.data[:,[1,2],:,:] = rand.normal(0,0.01,2048)
        self.Data.data[1::2,[1,2],:,:] = rand.normal(0,0.1,2048)
        ind = (5,1,1,1345)
        ind2 = (2,2,0,425)
        self.Data.data[ind] = 25.
        self.Data.data[ind2] = 0.5
        flag_data.apply_cuts(self.Data, -1, 3.0, 3, True, 0, 0)
        self.assertTrue(self.Data.data[ind] is ma.masked)
        self.assertTrue(self.Data.data[ind2] is ma.masked)
        # 8 = 2 bad data * 4 polarizations
        self.assertTrue(ma.count_masked(self.Data.data) == 8)


    def some_plots(self) :
        """Code I used to figure out my algorithms."""
        data = self.Data.data
        cal_ind = 1
        cross = (data[:,[1,2],cal_ind,:]**2
                 /data[:,[0],cal_ind,:]/data[:,[3],cal_ind,:])
        cross = ma.sum(cross, 1)
        plt.figure(0)
        #plt.semilogy(ma.mean(cross,0))
        plt.semilogy(cross[9,:])
        plt.figure(1)
        plt.semilogy(ma.std(cross,0))
        #plt.plot(cross[1,:])
        #print sp.sqrt(ma.median(cross**2))

        xx_ind = 0
        yy_ind = 3
        on_ind = 0
        off_ind = 1
        xy_inds = [1,2]

        #cal = (ma.median(data[:,[xx_ind,yy_ind],on_ind,:], 0) - 
        #       ma.median(data[:,[xx_ind,yy_ind],off_ind,:], 0))**2
        #cal = ma.sum(cal, 0)
        #cross = ma.sum(data[:,xy_inds,on_ind,:]**2, 1)
        #cross = cross - cal
        
        cal = ma.median(data[:,xy_inds,on_ind,:]-data[:,xy_inds,off_ind,:], 0) 
        cross = ma.sum((data[:,xy_inds,on_ind,:] - cal)**2, 1)
        
        
        cross = cross/data[:,xx_ind,on_ind,:]/data[:,yy_ind,on_ind,:]
        plt.figure(2)
        #plt.semilogy(ma.mean(cross,0))
        plt.semilogy(cross[9,:])
        plt.figure(3)
        plt.semilogy(ma.std(cross,0))


        plt.show()



    def tearDown(self) :
        del self.Data

if __name__ == '__main__' :
    unittest.main()

