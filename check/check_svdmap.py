True#! /usr/bin/env python 

import cPickle
import numpy as np
import scipy as sp
import scipy.ndimage.interpolation as ndimage_inter
import core.algebra as algebra
import matplotlib.pyplot as plt
import pyfits
import os
import sys
import copy
import gc

from map import beam
from parkes import cal_map

class plotmap(object):
    def __init__(self, mapdict, imap=None, qmap=None, umap=None, 
                 freq_cut=[], degrade_factor=0, plot_size=(10,8)):
        if imap is None:
            self.imap = algebra.make_vect(algebra.load(mapdict['imap']))
            if degrade_factor != 0:
                self.imap = self.degrade_resolution(self.imap, degrade_factor)
                self.degrade_factor=degrade_factor
            else:
                self.degrade_factor = 1.

            self.nmap = algebra.load(mapdict['nmap'])
            #self.imap *= self.nmap
        else:
            self.imap = imap


        self.name = mapdict['name']

        print self.name

        self.freq_list = tuple([ind for ind in range(self.imap.shape[0]) 
                                    if  ind not in freq_cut])
        #self.imap = self.imap[freq_list, ...]

        if 'mode' in mapdict.keys():
            self.mode = algebra.make_vect(algebra.load(mapdict['mode']))
        if 'vect' in mapdict.keys():
            self.svdmodpkl = cPickle.load(open(mapdict['vect']))
        #ra = self.imap.get_axis('ra')
        #print ra
        self.imap_secs = []
        if 'imap_sec' in mapdict.keys():
            for i in range(len(mapdict['imap_sec'])):
                imap_sec = algebra.make_vect(algebra.load(mapdict['imap_sec'][i]))
                if degrade_factor != 0:
                    imap_sec = self.degrade_resolution(imap_sec, degrade_factor)
                self.imap_secs.append(imap_sec)

        self.re_scale = None
        self.re_zerop = None
        self.plot_size = plot_size
        self.hipass_map_path = os.getenv('HIPASS_PATH')

    def test_fitting(self):
        calabrator = cal_map.CalibrateMap(self.imap, self.freq_list, self.nmap, 
                                          degrade_function = self.degrade_function())
        calabrator.cal_by_fitting(self.hipass_map_path)

    def degrade_function(self, degrade_factor=None):

        if degrade_factor == None:
            degrade_factor = self.degrade_factor

        freq_data = sp.array([1250, 1275, 1300, 1325, 1350, 1430], dtype=float)
        beam_data = sp.array([14.4, 14.4, 14.4, 14.4, 14.4, 14.4])/60. 
        beam_data = beam_data*1420/freq_data
        freq_data *= 1.0e6
        beam_diff = sp.sqrt(max(degrade_factor*beam_data)**2-(beam_data)**2)
        common_resolution = beam.GaussianBeam(beam_diff, freq_data)

        return common_resolution



if __name__=='__main__':

    mapdict = {}
    mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_clean_map_I_1315.npy'
    mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
    mapdict['name'] = 'test_allbeams_27n30_10by7_noise_inv_I'
    mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_A_noise_weight_I_1315.npy',
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_B_noise_weight_I_1315.npy', 
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_C_noise_weight_I_1315.npy', 
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_D_noise_weight_I_1315.npy', 
                          ]
    #degrade_factor = 1.1
    degrade_factor = 0
    if degrade_factor != 0:
        mapdict['name'] += '_1pt1_cov'

    a = plotmap(mapdict, freq_cut = [0,1,2,3,4,5,59,60,61,62,63], 
                degrade_factor=degrade_factor, plot_size=(15, 6))
    a.test_fitting()
    #a.mapping_coord(plot=True)
    #a.hipass_cat()
    #a.plot_map(with_nvss=True)
    #a.plot_map(with_nvss=True, diff=True)
    #a.check_nvss(flux_limit=[0.5, 100], rescale=False)
    #a.check_nvss(flux_limit=[0.5, 100])
    #a.plot_map(with_nvss=True)
    #a.plot_map(with_nvss=True, diff=True)

