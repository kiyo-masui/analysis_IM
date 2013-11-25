#! /usr/bin/env python 
import numpy as np
import scipy as sp
import scipy.ndimage.interpolation as ndimage_inter
import core.algebra as algebra
import pyfits
import os
import sys
import copy

def nvss_cat(nvss_catalog_path):
    nvss_cat = np.genfromtxt(nvss_catalog_path, delimiter='|', 
            usecols = tuple(range(1,8)), 
            dtype="S19, S12, S12, f4, f4, f4, f4", 
            names=['nvss','ra','dec','ra_err','dec_err', 'flux','flux_err'])
    
    ra_raw  = np.array([x.split(' ') for x in nvss_cat['ra']])
    dec_raw = np.array([x.split(' ') for x in nvss_cat['dec']])
    
    ra_symbol = np.array([-1. if x < 0 else 1. for x in ra_raw[:,0]])
    ra_degree = ra_symbol\
              *(np.abs(ra_raw[:,0].astype('float32'))\
              + ra_raw[:,1].astype('float32')/60.\
              + ra_raw[:,2].astype('float32')/60./60.)
    ra_degree *= 15.
    
    ra_degree[ra_degree>180] = ra_degree[ra_degree>180] - 360.

    dec_symbol= np.array([-1. if x < 0 else 1. for x in dec_raw[:,0].astype('float32')])
    dec_degree= dec_symbol\
              *(np.fabs(dec_raw[:,0].astype('float32'))\
              + dec_raw[:,1].astype('float32')/60.\
              + dec_raw[:,2].astype('float32')/60./60.)

    nvss_cat['ra'] = ra_degree
    nvss_cat['dec']= dec_degree

    return nvss_cat

def hipass_cat(hipass_map_path):

    hdulist = pyfits.open(hipass_map_path)

    hd = hdulist[0].header
    hipass_data = hdulist[0].data
    hipass_ra   = (np.arange(hipass_data.shape[1])-hd['CRPIX1'] + 1)*hd['CDELT1']
    hipass_ra  += hd['CRVAL1']
    #hipass_ra  += hd['CDELT1']
    hipass_ra[hipass_ra>180] = hipass_ra[hipass_ra>180] - 360.
    hipass_dec  = (np.arange(hipass_data.shape[0])-hd['CRPIX2'] + 1)*hd['CDELT2']
    hipass_dec += hd['CRVAL2']
    #hipass_dec += hd['CDELT2']

    hipass_delta_ra = hipass_ra[1] - hipass_ra[0]
    hipass_ra_bin_edges = np.append(hipass_ra - 0.5*hipass_delta_ra,
                                    hipass_ra[-1] + hipass_delta_ra)
    hipass_delta_dec = hipass_dec[1] - hipass_dec[0]
    hipass_dec_bin_edges = np.append(hipass_dec - 0.5*hipass_delta_dec,
                                     hipass_dec[-1] + hipass_delta_dec)

    return hipass_data, hipass_ra_bin_edges, hipass_dec_bin_edges, hipass_ra, hipass_dec

def mapping_coord(imap, hipass_map_path, degrade_function=None):

    hipass_data, ra_bin_edges, dec_bin_edges, ra, dec = hipass_cat(hipass_map_path)

    ra = imap.get_axis('ra')
    dec= imap.get_axis('dec')


    ra_map_index  = np.digitize(ra,  ra_bin_edges) - 1
    dec_map_index = np.digitize(dec, dec_bin_edges) - 1

    in_the_map = np.logical_and(ra_map_index>=0, 
                 ra_map_index<hipass_data.shape[1])
    ra_map_index = ra_map_index[in_the_map]

    in_the_map = np.logical_and(dec_map_index>=0,
                 dec_map_index<hipass_data.shape[0])
    dec_map_index = dec_map_index[in_the_map]

    DEC, RA = np.meshgrid(dec_map_index, ra_map_index)


    hipass_data_new = ndimage_inter.map_coordinates(hipass_data, 
                                                    [DEC, RA], 
                                                    order=0, 
                                                    mode='nearest')
    hipass_data_new = hipass_data_new[None, ...]
    hipass_data_new = algebra.make_vect(hipass_data_new, 
                                        axis_names=imap.info['axes'])
    hipass_data_new.info = copy.deepcopy(imap.info)
    hipass_data_new.info['freq_centre'] = 1420.e6

    if degrade_function != None:
        hipass_data_new = degrade_function.apply(hipass_data_new)
    hipass_data_new = hipass_data_new[0,...]

    return hipass_data_new

class CalibrateMap(object):
    def __init__(self, imap, freq_list, nmap=None, weight=None, degrade_map=True, 
            degrade_function=None):

        self.imap = imap
        self.nmap = nmap
        self.weight = weight
        self.freq_list = freq_list
        self.re_scale = None
        self.re_zerop = None

        # degrade_function should be a object of beam.GaussianBeam()
        self.degrade_function = degrade_function

        # convolve the map into common beam
        if not degrade_map:
            self.imap = self.degrade_resolution(self.imap)

    def cal_by_fitting(self, hipass_map_path, point_sources=[], mask_pixels=None):
        hipass = self.mapping_coord(hipass_map_path)
        parkes = copy.deepcopy(self.imap[self.freq_list, ...])
        parkes = np.mean(parkes, axis=0)
        if self.weight != None:
            weight = copy.deepcopy(self.weight[self.freq_list, ...])
            weight[weight==0] = np.inf
            weight = 1. / weight
            weight = np.mean(weight, axis=0)
            weight[weight==0] = np.inf
            weight = 1. / weight
        else:
            weight = np.ones(parkes.shape)

        if len(point_sources) == 0:
            hipass_vector = hipass.flatten()[:, None]
            parkes_vector = parkes.flatten()[:, None]
            weight_vector = weight.flatten()[:, None]
        elif mask_pixels == None:
            hipass_vector = hipass.flatten().take(point_sources)[:, None]
            parkes_vector = parkes.flatten().take(point_sources)[:, None]
            weight_vector = weight.flatten().take(point_sources)[:, None]
        else:
            mask = self.get_mask(point_sources, parkes.shape, mask_pixels)
            parkes[mask[0], mask[1]] = 0.
            weight[mask[0], mask[1]] = 0.
            hipass[mask[0], mask[1]] = 0.
            hipass_vector = hipass.flatten()[:, None]
            parkes_vector = parkes.flatten()[:, None]
            weight_vector = weight.flatten()[:, None]

        W = weight_vector * np.eye(weight_vector.shape[0])
        A = np.append(parkes_vector, np.ones(parkes_vector.shape), axis=1)
        P = np.linalg.inv(np.dot(np.dot(A.T, W), A))
        Q = np.dot(A.T, W)

        re_scale, re_zerop = np.dot(np.dot(P, Q), hipass_vector).flatten()

        return re_scale, re_zerop

    def get_mask(self, points, mapshape, pixels=2):

        point_ra, point_dec = np.unravel_index(points, mapshape)

        mask_length = 2 * pixels + 1
        d = np.arange(mask_length).astype(int) - 0.5 * (mask_length - 1)

        point_ra  = point_ra[:,None] + d[None, :]
        point_dec = point_dec[:,None] + d[None, :]

        mask = np.zeros([points.shape[0], mask_length, mask_length, 2])
        mask[...,0] = point_ra[:,:,None]
        mask[...,1] = point_dec[:,None,:]

        mask = mask.reshape(-1, 2)
        mask = mask.astype(int)
        mask[mask[:,0] > mapshape[0] - 1, 0] = mapshape[0] - 1
        mask[mask[:,0] < 0, 0] = 0
        mask[mask[:,1] > mapshape[1] - 1, 1] = mapshape[1] - 1
        mask[mask[:,1] < 0, 1] = 0

        return mask.T


    def cal_by_point_sources(self, nvss_catalog_path=None, hipass_map_path=None, 
                             flux_limit=None, ra_limit=None, dec_limit=None, 
                             by_fitting=True, mask_point=True, output=True):
        '''
            calibrate the map using the hipass flux of nvss sources
        '''
        if nvss_catalog_path == None:
            nvss_catalog_path = os.getenv('NVSS_PATH')
        if hipass_map_path == None:
            hipass_map_path = os.getenv('HIPASS_PATH')

        if output:
            print 'Calibrate the map using hipass flux of nvss sources'
        else:
            print '[Cal ',
            if by_fitting:
                print 'fitting ',
                if mask_point:
                    print 'mask point',
            else:
                print 'rescaling ',
            print ']',

        # get the ra and dec bin of our map
        ra_bin = self.imap.get_axis('ra')
        dec_bin= self.imap.get_axis('dec')

        # get the ra and dec bin edges of our map
        ra_bin_delta = ra_bin[1] - ra_bin[0]
        ra_bin_edges = np.append(ra_bin - 0.5*ra_bin_delta, 
                                 ra_bin[-1] + 0.5*ra_bin_delta)
        dec_bin_delta = dec_bin[1] - dec_bin[0]
        dec_bin_edges = np.append(dec_bin - 0.5*dec_bin_delta, 
                                  dec_bin[-1] + 0.5*dec_bin_delta)

        # set the ra and dec limit
        if ra_limit == None:
            ra_limit = [ra_bin_edges.min(), ra_bin_edges.max()]
        if dec_limit == None:
            dec_limit = [dec_bin_edges.min(), dec_bin_edges.max()]


        # read the nvss catalog, convert to Jy
        nvss = nvss_cat(nvss_catalog_path)
        nvss['flux'] *= 1.e-3
        nvss['flux_err'] *= 1.e-3

        # find the point sources 
        selected = np.array([True,]*nvss['ra'].shape[0])
        if flux_limit != None:
            if flux_limit[0] > nvss['flux'].min():
                if output:
                    print 'flux cut at lower limit'
            if flux_limit[1] < nvss['flux'].max():
                if output:
                    print 'flux cut at higher limit'
            selected = np.logical_and(selected,
                       np.logical_and(nvss['flux']>flux_limit[0], 
                                      nvss['flux']<flux_limit[1]))
        if ra_limit != None:
            if ra_limit[0] > nvss['ra'].astype('float32').min():
                if output:
                    print 'ra   cut at lower limit'
            if ra_limit[1] < nvss['ra'].astype('float32').max():
                if output:
                    print 'ra   cut at higher limit'
            selected = np.logical_and(selected,
                       np.logical_and(nvss['ra'].astype('float32')>ra_limit[0], 
                                      nvss['ra'].astype('float32')<ra_limit[1]))
        if dec_limit != None:
            if dec_limit[0] > nvss['dec'].astype('float32').min():
                if output:
                    print 'dec  cut at lower limit'
            if dec_limit[1] < nvss['dec'].astype('float32').max():
                if output:
                    print 'dec  cut at higher limit'
            selected = np.logical_and(selected,
                       np.logical_and(nvss['dec'].astype('float32')>dec_limit[0], 
                                      nvss['dec'].astype('float32')<dec_limit[1]))

        ra  = nvss['ra'][selected].astype('float32')
        dec = nvss['dec'][selected].astype('float32')
        flux= nvss['flux'][selected].astype('float32')
        flux_err= nvss['flux_err'][selected].astype('float32')
        name= nvss['nvss'][selected]

        # deep copy the imap, selecte the freqs, get the freqs averaged map
        imap = np.copy(self.imap)
        imap = np.ma.array(imap[self.freq_list, ...])
        imap[np.isnan(imap)] = np.ma.masked
        imap[np.isinf(imap)] = np.ma.masked
        imap = np.ma.mean(imap, axis=0)

        # find the point sources flux in our map
        ra_map_index  = np.digitize(ra,  ra_bin_edges) - 1
        dec_map_index = np.digitize(dec, dec_bin_edges) - 1
        in_the_map = np.logical_and( np.logical_and(ra_map_index>=0, 
                                                    ra_map_index<imap.shape[0]),
                                     np.logical_and(dec_map_index>=0,
                                                    dec_map_index<imap.shape[1]))
        ra_map_index = ra_map_index[in_the_map]
        dec_map_index = dec_map_index[in_the_map]

        flat_index = np.ravel_multi_index([ra_map_index.tolist(), 
                                           dec_map_index.tolist()], 
                                           imap.shape)
        
        flux_in_map = imap.flatten().take(flat_index)
        spec_in_map = self.imap[:,ra_map_index,dec_map_index]
        if self.nmap != None:
            nmap = np.copy(self.nmap)
            flux_err_in_map = nmap.flatten().take(flat_index) 

        flux_in_hipass = self.mapping_coord(hipass_map_path).flatten().take(flat_index)

        # get rescale and zerop
        if by_fitting:
            self.re_scale, self.re_zerop = self.cal_by_fitting(hipass_map_path,
                                                    point_sources = flat_index,
                                                    mask_pixels = 5)
        else:
            flux_range_hipass = flux_in_hipass[0] - flux_in_hipass[-1]
            flux_range_map = flux_in_map[0] - flux_in_map[-1]
            self.re_scale = flux_range_hipass/flux_range_map
            self.re_zerop = flux_in_hipass[0] - flux_in_map[0]* self.re_scale

        flux_in_map_rescaled = flux_in_map * self.re_scale
        spec_in_map_rescaled = spec_in_map * self.re_scale
        if self.nmap != None:
            flux_err_in_map_rescaled = flux_err_in_map*self.re_scale

        flux_in_map_rescaled = flux_in_map_rescaled + self.re_zerop
        spec_in_map_rescaled = spec_in_map_rescaled + self.re_zerop

        if output:
            print '--------------------------------------------------------'
            print " nvss | hipass | map | map rescale"
            for i in range(len(ra)):
                print '%s %7.3f, %7.3f, %7.4f | '%(name[i], ra[i], dec[i], flux[i]),
                print '%7.3f, %7.3f, %7.4f | '%(ra_bin[ra_map_index[i]], 
                                                dec_bin[dec_map_index[i]], 
                                                flux_in_hipass[i]),
                print '%7.3f, %7.3f, %7.4f | '%(ra_bin[ra_map_index[i]], 
                                                dec_bin[dec_map_index[i]], 
                                                flux_in_map[i]),
                print '%7.3f, %7.3f, %7.4f | '%(ra_bin[ra_map_index[i]], 
                                                dec_bin[dec_map_index[i]], 
                                                flux_in_map_rescaled[i])
            print '--------------------------------------------------------'
            print 'In the map, we got %d sources;'%flux_in_map.shape[0]
            print 'In NVSS,    we got %d sources.'%flux.shape[0]
            print 'In HIPASS,  we got %d sources.'%flux_in_hipass.shape[0]
            print '[%7.5f] x flux + [%7.5f]'%(self.re_scale, self.re_zerop)
            print '========================================================'
            print

        self.flux_in_nvss = flux
        self.flux_err_in_nvss = flux_err
        self.flux_in_hipass = flux_in_hipass
        self.flux_in_map = flux_in_map
        self.spec_in_map = spec_in_map
        self.flux_in_map_rescaled = flux_in_map_rescaled
        self.spec_in_map_rescaled = spec_in_map_rescaled
        if self.nmap != None:
            self.flux_err_in_map = flux_err_in_map
            self.flux_err_in_map_rescaled = flux_err_in_map_rescaled

    def mapping_coord(self, hipass_map_path):

        return mapping_coord(self.imap, hipass_map_path, self.degrade_function)

    def degrade_resolution(self, map):
        
        common_resolution = self.degrade_function

        map = common_resolution.apply(map)

        return map

