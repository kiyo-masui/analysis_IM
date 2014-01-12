#! /usr/bin/env python 

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

class PlotMap(object):
    def __init__(self, mapdict, imap=None, qmap=None, umap=None, freq_cut=[], 
                 degrade_factor=1.1, degrade_map=False, plot_size=(10,8),
                 nvss_catalog_path=None, hipass_map_path=None, clip_weight=False):

        if nvss_catalog_path == None:
            self.nvss_catalog_path = os.getenv('NVSS_PATH')
        else:
            self.nvss_catalog_path = nvss_catalog_path

        if hipass_map_path == None:
            self.hipass_map_path = os.getenv('HIPASS_PATH')
        else:
            self.hipass_map_path = hipass_map_path

        self.name = mapdict['name']

        if 'noise_name' in mapdict.keys():
            self.noise_name = mapdict['noise_name']
        else:
            self.noise_name = None

        self.degrade_factor=degrade_factor
        if imap is None:
            self.imap = algebra.make_vect(algebra.load(mapdict['imap']))
            self.nmap = algebra.make_vect(algebra.load(mapdict['nmap']))
            print self.imap.max(), self.imap.min()
            if not degrade_map:
                self.imap = self.degrade_resolution(self.imap, degrade_factor)

                self.nmap[self.nmap < 1.e-30] = 1.e-30
                noise = self.nmap
                noise = 1. / noise
                noise = self.degrade_resolution(noise, degrade_factor, noise=True)
                noise = self.degrade_resolution(noise, degrade_factor, noise=True)
                noise = 1. / noise
                self.nmap =  noise
                self.nmap[self.nmap < 1.e-20] = 0.
        else:
            self.imap = imap

        print self.name

        self.freq_list = tuple([ind for ind in range(self.imap.shape[0]) 
                                    if  ind not in freq_cut])
        #self.imap = self.imap[freq_list, ...]

        if 'mode' in mapdict.keys():
            self.mode = algebra.make_vect(algebra.load(mapdict['mode']))
        if 'vect' in mapdict.keys():
            self.svdmodpkl = cPickle.load(open(mapdict['vect']))

        self.imap_secs = []
        if 'imap_sec' in mapdict.keys():
            for i in range(len(mapdict['imap_sec'])):
                imap_sec = algebra.make_vect(algebra.load(mapdict['imap_sec'][i]))
                if degrade_factor != 0:
                    imap_sec = self.degrade_resolution(imap_sec, degrade_factor)
                self.imap_secs.append(imap_sec)

        self.nmap_secs = []
        if 'nmap_sec' in mapdict.keys():
            for i in range(len(mapdict['nmap_sec'])):
                nmap_sec = algebra.make_vect(algebra.load(mapdict['nmap_sec'][i]))
                if degrade_factor != 0:
                    nmap_sec[nmap_sec < 1.e-30] = 1.e-30
                    noise_sec = nmap_sec
                    noise_sec = 1. / noise_sec
                    noise_sec = self.degrade_resolution(noise_sec, degrade_factor, 
                                                        noise=True)
                    noise_sec = self.degrade_resolution(noise_sec, degrade_factor, 
                                                        noise=True)
                    noise_sec = 1. / noise_sec
                    nmap_sec = noise_sec
                    nmap_sec[nmap_sec < 1.e-20] = 0.
                self.nmap_secs.append(nmap_sec)

        self.re_scale = None
        self.re_zerop = None
        self.plot_size = plot_size
        self.clip_weight = clip_weight

    def saturate_weight(self, weight, percentile=50):
        weight_2d = np.mean(weight, axis=0)
        w_at_percentile = np.percentile(weight_2d, percentile)
        mask = weight_2d <= w_at_percentile
        for freq_index in range(weight.shape[0]):
            mweight = np.ma.array(weight[freq_index, ...], mask=mask)
            meanslice = mweight.mean()
            weight[freq_index, np.logical_not(mask)] = meanslice

        return weight

    def plot_mode_factor(self, num=0):

        print "plot fitting param"
        
        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        RA, DEC = np.meshgrid(ra, dec)

        mode = np.ma.array(self.mode[num, ...].T)
        mode[np.isnan(mode)] = np.ma.masked
        mode[np.isinf(mode)] = np.ma.masked

        max = np.ma.max(mode) * 0.5
        min = np.ma.min(mode) * 0.5
        if np.fabs(min) > np.fabs(max):
            max = np.fabs(min)
        else:
            min = -np.fabs(max)

        plt.figure(figsize=(8,6))
        plt.pcolormesh(RA, DEC, mode, vmax=max, vmin=min)
        plt.colorbar()
        plt.xlim(RA.min(), RA.max())
        plt.ylim(DEC.min(), DEC.max())
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec[deg]')
        plt.savefig('./png/%s_mode_factor_%d.png'%(self.name, num), format='png')

    def plot_mod(self, last_mode, num=0):

        print "plot modes"
        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        RA, DEC = np.meshgrid(ra, dec)
        
        if num >= self.mode.shape[0]:
            print "num error: should less than %d"%self.mode.shape[0]
            exit()

        mode_number = (last_mode + 1) - self.mode.shape[0] + num

        mode_vector = np.array(self.svdmodpkl[1])[mode_number-1]

        mode_factor = np.ma.array(self.mode[num, ...])
        mode_factor[np.isnan(mode_factor)] = np.ma.masked
        mode_factor[np.isinf(mode_factor)] = np.ma.masked

        mode = mode_vector[:, None, None] * mode_factor[None, :, :]
        mode = np.ma.mean(mode, axis=0)

        mode = mode.T
        max = np.ma.max(mode) * 0.5
        min = np.ma.min(mode) * 0.5

        plt.figure(figsize=self.plot_size)
        plt.pcolormesh(RA, DEC, mode, vmax=max, vmin=min)
        plt.colorbar()
        plt.xlim(RA.min(), RA.max())
        plt.ylim(DEC.min(), DEC.max())
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec[deg]')
        plt.savefig('./png/%s_mode_%d.png'%(self.name, mode_number), format='png')

        map = np.ma.array(self.imap[self.freq_list, ...])
        #map = np.ma.array(self.imap)
        map[np.isnan(map)] = np.ma.masked
        map[np.isinf(map)] = np.ma.masked
        map = np.ma.mean(map, axis=0)
        map = map.T

        mode = mode + map

        max = np.ma.max(mode) * 0.5
        min = np.ma.min(mode) * 0.5
        #if np.fabs(min) > np.fabs(max):
        #    max = np.fabs(min)
        #else:
        #    min = -np.fabs(max)

        plt.figure(figsize=self.plot_size)
        plt.pcolormesh(RA, DEC, mode, vmax=max, vmin=min)
        plt.colorbar()
        plt.xlim(RA.min(), RA.max())
        plt.ylim(DEC.min(), DEC.max())
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec[deg]')
        plt.savefig('./png/%s_mode_%d_plus_weighted_map.png'%(
                    self.name, mode_number), format='png')
   
    def make_mov(self):
        print "make mov"
        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        freq=self.imap.get_axis('freq')
        RA, DEC = np.meshgrid(ra, dec)
        freq = freq[list(self.freq_list)]
        map = np.ma.array(self.imap[self.freq_list, ...])
        max = np.ma.max(map) * 0.1
        min = np.ma.min(map) * 0.1
        if np.fabs(min) > np.fabs(max):
            max = np.fabs(min)
        else:
            min = -np.fabs(max)
        if not os.path.exists('./png/%s_mov/'%(self.name)):
            os.makedirs('./png/%s_mov/'%(self.name))
        print "frequence:",
        for i in range(map.shape[0]):
            print "%5.2fMHz -- "%(freq[i]/1.e6),
            sys.stdout.flush()
            plt.figure(figsize=self.plot_size)
            map = np.ma.array(self.imap[self.freq_list, ...])
            #map = np.ma.array(self.imap)
            map[np.isnan(map)] = np.ma.masked
            map[np.isinf(map)] = np.ma.masked
            map = map[i, ...]
            map = map.T
            plt.pcolormesh(RA, DEC, map, vmax=max, vmin=min)
            plt.colorbar()
            plt.xlim(RA.min(), RA.max())
            plt.ylim(DEC.min(), DEC.max())
            plt.xlabel('RA [deg]')
            plt.ylabel('Dec[deg]')
            plt.title('%s %5.2fMHz'%(self.name, freq[i]/1.e6))
            plt.savefig('./png/%s_mov/%s_f%03d.png'%(
                        self.name, self.name, i), format='png')
        print 'end'
    
    def hipass_cat(self, plot=True, cut=True):

        hipass_data, hipass_ra_bin_edges, hipass_dec_bin_edges, hipass_ra, hipass_dec =\
                cal_map.hipass_cat(self.hipass_map_path)

        if not cut:
            return hipass_data, hipass_ra_bin_edges, hipass_dec_bin_edges,\
                   hipass_ra, hipass_dec

        ra_bin = self.imap.get_axis('ra')
        dec_bin= self.imap.get_axis('dec')

        #print hipass_ra.min(), ra_bin.min()
        #print hipass_ra.max(), ra_bin.max()
        #print hipass_dec.min(), dec_bin.min()
        #print hipass_dec.max(), dec_bin.max()

        ra_range = np.digitize([ra_bin.min(), ra_bin.max()], 
                                hipass_ra_bin_edges) - 1
        dec_range = np.digitize([dec_bin.min(), dec_bin.max()], 
                                 hipass_dec_bin_edges) - 1
        #print ra_range
        #print dec_range

        hipass_ra_bin_edges = hipass_ra_bin_edges[ra_range[1]: ra_range[0]+1]
        hipass_dec_bin_edges = hipass_dec_bin_edges[dec_range[0]: dec_range[1]+1]
        hipass_data = hipass_data[dec_range[0]:dec_range[1],ra_range[1]:ra_range[0]]

        hipass_ra = hipass_ra[ra_range[1]: ra_range[0]]
        hipass_dec = hipass_dec[dec_range[0]: dec_range[1]]

        if plot:
        
            RA, DEC = np.meshgrid(hipass_ra, hipass_dec)
            plt.figure(figsize=self.plot_size)

            max = np.ma.max(hipass_data) #* 0.5
            min = np.ma.min(hipass_data) #* 0.5

            plt.pcolormesh(RA, DEC, hipass_data, vmax=max, vmin=min)
            plt.colorbar()

            nvss = cat_map.nvss_cat(self.nvss_catalog_path)
            nvss['flux'] *= 1.e-3
            selected = nvss['flux'] > .5
            ra  = nvss['ra'][selected]
            dec = nvss['dec'][selected]
            flux= nvss['flux'][selected]

            #for i in range(len(flux)):
            #    print flux[i],
            #    print ra[i],
            #    print dec[i]


            cmap = plt.get_cmap('jet')
            norm = plt.normalize(min, max)
            plt.scatter(ra, dec, s=40, c=cmap(norm(flux)), edgecolor='k' )


            plt.axes().set_aspect('equal')
            plt.xlim(ra_bin.min(), ra_bin.max())
            plt.ylim(dec_bin.min(), dec_bin.max())
            plt.xlabel('RA [deg]')
            plt.ylabel('Dec[deg]')
            plt.savefig('./png/%s_hipass_image.png'%self.name )
            #plt.show()
        else:
            return hipass_data, hipass_ra_bin_edges, hipass_dec_bin_edges,\
                   hipass_ra, hipass_dec

    def nvss_in_hipass(self, source_ra, source_dec):
        
        hipass_data, ra_bin_edges, dec_bin_edges, ra, dec\
            = self.hipass_cat(plot=False)

        ra_map_index  = np.digitize(source_ra,  ra_bin_edges) - 1
        dec_map_index = np.digitize(source_dec, dec_bin_edges) - 1
        in_the_map = np.logical_and( np.logical_and(ra_map_index>=0, 
                                     ra_map_index<hipass_data.shape[1]),
                                     np.logical_and(dec_map_index>=0,
                                     dec_map_index<hipass_data.shape[0]))
        ra_map_index = ra_map_index[in_the_map]
        dec_map_index = dec_map_index[in_the_map]


        flatten_index = np.ravel_multi_index([dec_map_index.tolist(), 
                                              ra_map_index.tolist(),], 
                                              hipass_data.shape)
        
        flux_in_map = hipass_data.flatten().take(flatten_index)

        return flux_in_map
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

    def degrade_resolution(self, map, degrade_factor=1.1, noise=False):

        common_resolution = self.degrade_function(degrade_factor=degrade_factor)

        if noise:
            map = common_resolution.apply(map, mode='constant', cval=1.e30)
        else:
            map = common_resolution.apply(map)

        return map

    def check_nvss(self, flux_limit=None, ra_limit=None, dec_limit=None, 
                   rescale=True):

        if len(self.imap_secs) != 0:
            # using var of section maps as the err
            imap_secs = np.ma.array(self.imap_secs)[:,self.freq_list, ...]
            err  = np.std(imap_secs, axis=0)
            err *= np.sqrt(imap_secs.shape[0]/(imap_secs.shape[0]-1.))
            err[np.isnan(err)] = np.ma.masked
            err[np.isinf(err)] = np.ma.masked
            err = np.sqrt(np.ma.mean(err**2., axis=0))
        else:
            #  using noise_inv as the err
            err = np.ma.array(self.nmap[self.freq_list, ...])
            err[np.isnan(err)] = np.ma.masked
            err[np.isinf(err)] = np.ma.masked
            err[err == 0.] = np.inf
            err = (1./err)
            err = np.sqrt(np.ma.mean(err, axis=0))

        calabrator = cal_map.CalibrateMap(self.imap, self.freq_list, nmap=err, 
                weight=self.nmap, degrade_map=True,
                degrade_function = self.degrade_function())
        calabrator.cal_by_point_sources(self.nvss_catalog_path, self.hipass_map_path,
                flux_limit=flux_limit, ra_limit=ra_limit, dec_limit=dec_limit)
        flux                     = calabrator.flux_in_nvss
        flux_err                 = calabrator.flux_err_in_nvss
        flux_in_hipass           = calabrator.flux_in_hipass
        flux_in_map              = calabrator.flux_in_map
        spec_in_map              = calabrator.spec_in_map
        flux_err_in_map          = calabrator.flux_err_in_map
        flux_in_map_rescaled     = calabrator.flux_in_map_rescaled
        spec_in_map_rescaled     = calabrator.spec_in_map_rescaled
        flux_err_in_map_rescaled = calabrator.flux_err_in_map_rescaled

        if rescale:
            # save the limit for section map rescaling in plotmap
            self.flux_limit = flux_limit
            self.ra_limit = ra_limit
            self.dec_limit = dec_limit

            self.re_scale = calabrator.re_scale
            self.re_zerop = calabrator.re_zerop

        fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=self.plot_size)
        plt.subplots_adjust(left=0.06, right=0.95, bottom=0.15, top=0.95, wspace=0.001)

        ax[0].errorbar(np.arange(flux_in_map.shape[0]), flux, flux_err, label='nvss',
                       fmt='ro')
        if rescale:
            ax[0].errorbar(np.arange(flux_in_map.shape[0]), flux_in_map_rescaled, 
                           flux_err_in_map_rescaled, label='map', fmt='gs')
        else:
            ax[0].errorbar(np.arange(flux_in_map.shape[0]), flux_in_map, 
                           flux_err_in_map, label='map', fmt='gs')
        ax[0].scatter(np.arange(flux_in_hipass.shape[0]), flux_in_hipass,
                     label='hipass', s=40, c='b', )
        if rescale:
            ax[0].text(0.1, 0.8,
                       r'rescaled by [%7.5f] $\times$ flux + [%7.5f]'\
                       %(self.re_scale, self.re_zerop),
                       horizontalalignment='left',
                       verticalalignment='center',
                       transform=ax[0].transAxes)
        ax[0].set_xlabel('Source No.')
        ax[0].set_ylabel('Flux [Jy]')
        ax[0].legend(frameon=False, scatterpoints=1)
        ax[0].set_xlim(xmin=-0.5)

        spec_axis = self.imap.get_axis('freq')
        for i in range(flux_in_map.shape[0]):
            if rescale:
                ax[1].plot(spec_axis/1.e6, spec_in_map_rescaled[:, i])
            else:
                ax[1].plot(spec_axis/1.e6, spec_in_map[:, i])
        ax[1].set_xlabel('Frequency [MHz]')
        ax[1].set_xlim(xmin=spec_axis.min()/1.e6, xmax=spec_axis.max()/1.e6)

        if rescale:
            plt.savefig('./png/%s_check_nvss_rescaled.png'%(self.name), format='png')
        else:
            plt.savefig('./png/%s_check_nvss.png'%(self.name), format='png')

    def rescale_by_hipass(self):

        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')

        map = np.copy(self.imap)
        map = np.ma.array(map[self.freq_list, ...])
        map[np.isnan(map)] = np.ma.masked
        map[np.isinf(map)] = np.ma.masked
        map = np.ma.mean(map, axis=0)

        # get the same field from hipass
        hipass_data, ra_bin_edges, dec_bin_edges, ra, dec\
            = self.hipass_cat(plot=False)
        scale_hipass = np.max(hipass_data)-np.min(hipass_data)
        scale_map = np.max(map) - np.min(map)
        self.re_scale = scale_hipass/scale_map
        self.re_zerop = np.max(hipass_data) - np.max(map*self.re_scale)

    def get_edges(self, x):

        dx = x[1] - x[0]
        x = x - 0.5 * dx
        x = np.append(x, x[-1] + dx)

        return x

    def mapping_coord(self, plot=False):

        hipass_data_new = cal_map.mapping_coord(self.imap, self.hipass_map_path, 
                degrade_function = self.degrade_function())

        if plot:
            plt.figure(figsize=self.plot_size)
            max = 0.5*np.max(hipass_data_new)
            min = 0.5*np.min(hipass_data_new)
            ra = self.imap.get_axis('ra')
            dec= self.imap.get_axis('dec')
            RA, DEC = np.meshgrid( self.get_edges(ra),  self.get_edges(dec))
            plt.pcolormesh(RA, DEC, hipass_data_new.T, vmax=max, vmin=min)
            plt.colorbar()

            nvss = cal_map.nvss_cat(self.nvss_catalog_path)
            nvss['flux'] *= 1.e-3
            selected = nvss['flux'] > .5
            ra  = nvss['ra'][selected]
            dec = nvss['dec'][selected]
            flux= nvss['flux'][selected]

            cmap = plt.get_cmap('jet')
            norm = plt.normalize(min, max)
            plt.scatter(ra, dec, s=40, c=cmap(norm(flux)), edgecolor='k' )

            plt.axes().set_aspect('equal')
            plt.xlim(RA.min(), RA.max())
            plt.ylim(DEC.min(), DEC.max())
            plt.xlabel('RA [deg]')
            plt.ylabel('Dec[deg]')
            plt.savefig('./png/%s_hipass_mapping.png'%(self.name), format='png')

        return hipass_data_new

    def plot_map(self, num=None, with_nvss=False, diff=False):
        print "plot maps"


        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        RA, DEC = np.meshgrid( self.get_edges(ra),  self.get_edges(dec))

        map_list = [self.imap, ]  + self.imap_secs
        sec_list = ['',]
        if len(self.imap_secs) == 1:
            sec_list += ['_A']
        elif len(self.imap_secs) == 4:
            sec_list += ['_A', '_B', '_C', '_D']
        elif len(self.imap_secs) == 5:
            sec_list += ['_A', '_B', '_C', '_D', '_E']
        elif len(self.imap_secs) == 3:
            sec_list += ['_A', '_B', '_C']
        elif len(self.imap_secs) == 6:
            sec_list += ['_AB', '_AC', '_BA', '_BC', '_CA', '_CB']
        elif len(self.imap_secs) == 12:
            sec_list += ['_AB', '_AC', '_AD', '_BA', '_BC', '_BD', 
                         '_CA', '_CB', '_CD', '_DA', '_DB', '_DC']
        elif len(self.imap_secs) == 20:
            sec_list += ['_AB', '_AC', '_AD', '_AE', 
                         '_BA', '_BC', '_BD', '_BE',
                         '_CA', '_CB', '_CD', '_CE',
                         '_DA', '_DB', '_DC', '_DE',
                         '_EA', '_EB', '_EC', '_ED',]
        if self.noise_name != None:
            map_list = map_list + [self.nmap, ] + self.nmap_secs
            sec_list = sec_list + sec_list

        for i in range(len(map_list)):

            if i < len(self.imap_secs) + 1:
                print '\t', self.name
            else:
                print '\t', self.noise_name

            suffix = ''

            map = np.copy(map_list[i])
            if i >= len(self.imap_secs) + 1 and self.clip_weight:
                map = self.saturate_weight(map)
            map = np.ma.array(map[self.freq_list, ...])
            map[np.isnan(map)] = np.ma.masked
            map[np.isinf(map)] = np.ma.masked


            if num==None:
                map = np.ma.mean(map, axis=0)
                suffix = '_average'
                print '\taverage',
            else:
                map = map[num, ...]
                suffix = '_freq%03d'%num
                print '\tfreq%03d'%num,
            map = map.T

            suffix += '%s'%sec_list[i]
            print '%s'%sec_list[i].replace('_', ' '),

            if i >= len(self.imap_secs) + 1 and self.clip_weight:
                suffix += '_clip50'

            if self.re_scale != None:
                if i < len(self.imap_secs) + 1:
                    calabrator = cal_map.CalibrateMap(map_list[i], self.freq_list, 
                            weight=map_list[i + len(self.imap_secs) + 1], 
                            degrade_map=True, degrade_function = self.degrade_function())
                    calabrator.cal_by_point_sources(self.nvss_catalog_path, 
                            self.hipass_map_path, flux_limit=self.flux_limit, 
                            ra_limit=self.ra_limit, dec_limit=self.dec_limit,
                            output=False)
                    re_scale = calabrator.re_scale
                    re_zerop = calabrator.re_zerop
                else:
                    re_scale = 1.
                    re_zerop = 0.

                map *= re_scale
                map += re_zerop
                suffix += '_rescaled'
                print ' rescaled x[%5.3f], +[%5.3f]'%(re_scale, re_zerop),

            if diff: 
                if i < len(self.imap_secs) + 1:
                    hipass = self.mapping_coord()
                    map -= hipass.T
                    suffix += '_diff'
                    print ' diff'
                else:
                    # do not plot the diff map for noise
                    print 'do not plot diff'
                    print
                    continue
            else:
                print 
            print

            max = np.max(map) * 0.5
            min = np.min(map) * 0.5

            map = np.ma.array(map)
            map[map == 0] = np.ma.masked
            plt.figure(figsize=self.plot_size)
            plt.pcolormesh(RA, DEC, map, vmax=max, vmin=min)
            plt.colorbar()

            if with_nvss:
                nvss = cal_map.nvss_cat(self.nvss_catalog_path)
                nvss['flux'] *= 1.e-3
                selected = nvss['flux'] > .3
                ra  = nvss['ra'][selected]
                dec = nvss['dec'][selected]
                flux= nvss['flux'][selected]

                cmap = plt.get_cmap('jet')
                norm = plt.normalize(min, max)
                if diff:
                    plt.scatter(ra, dec, s=40, c='none', edgecolor='k' )
                else:
                    #plt.scatter(ra, dec, s=40, c=cmap(norm(flux)), edgecolor='k' )
                    plt.scatter(ra, dec, s=40, c='none', edgecolor='0.5' )

            plt.axes().set_aspect('equal')
            plt.xlim(RA.min(), RA.max())
            plt.ylim(DEC.min(), DEC.max())
            plt.ylabel('Dec[deg]')
            plt.xlabel('RA [deg]')
            if i < len(self.imap_secs) + 1:
                plt.savefig('./png/%s%s.png'%(self.name, suffix), format='png')
                plt.savefig('./png/%s%s.eps'%(self.name, suffix), format='eps')
            else:
                plt.savefig('./png/%s%s.png'%(self.noise_name, suffix), format='png')
                plt.savefig('./png/%s%s.eps'%(self.noise_name, suffix), format='eps')
            plt.close()
            del map
            if with_nvss:
                del ra, dec
            gc.collect()
        

if __name__=='__main__':
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_themselves/"
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_1gwj/IQUmap_clean_withIxIsvd/"
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQUV_extend_legendre_modes_0gwj/Emap_clean_themselves/"
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQUV_extend_legendre_modes_0gwj/Emap_clean_themselves/"

    mapdict = {}
    #maproot = "/home/ycli/data/map_result/maps/fir_RA+27_Dec-30_parkes_2010_10_24-28_over_flagged/"
    #mapdict['imap'] = maproot + "fir_RA+27_Dec-30_parkes_2010_10_24-28_dirty_map_I_1316.npy"
    #mapdict['name'] = "fir_RA+27_Dec-30_parkes_2010_10_24-28_dirty_map"

    #mapdict['imap'] = maproot + "fir_RA+27_Dec-30_parkes_2010_10_24-28_clean_map_I_1316.npy"
    #mapdict['name'] = "fir_RA+27_Dec-30_parkes_2010_10_24-28_clean_map"

    #maproot = "/Users/ycli/DATA/parkes/maps_by_chris/pks_27n30_10by7/"
    #mapdict['imap'] = maproot + "test_allbeams_27n30_10by7_clean_map_I_1315.npy"
    #mapdict['name'] = "test_allbeams_27n30_10by7_clean_map"

    #maproot = "/home/ycli/data/map_result/maps/parkes/"
    #mapdict['imap'] = maproot + "fir_RA+27_Dec-30_parkes_2010_10_24-28_dirty_map_I_1316.npy"
    #mapdict['name'] = "fir_RA+27_Dec-30_parkes_2010_10_24-28_dirty_map"
    #mapdict['imap'] = maproot + "fir_RA+27_Dec-30_parkes_2010_10_24-28_clean_map_I_1316.npy"
    #mapdict['name'] = "fir_RA+27_Dec-30_parkes_2010_10_24-28_clean_map"
    #a = plotmap(mapdict)
    #a.plot_map()

    #exit()

    #----------------------

    #maproot = "/home/ycli/data/cln_result/15hr_AA_fine_freq_11conv/Emap_clean_themselves/"

    #for i in range(5):
    #    mapdict['imap'] = maproot + "sec_A_cleaned_clean_map_I_with_A_%dmodes.npy"%i
    #    mapdict['mode'] = maproot + "sec_A_modes_clean_map_I_with_A_%dmodes.npy"%i
    #    mapdict['name'] = 'AA_cleaned_clean_%d'%i
    #    mapdict['vect'] =  maproot + "SVD_pair_A_with_A.pkl"
    #    a = plotmap(mapdict,)
    #    a.plot_map()
    #    if i ==0:
    #        continue
    #    a.plot_mod(i)

    #exit()


    #----------------------

    #maproot = "/home/ycli/data/cln_result/DRACO_AA_dirty_map/Emap_clean_themselves/"
    #maproot = "/home/ycli/data/cln_result/3C286/Emap_clean_themselves/"

    #for i in range(5):
    #    mapdict['imap'] = maproot + "sec_A_cleaned_clean_map_I_with_A_%dmodes.npy"%i
    #    mapdict['mode'] = maproot + "sec_A_modes_clean_map_I_with_A_%dmodes.npy"%i
    #    mapdict['name'] = 'AA_cleaned_clean_%d'%i
    #    a = plotmap(mapdict, freq_cut = [1,2,3,4,5,59,60,61,62,63])
    #    a.plot_map()
    #    if i ==0:
    #        continue
    #    a.plot_mod()

    #exit()


    #----------------------

    #maproot = "/home/ycli/data/cln_result/PKS_ABC_allbeams_10by7_freq_cut/Emap_clean_themselves/"

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_27n30_10by7_clean_map_I'
    #mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_A_clean_map_I_1315.npy',
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_B_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_C_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_D_clean_map_I_1315.npy', 
    #                      ]

    #mode = 2
    #mapdict['imap'] = '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/combined_clean_map_%dmodes.npy'%mode
    #mapdict['nmap'] = '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/combined_clean_weight_%dmodes.npy'%mode
    #mapdict['name'] = 'test_allbeams_27n30_10by7_clean_map_I_%dmode_1pt1_cov'%mode
    #mapdict['imap_sec'] = ['/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_A_cleaned_clean_map_I_with_B_%dmodes.npy'%mode,
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_A_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_A_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_B_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_B_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_B_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_C_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_C_cleaned_clean_map_I_with_B_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_C_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_D_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_D_cleaned_clean_map_I_with_B_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_27n30_10by7_ABCD_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_D_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
    #                      ]

    mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
    mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
    mapdict['name'] = 'test_allbeams_27n30_10by7_noise_inv_I'
    mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_A_noise_weight_I_1315.npy',
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_B_noise_weight_I_1315.npy', 
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_C_noise_weight_I_1315.npy', 
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_D_noise_weight_I_1315.npy', 
                          ]

    #mapdict['imap'] = '/home/ycli/data/map_result/maps/parkes/fir_p2750n2950_parkes_2010_10_24-28_clean_map_I_1316.npy'
    #mapdict['nmap'] = '/home/ycli/data/map_result/maps/parkes/fir_p2750n2950_parkes_2010_10_24-28_noise_weight_I_1316.npy'
    #mapdict['name'] = 'p2750n2950_parkes_2010_10_24-28_clean_map_I_1316_1pt1_cov'

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_n03n30_10by7_0627pixel/test_allbeams_n03n30_10by7_0627pixel_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_n03n30_10by7_0627pixel/test_allbeams_n03n30_10by7_0627pixel_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_n03n30_10by7_0627pixel_clean_map_I_1pt1_cov'
    #mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_n03n30_10by7_0627pixel/test_allbeams_n03n30_10by7_0627pixel_A_clean_map_I_1315.npy',
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_n03n30_10by7_0627pixel/test_allbeams_n03n30_10by7_0627pixel_B_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_n03n30_10by7_0627pixel/test_allbeams_n03n30_10by7_0627pixel_C_clean_map_I_1315.npy', 
    #                      ]

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_07n30_10by7_0627pixel/test_allbeams_07n30_10by7_0627pixel_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_07n30_10by7_0627pixel/test_allbeams_07n30_10by7_0627pixel_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_07n30_10by7_0627pixel_clean_map_I_1pt1_cov'
    #mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_07n30_10by7_0627pixel/test_allbeams_07n30_10by7_0627pixel_A_clean_map_I_1315.npy',
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_07n30_10by7_0627pixel/test_allbeams_07n30_10by7_0627pixel_B_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_07n30_10by7_0627pixel/test_allbeams_07n30_10by7_0627pixel_C_clean_map_I_1315.npy', 
    #                      ]

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_17n30_10by7_0627pixel_clean_map_I_1pt1_cov'
    #mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_A_clean_map_I_1315.npy',
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_B_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_C_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_D_clean_map_I_1315.npy', 
    #                      ]

    #mapdict['imap'] = '/home/ycli/data/parkes/maps/parkes/fir_RA+27_Dec-30_parkes_2010_10_24-28_clean_map_I_1316.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps/parkes/fir_RA+27_Dec-30_parkes_2010_10_24-28_noise_weight_I_1316.npy'
    #mapdict['name'] = 'fir_RA+27_Dec-30_parkes_2010_10_24-28_clean_map_I_1316_1pt1_cov'
    #mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_A_clean_map_I_1315.npy',
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_B_clean_map_I_1315.npy', 
    #                       '/home/ycli/data/parkes/maps_by_chris/pks_17n30_10by7_0627pixel/test_allbeams_17n30_10by7_0627pixel_C_clean_map_I_1315.npy', 
    #                      ]

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_A_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_A_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_27n30_10by7_A_clean_map_I_1pt1_cov'

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_B_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_B_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_27n30_10by7_B_clean_map_I_1pt1_cov'

    #mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_C_clean_map_I_1315.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_C_noise_weight_I_1315.npy'
    #mapdict['name'] = 'test_allbeams_27n30_10by7_C_clean_map_I_1pt1_cov'

    #mapdict['imap'] = '/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_24-28_clean_map_I_1316.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_24-28_noise_weight_I_1316.npy'
    #mapdict['name'] = 'fir_p3500n3000_parkes_2010_10_24-28_clean_map'
    #mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_24_clean_map_I_1316.npy',
    #                       '/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_25_clean_map_I_1316.npy', 
    #                       '/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_26_clean_map_I_1316.npy', 
    #                      ]

    #mode = 5
    #mapdict['imap'] = '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/combined_clean_map_%dmodes.npy'%mode
    #mapdict['nmap'] = '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/combined_clean_weight_%dmodes.npy'%mode
    #mapdict['name'] = 'fir_p3500n3000_parkes_2010_10_24-28_clean_map_%dmode'%mode
    #mapdict['imap_sec'] = ['/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_A_cleaned_clean_map_I_with_B_%dmodes.npy'%mode,
    #                       '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_A_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_B_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_B_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_C_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
    #                       '/home/ycli/data/cln_result/PKS_p3500n3000_parkes_2010_10_24-28_ABC_freq_cut_1pt1_cov_cal/Emap_clean_themselves/sec_C_cleaned_clean_map_I_with_B_%dmodes.npy'%mode, 
    #                      ]

    #degrade_factor = 1.1
    degrade_factor = 0
    if degrade_factor != 0:
        mapdict['name'] += '_1pt1_cov'

    a = plotmap(mapdict, freq_cut = [0,1,2,3,4,5,59,60,61,62,63], 
                degrade_factor=degrade_factor, plot_size=(15, 6))
    #a.mapping_coord(plot=True)
    #a.hipass_cat()
    a.plot_map(with_nvss=True)
    #a.plot_map(with_nvss=True, diff=True)
    #a.check_nvss(flux_limit=[0.5, 100], rescale=False)
    #a.check_nvss(flux_limit=[0.5, 100])
    #a.plot_map(with_nvss=True)
    #a.plot_map(with_nvss=True, diff=True)

    exit()

    for i in range(13):
        if i == 3: continue
        maproot = '/home/ycli/data/map_result/maps/parkes/'
        mapdict['imap'] = maproot + "fir_p2750n2950_parkes_2010_10_24-28_beam_%d_clean_map_I_1316.npy"%i
        mapdict['nmap'] = maproot + "fir_p2750n2950_parkes_2010_10_24-28_beam_%d_noise_weight_I_1316.npy"%i
        mapdict['name'] = 'p2750n2950_parkes_2010_10_24-28_beam_%d_clean_map_I_1316.npy'%i
        a = plotmap(mapdict, freq_cut = [0,1,2,3,4,5,59,60,61,62,63])
        a.check_nvss(flux_limit=[1., 10], degrade=True)
        a.plot_map()

    exit()

    for i in range(0, 5):
        mapdict['imap'] = maproot + "sec_B_cleaned_clean_map_I_with_C_%dmodes.npy"%i
        mapdict['nmap'] = maproot + "sec_B_cleaned_noise_inv_I_with_C_%dmodes.npy"%i
        #mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
        mapdict['mode'] = maproot + "sec_B_modes_clean_map_I_with_C_%dmodes.npy"%i
        mapdict['name'] = 'BC_cleaned_clean_%d'%i
        mapdict['vect'] =  maproot + "SVD_pair_B_with_C.pkl"
        a = plotmap(mapdict, freq_cut = [0,1,2,3,4,5,59,60,61,62,63])
        #a.plot_map(num=10)
        #a.make_mov()
        a.check_nvss(flux_limit=[1., 10])
        a.plot_map()
        #a.plot_mod(i)

    exit()

    #----------------------

    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_1modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_1modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_1modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_2'
    #a = plotIQU(mapdict)
    #imap1 = a.imap
    #qmap1 = a.qmap
    #umap1 = a.umap
    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_2modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_2modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_2modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_2'
    #b = plotIQU(mapdict)
    #imap2 = b.imap
    #qmap2 = b.qmap
    #umap2 = b.umap


    #imap = imap1 - imap2
    #qmap = qmap1 - qmap2
    #umap = umap1 - umap2

    #mapdict['name'] = 'iqu_cleaned_clean_2-1'
    #c = plotIQU(mapdict, imap, qmap, umap)

    #c.plotiqu(150)


