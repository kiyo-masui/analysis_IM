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

from map import beam

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
        plt.savefig('./png/%s_mode_%d.png'%(self.name, num), format='png')
        #plt.show()

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
        #mode = mode[20,...]

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
        plt.savefig('./png/%s_mode_%d_plus_weighted_map.png'%(self.name, mode_number), format='png')

   
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
    
    def nvss_cat(self):

        nvss_cat = np.genfromtxt('/home/ycli/data/nvss/nvss_catlog.dat', 
                   delimiter='|', 
                   usecols = tuple(range(1,8)),
                   dtype="S19, S12, S12, f4, f4, f4, f4",
                   names=['nvss','ra','dec','ra_err','dec_err',
                          'flux','flux_err'])
        
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

    def hipass_cat(self, plot=True, cut=True):

        hdulist = pyfits.open('/home/ycli/data/hipass/CHIPASS_0040-30_WGTMED1.continuum.fits')

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

            max = np.ma.max(hipass_data) * 0.5
            min = np.ma.min(hipass_data) * 0.5

            plt.pcolormesh(RA, DEC, hipass_data, vmax=max, vmin=min)
            plt.colorbar()

            nvss = self.nvss_cat()
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


        #print hipass_ra.shape
        #print hipass_dec.shape
        #print hipass_data.shape


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

        print 
        print 'hipass'
        for i in range(len(ra_map_index)):
            print ra[ra_map_index[i]], dec[dec_map_index[i]], flux_in_map[i]
        print 

        return flux_in_map

    def degrade_resolution(self, map, degrade_factor=1.1):
        
        freq_data = sp.array([1250, 1275, 1300, 1325, 1350, 1430], dtype=float)
        beam_data = sp.array([14.4, 14.4, 14.4, 14.4, 14.4, 14.4])/60. 
        beam_data = beam_data*1420/freq_data
        freq_data *= 1.0e6
        beam_diff = sp.sqrt(max(degrade_factor*beam_data)**2-(beam_data)**2)
        common_resolution = beam.GaussianBeam(beam_diff, freq_data)

        map = common_resolution.apply(map)

        return map

    def check_nvss(self, flux_limit=None, ra_limit=None, dec_limit=None, 
                   rescale=True):

        ra_bin = self.imap.get_axis('ra')
        dec_bin= self.imap.get_axis('dec')

        ra_bin_delta = ra_bin[1] - ra_bin[0]
        ra_bin_edges = np.append(ra_bin - 0.5*ra_bin_delta, 
                                 ra_bin[-1] + 0.5*ra_bin_delta)
        dec_bin_delta = dec_bin[1] - dec_bin[0]
        dec_bin_edges = np.append(dec_bin - 0.5*dec_bin_delta, 
                                  dec_bin[-1] + 0.5*dec_bin_delta)

        if ra_limit == None:
            ra_limit = [ra_bin_edges.min(), ra_bin_edges.max()]
        if dec_limit == None:
            dec_limit = [dec_bin_edges.min(), dec_bin_edges.max()]


        nvss = self.nvss_cat()
        nvss['flux'] *= 1.e-3
        nvss['flux_err'] *= 1.e-3

        selected = np.array([True,]*nvss['ra'].shape[0])
        if flux_limit != None:
            if flux_limit[0] > nvss['flux'].min():
                print 'flux cut at lower limit'
            if flux_limit[1] < nvss['flux'].max():
                print 'flux cut at higher limit'
            selected = np.logical_and(selected,
                       np.logical_and(nvss['flux']>flux_limit[0], 
                                      nvss['flux']<flux_limit[1]))
        if ra_limit != None:
            if ra_limit[0] > nvss['ra'].astype('float32').min():
                print 'ra   cut at lower limit'
            if ra_limit[1] < nvss['ra'].astype('float32').max():
                print 'ra   cut at higher limit'
            selected = np.logical_and(selected,
                       np.logical_and(nvss['ra'].astype('float32')>ra_limit[0], 
                                      nvss['ra'].astype('float32')<ra_limit[1]))
        if dec_limit != None:
            if dec_limit[0] > nvss['dec'].astype('float32').min():
                print 'dec  cut at lower limit'
            if dec_limit[1] < nvss['dec'].astype('float32').max():
                print 'dec  cut at higher limit'
            selected = np.logical_and(selected,
                       np.logical_and(nvss['dec'].astype('float32')>dec_limit[0], 
                                      nvss['dec'].astype('float32')<dec_limit[1]))

        ra  = nvss['ra'][selected].astype('float32')
        dec = nvss['dec'][selected].astype('float32')
        flux= nvss['flux'][selected].astype('float32')
        flux_err= nvss['flux_err'][selected].astype('float32')
        name= nvss['nvss'][selected]

        map = np.copy(self.imap)
        map = np.ma.array(map[self.freq_list, ...])
        map[np.isnan(map)] = np.ma.masked
        map[np.isinf(map)] = np.ma.masked
        map = np.ma.mean(map, axis=0)

        if len(self.imap_secs) != 0:
            # using var of section maps as the err
            imap_secs = np.ma.array(self.imap_secs)[:,self.freq_list, ...]
            err  = np.std(imap_secs, axis=0)
            err *= np.sqrt(imap_secs.shape[0]/(imap_secs.shape[0]-1.))
            print "err shape:", err.shape
            err[np.isnan(map)] = np.ma.masked
            err[np.isinf(map)] = np.ma.masked
            err = np.sqrt(np.ma.mean(err**2., axis=0))
        else:
            #  using noise_inv as the err
            err = np.ma.array(self.nmap[self.freq_list, ...])
            err[np.isnan(map)] = np.ma.masked
            err[np.isinf(map)] = np.ma.masked
            err[err == 0.] = np.inf
            err = (1./err)**2
            err = np.sqrt(np.ma.mean(err, axis=0))

        ra_map_index  = np.digitize(ra,  ra_bin_edges) - 1
        dec_map_index = np.digitize(dec, dec_bin_edges) - 1
        in_the_map = np.logical_and( np.logical_and(ra_map_index>=0, 
                                                    ra_map_index<map.shape[0]),
                                     np.logical_and(dec_map_index>=0,
                                                    dec_map_index<map.shape[1]))
        ra_map_index = ra_map_index[in_the_map]
        dec_map_index = dec_map_index[in_the_map]

        flatten_index = np.ravel_multi_index([ra_map_index.tolist(), 
                                              dec_map_index.tolist()], 
                                              map.shape)
        
        flux_in_map = map.flatten().take(flatten_index)
        flux_err_in_map = err.flatten().take(flatten_index)

        spec_in_map = self.imap[:,ra_map_index,dec_map_index]
        spec_axis = self.imap.get_axis('freq')

        #flux_in_hipass = self.nvss_in_hipass(ra, dec)
        flux_in_hipass = self.mapping_coord().flatten().take(flatten_index)

        # re-scale map to hipass scale
        if self.re_scale != None:
            flux_in_map_rescaled = flux_in_map * self.re_scale
            flux_in_map_rescaled = flux_in_map + self.re_zerop
            spec_in_map_rescaled = spec_in_map * self.re_scale
            spec_in_map_rescaled = spec_in_map + self.re_zerop
            rescale = True
        elif rescale:
            print "rescale by bright sources"
            #flux_range_hipass = flux_in_hipass.max() - flux_in_hipass.min()
            #flux_range_map = flux_in_map.max() - flux_in_map.min()
            flux_range_hipass = flux_in_hipass[0] - flux_in_hipass[-1]
            flux_range_map = flux_in_map[0] - flux_in_map[-1]
            self.re_scale = flux_range_hipass/flux_range_map
            flux_in_map_rescaled = flux_in_map * self.re_scale
            spec_in_map_rescaled = spec_in_map * self.re_scale
            #self.re_zerop = flux_in_hipass.min() - flux_in_map.min()
            self.re_zerop = flux_in_hipass[0] - flux_in_map_rescaled[0]
            flux_in_map_rescaled = flux_in_map_rescaled + self.re_zerop
            spec_in_map_rescaled = spec_in_map_rescaled + self.re_zerop
            flux_err_in_map_rescaled = flux_err_in_map*self.re_scale

        print '--------------------------------------------------------'
        if rescale:
            print " nvss | hipass | map | map rescale"
        else:
            print " nvss | hipass | map "
        for i in range(len(ra)):
            print '%s %7.3f, %7.3f, %7.4f | '%(name[i], ra[i], dec[i], flux[i]),
            print '%7.3f, %7.3f, %7.4f | '%(ra_bin[ra_map_index[i]], 
                                            dec_bin[dec_map_index[i]], 
                                            flux_in_hipass[i]),
            print '%7.3f, %7.3f, %7.4f | '%(ra_bin[ra_map_index[i]], 
                                            dec_bin[dec_map_index[i]], 
                                            flux_in_map[i]),
            if rescale:
                print '%7.3f, %7.3f, %7.4f | '%(ra_bin[ra_map_index[i]], 
                                                dec_bin[dec_map_index[i]], 
                                                flux_in_map_rescaled[i])
            else:
                print
        print '--------------------------------------------------------'
        print 'In the map, we got %d sources;'%flux_in_map.shape[0]
        print 'In NVSS,    we got %d sources.'%flux.shape[0]
        print 'In HIPASS,  we got %d sources.'%flux_in_hipass.shape[0]
        if rescale:
            print '[%7.5f] x flux + [%7.5f]'%(self.re_scale, self.re_zerop)
        print '========================================================'
        print

        #plt.figure(figsize=self.plot_size)
        fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=self.plot_size)
        plt.subplots_adjust(left=0.06, right=0.95, bottom=0.15, top=0.95, wspace=0.001)

        #plt.subplot(121)
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

        #plt.subplot(122)
        for i in range(len(ra_map_index)):
            if rescale:
                ax[1].plot(spec_axis/1.e6, spec_in_map_rescaled[:, i])
            else:
                ax[1].plot(spec_axis/1.e6, spec_in_map[:, i])
        #ax[1].plot(spec_axis/1.e6, np.mean(np.mean(self.imap, axis=-1), axis=-1), 
        #           c='0.5', linewidth=2)
        ax[1].set_xlabel('Frequency [MHz]')
        #ax[1].set_ylabel('Flux [Jy]')
        ax[1].set_xlim(xmin=spec_axis.min()/1.e6, xmax=spec_axis.max()/1.e6)

        if rescale:
            plt.savefig('./png/%s_check_nvss_rescaled.png'%(self.name), format='png')
        else:
            plt.savefig('./png/%s_check_nvss.png'%(self.name), format='png')
        #plt.show()

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

    def mapping_coord(self, plot=False):

        hipass_data, ra_bin_edges, dec_bin_edges, ra, dec\
            = self.hipass_cat(plot=False, cut=False)

        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')


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
                                            axis_names=self.imap.info['axes'])
        hipass_data_new.info = copy.deepcopy(self.imap.info)
        hipass_data_new.info['freq_centre'] = 1420.e6

        hipass_data_new = self.degrade_resolution(hipass_data_new, 
                               degrade_factor=self.degrade_factor)
        hipass_data_new = hipass_data_new[0,...]

        if plot:
            plt.figure(figsize=self.plot_size)
            max = 0.5*np.max(hipass_data_new)
            min = 0.5*np.min(hipass_data_new)
            RA, DEC = np.meshgrid(ra, dec)
            plt.pcolormesh(RA, DEC, hipass_data_new.T, vmax=max, vmin=min)
            plt.colorbar()

            nvss = self.nvss_cat()
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

        suffix = ''

        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        RA, DEC = np.meshgrid(ra, dec)
        #print RA.shape
        #print DEC.shape
        if num==None:
            #print "average map"
            plt.figure(figsize=self.plot_size)
            map = np.copy(self.imap)
            map = np.ma.array(map[self.freq_list, ...])
            #map = np.ma.array(self.imap)
            map[np.isnan(map)] = np.ma.masked
            map[np.isinf(map)] = np.ma.masked
            map = np.ma.mean(map, axis=0)
            #map = map[20, ...]
            map = map.T

            #if np.fabs(min) > np.fabs(max):
            #    max = np.fabs(min)
            #else:
            #    min = -np.fabs(max)
            #print map.shape

            if self.re_scale != None:
                suffix += '_rescaled'
                map *= self.re_scale
                map += self.re_zerop

            if diff:
                suffix += '_diff'
                hipass = self.mapping_coord()
                map -= hipass.T

            max = np.max(map) * 0.5
            min = np.min(map) * 0.5

            plt.pcolormesh(RA, DEC, map, vmax=max, vmin=min)
            plt.colorbar()


            if with_nvss:
                nvss = self.nvss_cat()
                nvss['flux'] *= 1.e-3
                selected = nvss['flux'] > .5
                ra  = nvss['ra'][selected]
                dec = nvss['dec'][selected]
                flux= nvss['flux'][selected]

                cmap = plt.get_cmap('jet')
                norm = plt.normalize(min, max)
                if diff:
                    plt.scatter(ra, dec, s=40, c='none', edgecolor='k' )
                else:
                    plt.scatter(ra, dec, s=40, c=cmap(norm(flux)), edgecolor='k' )

            plt.axes().set_aspect('equal')
            plt.xlim(RA.min(), RA.max())
            plt.ylim(DEC.min(), DEC.max())
            plt.xlabel('RA [deg]')
            plt.ylabel('Dec[deg]')
            plt.savefig('./png/%s_average%s.png'%(self.name, suffix), format='png')
            #plt.show()
        else:
            print "map freq %d"%num
            i = num
            plt.figure(figsize=self.plot_size)
            map = np.copy(self.imap)
            map = np.ma.array(map[self.freq_list, ...])
            #map = np.ma.array(self.imap)
            map[np.isnan(map)] = np.ma.masked
            map[np.isinf(map)] = np.ma.masked

            #map = np.ma.mean(map, axis=0)
            #map = map[i, ...] - np.ma.mean(map, axis=0)
            map = map[i, ...]

            map = map.T
            max = np.ma.max(map) * 0.5
            min = np.ma.min(map) * 0.5
            plt.pcolormesh(RA, DEC, map, vmax=max, vmin=min)
            plt.colorbar()
            plt.xlim(RA.min(), RA.max())
            plt.ylim(DEC.min(), DEC.max())
            plt.xlabel('RA [deg]')
            plt.ylabel('Dec[deg]')
            plt.axis('equal')
            #plt.savefig('./png/%s_%d.png'%(self.name, i), format='png')
            #plt.savefig('./png/%s_%d-mean.png'%(self.name, i), format='png')
            plt.show()

        

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

    #maproot = "/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/"
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

    mapdict['imap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_clean_map_I_1315.npy'
    mapdict['nmap'] = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
    mapdict['name'] = 'test_allbeams_27n30_10by7_clean_map_I'
    #mapdict['name'] = 'test_allbeams_27n30_10by7_clean_map_I_1pt1_cov'
    mapdict['imap_sec'] = ['/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_A_clean_map_I_1315.npy',
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_B_clean_map_I_1315.npy', 
                           '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/test_allbeams_27n30_10by7_C_clean_map_I_1315.npy', 
                          ]

    #mapdict['imap'] = '/home/ycli/data/map_result/maps/parkes/fir_p2750n2950_parkes_2010_10_24-28_clean_map_I_1316.npy'
    #mapdict['nmap'] = '/home/ycli/data/map_result/maps/parkes/fir_p2750n2950_parkes_2010_10_24-28_noise_weight_I_1316.npy'
    #mapdict['name'] = 'p2750n2950_parkes_2010_10_24-28_clean_map_I_1316_1pt1_cov'

    #mapdict['imap'] = '/home/ycli/data/map_result/maps/parkes/fir_p3500n3000_parkes_2010_10_24-28_dirty_map_I_1316.npy'
    #mapdict['nmap'] = '/home/ycli/data/map_result/maps/parkes/fir_p3500n3000_parkes_2010_10_24-28_dirty_map_I_1316.npy'
    #mapdict['name'] = 'p3500n3000_parkes_2010_10_24-28_clean_map_I_1316_1pt1_cov'

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

    #mapdict['imap'] = '/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_24_clean_map_I_1316.npy'
    #mapdict['nmap'] = '/home/ycli/data/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/fir_p3500n3000_parkes_2010_10_24_noise_weight_I_1316.npy'
    #mapdict['name'] = 'fir_p3500n3000_parkes_2010_10_24_clean_map'

    #degrade_factor = 1.1
    degrade_factor = 0
    if degrade_factor != 0:
        mapdict['name'] += '_1pt1_cov'

    a = plotmap(mapdict, freq_cut = [0,1,2,3,4,5,59,60,61,62,63], 
                degrade_factor=degrade_factor, plot_size=(15, 6))
    a.mapping_coord(plot=True)
    # #a.plot_map(num=10)
    a.hipass_cat()
    a.plot_map(with_nvss=True)
    a.plot_map(with_nvss=True, diff=True)
    a.check_nvss(flux_limit=[0.5, 100], rescale=False)
    # #a.rescale_by_hipass()
    a.check_nvss(flux_limit=[0.5, 100])
    a.plot_map(with_nvss=True)
    a.plot_map(with_nvss=True, diff=True)
    #a.plot_map()
    #a.make_mov()

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


