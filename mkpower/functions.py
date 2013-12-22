"""
    Public Function are defined here


"""
# public module
import numpy as np
import scipy as sp
import os
import gc
import string
from scipy import integrate
from math import *
from map import physical_gridding as gridding
import copy
#import fftw3 as FFTW
#import pyfftw.FFTW as FFTW

# kiyo module
from core import algebra
from utils import distance
from utils import fftutil

# li module
#import MakePower

class BOX(object):
    def __init__(self, imap1, imap2, weight1, weight2):
        self.imap1 = imap1
        self.imap2 = imap2
        self.weight1 = weight1
        self.weight2 = weight2

    def mapping_to_xyz(self):
        #self.box_bin, self.boxunit = get_box_xyz(self.imap1, self.boxshape)

        #self.ibox1, self.nbox1 = get_box(self.box_bin, self.imap1, self.weight1)
        #self.ibox2, self.nbox2 = get_box(self.box_bin, self.imap2, self.weight2)

        self.flag_nan(target="map")
        gridding_method = gridding.physical_grid_largeangle
        #gridding_method = gridding.physical_grid

        # shift test 
        print "Shift 2 pixels for testing "
        self.imap2[:, 2:, :] = self.imap2[:, :-2, :]
        self.imap2[:, :2, :] = 0.
        self.weight2[:, 2:, :] = self.weight2[:, :-2, :]
        self.weight2[:, :2, :] = 0.

        self.ibox1, ibox1_info = gridding_method(self.imap1, refinement=1, order=1)
        self.ibox2, ibox2_info = gridding_method(self.imap2, refinement=1, order=1)
        self.nbox1, nbox1_info = gridding_method(self.weight1, refinement=1, order=1)
        self.nbox2, nbox2_info = gridding_method(self.weight2, refinement=1, order=1)

        self.ibox1 = algebra.make_vect(self.ibox1, ibox1_info['axes'])
        self.ibox2 = algebra.make_vect(self.ibox2, ibox2_info['axes'])
        self.nbox1 = algebra.make_vect(self.nbox1, nbox1_info['axes'])
        self.nbox2 = algebra.make_vect(self.nbox2, nbox2_info['axes'])

        self.ibox1.info = ibox1_info
        self.ibox2.info = ibox2_info
        self.nbox1.info = nbox1_info
        self.nbox2.info = nbox2_info

        self.boxshape = []
        self.boxunit = []
        for i in range(self.ibox1.ndim):
            self.boxshape.append(self.ibox1.shape[i])
            axis_name = self.ibox1.axes[i]
            axis = self.ibox1.get_axis(axis_name)
            delta_axis = np.fabs(axis[1] - axis[0])
            self.boxunit.append(delta_axis)
        #axis = self.ibox1.get_axis()
        #self.boxunit = [self.ibox1.info['freq_delta'],
        #                self.ibox1.info['ra_delta'], 
        #                self.ibox1.info['dec_delta']]

    def flag_nan(self, target="box"):
        if target == "box":
            self.ibox1[np.isnan(self.ibox1)] = 0.
            self.ibox2[np.isnan(self.ibox2)] = 0.
            self.nbox1[np.isnan(self.ibox1)] = 0.
            self.nbox1[np.isnan(self.nbox1)] = 0.
            self.nbox2[np.isnan(self.ibox2)] = 0.
            self.nbox2[np.isnan(self.nbox2)] = 0.
            self.ibox1[np.isinf(self.ibox1)] = 0.
            self.ibox2[np.isinf(self.ibox2)] = 0.
            self.nbox1[np.isinf(self.ibox1)] = 0.
            self.nbox1[np.isinf(self.nbox1)] = 0.
            self.nbox2[np.isinf(self.ibox2)] = 0.
            self.nbox2[np.isinf(self.nbox2)] = 0.
        elif target == "map":
            self.imap1[np.isnan(self.imap1)] = 0.
            self.imap2[np.isnan(self.imap2)] = 0.
            self.weight1[np.isnan(self.imap1)] = 0.
            self.weight1[np.isnan(self.weight1)] = 0.
            self.weight2[np.isnan(self.imap2)] = 0.
            self.weight2[np.isnan(self.weight2)] = 0.
            self.imap1[np.isinf(self.imap1)] = 0.
            self.imap2[np.isinf(self.imap2)] = 0.
            self.weight1[np.isinf(self.imap1)] = 0.
            self.weight1[np.isinf(self.weight1)] = 0.
            self.weight2[np.isinf(self.imap2)] = 0.
            self.weight2[np.isinf(self.weight2)] = 0.

            self.weight1[self.weight1 < 1.e-20] = 0.
            self.weight2[self.weight2 < 1.e-20] = 0.

    def subtract_mean(self):
        self.ibox1 = self.ibox1 - self.ibox1.flatten().mean()
        self.ibox2 = self.ibox2 - self.ibox2.flatten().mean()

    def subtract_weighted_mean(self):
        mean1 = np.sum(np.sum(self.ibox1, -1), -1)
        mean1_weight = np.sum(np.sum(self.nbox1, -1), -1)
        mean1_weight[mean1_weight==0] = np.inf
        mean1 /= mean1_weight
        #print "Max and Min of weighted mean1 : ", mean1.max(), mean1.min()
        self.ibox1 -= mean1[:, None, None]

        mean2 = np.sum(np.sum(self.ibox2, -1), -1)
        mean2_weight = np.sum(np.sum(self.nbox2, -1), -1)
        mean2_weight[mean2_weight==0] = np.inf
        mean2 /= mean2_weight
        self.ibox2 -= mean2[:, None, None]
        #print "Max and Min of weighted mean2 : ", mean1.max(), mean2.min()

    def get_overdensity(self, map, sel):

        sel[sel==0] = np.inf
        map = map/sel - 1.
        sel[sel == np.inf] = 0.

        return map
        

    def estimate_ps_3d(self, window="blackman"):

        #self.flag_nan(target="box")

        #window_function = fftutil.window_nd(self.nbox1.shape, name=window)
        #window_function = fftutil.window_nd((self.nbox1.shape[0],), name=window)
        window_func = getattr(np, window)
        window_function = window_func(self.nbox1.shape[0])
        #self.nbox1 *= window_function[:, None, None]
        #self.nbox2 *= window_function[:, None, None]

        self.ibox1 *= self.nbox1
        self.ibox2 *= self.nbox2


        #self.subtract_mean()
        #self.subtract_weighted_mean()

        normal = np.sum(self.nbox1 * self.nbox2)
        delta_v = np.prod(self.boxunit)

        #  iput_1 = np.zeros(self.boxshape, dtype=complex)
        #  oput_1 = np.zeros(self.boxshape, dtype=complex)
        #  iput_1.imag = 0.
        #  iput_1.real = self.ibox1
        #  #plan_1 = FFTW.Plan(iput_1, oput_1, direction='forward', flags=['measure'])
        #  #FFTW.execute(plan_1)
        #  oput_1 = np.fft.fftn(iput_1)

        #  iput_2 = np.zeros(self.boxshape, dtype=complex)
        #  oput_2 = np.zeros(self.boxshape, dtype=complex)
        #  iput_2.imag = 0.
        #  iput_2.real = self.ibox2
        #  #plan_2 = FFTW.Plan(iput_2, oput_2, direction='forward', flags=['measure'])
        #  #FFTW.execute(plan_2)
        #  oput_2 = np.fft.fftn(iput_2)

        #  oput_1 = np.fft.fftshift(oput_1)
        #  oput_2 = np.fft.fftshift(oput_2)

        oput_1 = np.fft.fftshift(np.fft.fftn(self.ibox1))
        oput_2 = np.fft.fftshift(np.fft.fftn(self.ibox2))

        self.ps_3d  = (oput_1 * oput_2.conj()).real
        self.ps_3d *= delta_v
        self.ps_3d /= normal

        del oput_1
        del oput_2
        gc.collect()

    def get_k_bin_centre(self):
        #print self.boxunit
        k_bin_x = 2. * np.pi * np.fft.fftshift(np.fft.fftfreq(self.boxshape[0],
                                                              self.boxunit[0]))
        k_bin_y = 2. * np.pi * np.fft.fftshift(np.fft.fftfreq(self.boxshape[1],
                                                              self.boxunit[1]))
        k_bin_z = 2. * np.pi * np.fft.fftshift(np.fft.fftfreq(self.boxshape[2],
                                                              self.boxunit[2]))
        return k_bin_x, k_bin_y, k_bin_z

    def convert_ps_to_unitless(self):

        print "convert power to unitless", 
        print self.boxshape
        k_bin_x, k_bin_y, k_bin_z = self.get_k_bin_centre()

        k_bin_r = np.sqrt( (k_bin_x**2)[:, None, None] + 
                           (k_bin_y**2)[None, :, None] + 
                           (k_bin_z**2)[None, None, :] )

        import scipy.special
        ndim = 3.
        factor = 2. * np.pi ** (ndim / 2.) / scipy.special.gamma(ndim / 2.)
        factor /= (2. * np.pi) ** ndim

        #self.ps_3d = self.ps_3d * k_bin_r**3 / 2. / np.pi**2
        self.ps_3d = self.ps_3d * k_bin_r ** ndim * factor

    def convert_3dps_to_1dps(self, k_edges):

        print "convert 3d power ot 1d power",
        print self.boxshape
        k_bin_x, k_bin_y, k_bin_z = self.get_k_bin_centre()

        k_bin_r = np.sqrt( (k_bin_x**2)[:, None, None] + 
                           (k_bin_y**2)[None, :, None] + 
                           (k_bin_z**2)[None, None, :] )

        ps_3d_flatten = copy.deepcopy(self.ps_3d.flatten())
        k_bin_r = k_bin_r.flatten()[np.isfinite(ps_3d_flatten)]
        ps_3d_flatten = ps_3d_flatten[np.isfinite(ps_3d_flatten)]

        kn_1d, edges = np.histogram(k_bin_r, k_edges)
        ps_1d, edges = np.histogram(k_bin_r, k_edges, weights=ps_3d_flatten)

        kn_1d = kn_1d.astype(float)
        #kn_1d[kn_1d==0] = np.inf
        ps_1d[kn_1d != 0] /= kn_1d[kn_1d != 0] 
        ps_1d[kn_1d == 0] = 0.
        #kn_1d[kn_1d==np.inf] = 0.

        self.kn_1d = kn_1d
        self.ps_1d = ps_1d

    def convert_3dps_to_2dps(self, k_edges_p, k_edges_v):
        
        print "convert 3d power ot 2d power",
        print self.boxshape
        k_bin_x, k_bin_y, k_bin_z = self.get_k_bin_centre()

        k_bin_p = np.sqrt( k_bin_x**2)
        k_bin_v = np.sqrt( (k_bin_y**2)[:,None] + (k_bin_z**2)[None,:] )

        k_bin_2d = np.zeros(shape=(2,)+ k_bin_p.shape + k_bin_v.shape)
        k_bin_2d[0] = k_bin_p[:, None, None]
        k_bin_2d[1] = k_bin_v[None, :, :]

        ps_3d_flatten = copy.deepcopy(self.ps_3d.flatten())
        k_bin_2d = k_bin_2d.reshape(2, -1)
        k_bin_2d = k_bin_2d[:, np.isfinite(ps_3d_flatten)]
        ps_3d_flatten = ps_3d_flatten[np.isfinite(ps_3d_flatten)]

        #kn_2d, xedges, yedges = np.histogram2d(k_bin_2d[0], k_bin_2d[1],
        #                                       (k_edges_p, k_edges_v),)

        #ps_2d, xedges, yedges = np.histogram2d(k_bin_2d[0], k_bin_2d[1],
        #                                       (k_edges_p, k_edges_v),
        #                                       weights=ps_3d_flatten)
        kn_2d, xedges, yedges = np.histogram2d(k_bin_2d[1], k_bin_2d[0],
                                               (k_edges_p, k_edges_v),)

        ps_2d, xedges, yedges = np.histogram2d(k_bin_2d[1], k_bin_2d[0],
                                               (k_edges_p, k_edges_v),
                                               weights=ps_3d_flatten)

        kn_2d = kn_2d.astype(float)
        #kn_2d[kn_2d==0] = np.inf
        ps_2d[kn_2d!=0] /= kn_2d[kn_2d!=0]
        ps_2d[kn_2d==0] = 0.
        #kn_2d[kn_2d==np.inf] = 0. 

        self.kn_2d = kn_2d
        self.ps_2d = ps_2d



pi = np.pi
deg2rad = pi/180.
#----------------------------------------------------------------------#
def xyzv(ra, de, r, ra0=0.):
    x = r*sin(0.5*pi-de)*cos(ra-ra0)
    y = r*sin(0.5*pi-de)*sin(ra-ra0)
    z = r*cos(0.5*pi-de)
    v = r**2*sin(0.5*pi-de)
    return x, y, z, v

def fq2r(freq, freq0=1.42e9, Omegam=0.27, Omegal=0.73, h_0=0.72):
    """change the freq to distence"""
    zz =  freq0/freq - 1.
    cosmo = {'omega_M_0' : Omegam, 'omega_lambda_0' : Omegal, 'h' : h_0}
    cosmo = distance.set_omega_k_0(cosmo)
    zz = h_0*distance.comoving_distance_transverse(zz, **cosmo)
    return zz

def discrete(params, array):
    """discrete the data pixel into small size"""
    newarray = sp.zeros(params['discrete']*(array.shape[0]-1)+1)
    for i in range(0, array.shape[0]-1):
        delta = (array[i+1]-array[i])/float(params['discrete'])
        for j in range(0, params['discrete']):
            newarray[i*params['discrete']+j] = array[i] + j*delta
    newarray[-1] = array[-1]
    return newarray
#----------------------------------------------------------------------#
def get_radecr_bin_edges(map_temp):
    # find the ra bin edges
    ra   = map_temp.get_axis('ra')
    ra   = ra - map_temp.info['ra_centre']
    ra   = np.append(ra, ra[-1] + map_temp.info['ra_delta'])
    ra   = (ra - 0.5 * map_temp.info['ra_delta']) * deg2rad

    # find the dec bin edges
    dec  = map_temp.get_axis('dec')
    dec  = dec - map_temp.info['dec_centre']
    dec  = np.append(dec, dec[-1] + map_temp.info['dec_delta'])
    dec  = (dec - 0.5 * map_temp.info['dec_delta'])*deg2rad

    # find the r bin edges
    freq = map_temp.get_axis('freq')
    freq = np.append(freq, freq[-1] + map_temp.info['freq_delta'])
    freq = freq - 0.5 * map_temp.info['freq_delta']
    r    = fq2r(freq)

    return ra, dec, r

def getedge(map_temp):
    r'''
        This function is used to find the suitable box edges in cartesion coordinate 
        system of comoving distance for a map in ra-dec coordinate. 

        input:
            map_temp: one map as temp
        output:
            x_range, y_range, z_range
    '''
    ra, dec, r = get_radecr_bin_edges(map_temp)

    radec_map = np.zeros((3,) + r.shape + ra.shape + dec.shape)
    radec_map[0,...] = r[:, None, None]
    radec_map[1,...] = ra[None, :, None]
    radec_map[2,...] = dec[None, None, :]
    xyz_map = np.zeros((3,) + r.shape + ra.shape + dec.shape)
    xyz_map[0,...]=radec_map[0,...]*np.cos(radec_map[2,...])*np.cos(radec_map[1,...])
    xyz_map[1,...]=radec_map[0,...]*np.cos(radec_map[2,...])*np.sin(radec_map[1,...])
    xyz_map[2,...]=radec_map[0,...]*np.sin(radec_map[2,...])

    x = xyz_map[0].flatten()
    y = xyz_map[1].flatten()
    z = xyz_map[2].flatten()

    x_range = (x.min(), x.max())
    y_range = (y.min(), y.max())
    z_range = (z.min(), z.max())

    print x_range, y_range, z_range

    return x_range, y_range, z_range

#----------------------------------------------------------------------#
def get_box(box_bin, imap, nmap, mmap=None):

    (x, y, z) = box_bin
    xyz_box = np.zeros((3,) + x.shape + y.shape + z.shape)
    xyz_box[0,...] = x[:, None, None]
    xyz_box[1,...] = y[None, :, None]
    xyz_box[2,...] = z[None, None, :]
    xyz_box = xyz_box.reshape(-1, x.shape[0] * y.shape[0] * z.shape[0])

    radec_box = np.zeros((3, x.shape[0] * y.shape[0] * z.shape[0]))
    radec_box[0] = np.sqrt(xyz_box[0]**2 + xyz_box[1]**2 + xyz_box[2]**2)
    radec_box[1] = np.arctan2(-xyz_box[1], xyz_box[0])
    radec_box[2] = np.arcsin(xyz_box[2]/radec_box[0])

    ra, dec, r = get_radecr_bin_edges(imap)

    index_box = np.zeros((3, x.shape[0] * y.shape[0] * z.shape[0]), dtype=int)
    index_box[0] = np.digitize(radec_box[0], r)   - 1
    index_box[1] = np.digitize(radec_box[1], ra)  - 1
    index_box[2] = np.digitize(radec_box[2], dec) - 1

    del radec_box
    del xyz_box

    mask_box  = (index_box[0] == -1)
    mask_box |= (index_box[1] == -1)
    mask_box |= (index_box[2] == -1)
    mask_box |= (index_box[0] == imap.shape[0])
    mask_box |= (index_box[1] == imap.shape[1])
    mask_box |= (index_box[2] == imap.shape[2])

    index_box[0,mask_box] = imap.shape[0]
    index_box[1,mask_box] = 0
    index_box[2,mask_box] = 0

    ext_map = np.zeros(imap.shape[1:])
    print imap.dtype
    imap = np.concatenate((imap, ext_map[None,:,:]), axis=0)
    nmap = np.concatenate((nmap, ext_map[None,:,:]), axis=0)

    ibox = imap[index_box[0], index_box[1], index_box[2]]
    nbox = nmap[index_box[0], index_box[1], index_box[2]]

    ibox = ibox.reshape(x.shape + y.shape + z.shape)
    nbox = nbox.reshape(x.shape + y.shape + z.shape)

    del index_box
    gc.collect()

    return ibox, nbox

def get_box_xyz(map_temp, box_shape, position='centre'):
    '''
        position = 'centre' : return the value of bin centre
        position = 'edges' : return the value of bin edges, with shape + 1
    '''
    x_range, y_range, z_range = getedge(map_temp)
    boxunit = max((x_range[1]-x_range[0])/box_shape[0],
                  (y_range[1]-y_range[0])/box_shape[1],
                  (z_range[1]-z_range[0])/box_shape[2])
    #print boxunit
    x_centre = 0.5 * ( x_range[1] - x_range[0] ) + x_range[0]
    y_centre = 0.5 * ( y_range[1] - y_range[0] ) + y_range[0]
    z_centre = 0.5 * ( z_range[1] - z_range[0] ) + z_range[0]

    if position == 'centre':
        x_centre_idx = 0.5 * float(box_shape[0] - 1)
        y_centre_idx = 0.5 * float(box_shape[1] - 1)
        z_centre_idx = 0.5 * float(box_shape[2] - 1)

        x_bin = (np.arange(box_shape[0]) - x_centre_idx) * boxunit + x_centre
        y_bin = (np.arange(box_shape[1]) - y_centre_idx) * boxunit + y_centre
        z_bin = (np.arange(box_shape[2]) - z_centre_idx) * boxunit + z_centre

        return (x_bin, y_bin, z_bin), boxunit

    elif position == 'edges':
        x_centre_idx = 0.5 * float(box_shape[0])
        y_centre_idx = 0.5 * float(box_shape[1])
        z_centre_idx = 0.5 * float(box_shape[2])

        x_bin = (np.arange(box_shape[0] + 1) - x_centre_idx) * boxunit + x_centre
        y_bin = (np.arange(box_shape[1] + 1) - y_centre_idx) * boxunit + y_centre
        z_bin = (np.arange(box_shape[2] + 1) - z_centre_idx) * boxunit + z_centre

        return x_bin, y_bin, z_bin, boxunit

#----------------------------------------------------------------------#
def get_mapdict(dir, selection=None):
    r"""
    Generate a map dict according to the map file in a dir
    """

    if dir==None:
        return None

    maplist = os.listdir(dir)
    mapdict = {}
    for map in maplist:
        if map.split('.')[-1]=='npy':
            mapsplit = map.split('.')[0].split('_')
            if mapsplit[0] == 'sec':
                #print map
                key1 = mapsplit[1] + '_with_' + mapsplit[7]
                if mapsplit[2] == 'modes':
                    key2 = mapsplit[2]
                else:
                    key2 = mapsplit[4]
                if key2 == 'inv' or key2 == 'diag':
                    key2 = mapsplit[3] + '_' + key2
                key3 = mapsplit[-1]

                mapdict['%s;%s;%s'%(key1, key2, key3)] = dir + map
            if mapsplit[0] == 'combined':
                key1 = mapsplit[2]
                key2 = mapsplit[3]

                mapdict['%s;%s'%(key1, key2)] = dir + map

            if mapsplit[0] == 'secA' or\
               mapsplit[0] == 'secB' or\
               mapsplit[0] == 'secC' or\
               mapsplit[0] == 'secD':
                key1 = mapsplit[0][-1]
                key2 = mapsplit[3] + '_' + mapsplit[4]
                #if key2=='noise_inv' and mapsplit[5] == 'diag':
                #    key2='noise_weight'

                mapdict['%s;%s'%(key1, key2)] = dir + map

            if mapsplit[0] == 'sim' and mapsplit[1] == selection:
                key1 = int(mapsplit[2])

                mapdict['%d'%key1] = dir + map

            if len(mapsplit)>=3:
                if mapsplit[2] == '2df' and mapsplit[0] == selection:
                    if selection == 'sele' and mapsplit[-1] == 'separable':
                        return dir + map
                    #if selection == 'real' and len(mapsplit) == 3:
                    #    return dir + map
                    #if selection == 'mock' and len(mapsplit) == 4:
                    #    key1 = int(mapsplit[3])
                    #    mapdict['%d'%key1] = dir + map
                    #
                    if len(mapsplit)>=4 and mapsplit[3] == 'delta':
                        if selection == 'real':
                            return dir + map
                        elif selection == 'mock':
                            key1 = int(mapsplit[4])
                            mapdict['%d'%key1] = dir + map

            if mapsplit[0][:3] == 'reg':
                if mapsplit[0].split('.')[0][5:] == selection:
                    return dir + map
                elif len(mapsplit[0].split('.')[0]) == 12 and\
                     mapsplit[0].split('.')[0][5:-3] == 'rand':
                    key1 = string.atoi(mapsplit[0].split('.')[0][-3:])
                    mapdict['%s'%key1] = dir + map

        if os.path.isfile(dir+map) and map.split('.')[-1]=='pkl':
            mapsplit = map.split('.')[0].split('_')
            if mapsplit[0] == 'SVD':
                #print map
                key1 = mapsplit[2] + '_with_' + mapsplit[4]
                key2 = mapsplit[0]
                #print key1, key2

                mapdict['%s;%s'%(key1, key2)] = dir + map

    maps = [mapdict.keys(), mapdict]
    return maps

def fill(params, imap, nmap, mmap=None):
    """
    Function that used to fill the fftbox with the intensity map
    
    params : the params dict for each module
    imap : the direction to the intensity maps
    nmap : the direction to the noise maps
    mmap : the direction to the mock maps

    It will return the fft box and nbox.
    box is for the intensity maps, while
    nbox is for the noise intensity maps.
    If the mmap!=None, it also return mbox
    which for the mock maps
    """
    #params = self.params
    
    mapshape = np.array(imap.shape)

    r  = fq2r(imap.get_axis('freq'))
    ra = imap.get_axis('ra')*deg2rad
    de = imap.get_axis('dec')*deg2rad
    ra0= ra[int(ra.shape[0]/2)]
    ra = ra - ra0
    dra= ra.ptp()/ra.shape[0]
    dde= de.ptp()/de.shape[0]


    #ra_far = 0
    #de_far = 0
    #if ra.min()*ra.max()>0:
    #   if fabs(ra.min())<fabs(ra.max()):
    #       ra_far = ra.min()
    #   else:
    #       ra_far = ra.max()
    #if de.min()*de.max()>0:
    #   if fabs(de.min())<fabs(de.max()):
    #       de_far = de.min()
    #   else:
    #       de_far = de.max()
    #point = []
    #for i in range(3):
    #   point.append([xyzv(ra.min(), de.min(), r.min())[i], 
   #                 xyzv(ra.max(), de.min(), r.min())[i],
   #                 xyzv(ra.min(), de.max(), r.min())[i],
   #                 xyzv(ra.max(), de.max(), r.min())[i],
   #                 xyzv(ra.min(), de.min(), r.max())[i],
   #                 xyzv(ra.max(), de.min(), r.max())[i],
   #                 xyzv(ra.min(), de.max(), r.max())[i],
   #                 xyzv(ra.max(), de.max(), r.max())[i],
   #                 xyzv(ra_far,   de_far,   r.max())[i],
    #                     ])
    #   point[i].sort()
    #print point
    #print params['boxshape']
    #print (point[0][-1]-point[0][0])/params['boxshape'][0]
    #print (point[1][-1]-point[1][0])/params['boxshape'][1]
    #print (point[2][-1]-point[2][0])/params['boxshape'][2]
    #return 0

    #print r.min(), r.max()
    #print xyz(ra.min(), de.min(), r.min())
    #print xyz(ra.max(), de.min(), r.min())
    #print xyz(ra.min(), de.max(), r.min())
    #print xyz(ra.max(), de.max(), r.min())
    #print xyz(ra.min(), de.min(), r.max())
    #print xyz(ra.max(), de.min(), r.max())
    #print xyz(ra.min(), de.max(), r.max())
    #print xyz(ra.max(), de.max(), r.max())

    ###return 0

    mapinf = [ra.min(), dra, de.min(), dde]
    mapinf = np.array(mapinf)

    box = algebra.info_array(sp.zeros(params['boxshape']))
    box.axes = ('x','y','z')
    box = algebra.make_vect(box)
    
    box_xrange = params['Xrange']
    box_yrange = params['Yrange']
    box_zrange = params['Zrange']
    box_unit = params['boxunit']
    box_disc = params['discrete']

    print box_xrange, box_yrange, box_zrange, box_unit, 

    box_x = np.arange(box_xrange[0], box_xrange[1], box_unit/box_disc)
    box_y = np.arange(box_yrange[0], box_yrange[1], box_unit/box_disc)
    box_z = np.arange(box_zrange[0], box_zrange[1], box_unit/box_disc)

    #print box_x.shape
    #print box_y.shape
    #print box_z.shape

    boxshape = np.array(box.shape)*int(box_disc)

    boxinf0 = [0, 0, 0]
    boxinf0 = np.array(boxinf0)
    boxinf1 = [boxshape[0], boxshape[1], boxshape[2]]
    boxinf1 = np.array(boxinf1)

    print "MapPrepare: Filling the FFT BOX"
    MakePower.Filling(
        imap, r, mapinf, box, boxinf0, boxinf1, box_x, box_y, box_z)


    nbox = algebra.info_array(sp.zeros(params['boxshape']))
    nbox.axes = ('x','y','z')
    nbox = algebra.make_vect(nbox)

    #nbox = algebra.info_array(sp.ones(params['boxshape']))
    #nbox.axes = ('x','y','z')
    #nbox = algebra.make_vect(nbox)

    #print boxinf1
    MakePower.Filling(
        nmap, r, mapinf, nbox, boxinf0, boxinf1, box_x, box_y, box_z)


    if mmap != None:
        mbox = algebra.info_array(sp.zeros(params['boxshape']))
        mbox.axes = ('x','y','z')
        mbox = algebra.make_vect(mbox)
        MakePower.Filling(
            mmap, r, mapinf, mbox, boxinf0, boxinf1, box_x, box_y, box_z)
        return box, nbox, mbox
    else:
        return box, nbox


#----------------------------------------------------------------------#
def getresultf(params):
    """
    get the resultf 
    """
    resultf = params['resultf']
    if resultf == '':                               
        resultf = 'testresult'
        #resultf = params['hr'][0]                  
        #if len(params['last']) != 0:
        #   resultf = resultf + params['last'][0]
        #resultf = resultf + '-' + params['hr'][1]
        #if len(params['last']) != 0:
        #   resultf = resultf + params['last'][1]
    return resultf

#----------------------------------------------------------------------#
def getmap_halfz(map, half):
    freq   = map.get_axis('freq')
    z      =  1.42e9/freq - 1.
    z_half = (z.max()+z.min())/2
    half_indx = 0
    for i in range(len(z)):
        if z[i]>z_half:
            half_indx = i - 1
            break
    if half_indx == -1:
        print 'Error: z wrong order'
        exit()
    if half == 'lower':
        map = map[:half_indx]
        map.set_axis_info('freq',freq[map.shape[0]//2],map.info['freq_delta'])
    elif half=='upper':
        map = map[half_indx:]
        map.set_axis_info('freq',freq[map.shape[0]//2+half_indx],map.info['freq_delta'])
    return map

def getmap(imap_fname, nmap_fname, mmap_fname=None, half=None):
    """
    get the matrix of intensity map and noise map
    """
    #in_root = params['input_root']

    imap = algebra.load(imap_fname)
    imap = algebra.make_vect(imap)

    if half!=None:
        imap = getmap_halfz(imap, half)
    #print "--The neam value for imap is:",imap.flatten().mean(),"--"
    #imap = imap - imap.flatten().mean()
    if imap.axes != ('freq', 'ra', 'dec') :
        raise ce.DataError('AXES ERROR!')

    try:
        nmap = algebra.load(nmap_fname)
        nmap = algebra.make_vect(nmap)

        if half!=None:
            nmap = getmap_halfz(nmap, half)

        bad = nmap<1.e-5*nmap.flatten().max()
        nmap[bad] = 0.
        non0 = nmap.nonzero()
        #imap[non0] = imap[non0]/nmap[non0]
    except IOError:
        print 'NO Noise File :: Set Noise to One'
        nmap = algebra.info_array(sp.ones(imap.shape))
        nmap.axes = imap.axes
        nmap = algebra.make_vect(nmap)
    nmap.info = imap.info
    if nmap.axes != ('freq', 'ra', 'dec') :
        raise ce.DataError('AXES ERROR!')

    if mmap_fname != None:
        try:
            mmap = algebra.load(mmap_fname)
            mmap = algebra.make_vect(mmap)
            if half!=None:
                mmap = getmap_halfz(mmap, half)
        except IOError:
            print 'NO Mock File :: Make it!'
            mmap = algebra.info_array(
                2.*np.random.rand(imap.shape[0],imap.shape[1], imap.shape[2])-0.5)
            mmap.axes = imap.axes
            mmap = algebra.make_vect(mmap)
        
        return imap, nmap, mmap
    else:
        return imap, nmap
    
#----------------------------------------------------------------------#
def grothfunction(z, omegam, omegal):
    func = lambda z, omegam, omegal: (1+z)*(E(z, omegam, omegal))**(-3)
    D, Derr = integrate.quad(func, z, np.inf, args=(omegam, omegal))

    D = D*(5/2)*omegam*E(z, omegam, omegal)

    return D
#----------------------------------------------------------------------#
def getpower_th(params):
    """
        Get the theoretial power spectrum result given by camb
        Return : Tb*Tb * G * P_th
    """
    PKcamb = sp.load(params['input_root']+'pk_camb.npy')
    #PKnonl = sp.load(params['input_root']+'nonlPK_'+resultf+'.npy')
    #PKnonl_k = sp.load(params['input_root']+'k_nonlPK_'+resultf+'.npy')
    OmegaHI= params['OmegaHI']
    Omegam = params['Omegam']
    OmegaL = params['OmegaL']
    z = params['z']

    def get_distortion(b):
        f = (Omegam*(1+z)**3)**0.55
        t = (1.+(2./3.)*(f/b)+(1./5.)*(f/b)**2)
        return t

    #b_opt = 1.2
    #t_opt = get_distortion(b_opt)

    b_opt = 1.2 
    t_opt = 1.
        
    if params['optical']:
        print '\tb=%f'%b_opt, '\tt=%f'%t_opt
        PKcamb[1] = PKcamb[1]*b_opt**2*t_opt**2
        PKcamb[1] = PKcamb[1]*(PKcamb[0]**3)/2./3.1415926/3.1415926
        return PKcamb

    #b_gbt = 1.7
    #t_gbt = get_distortion(b_gbt)

    b_gbt = 1.35
    t_gbt = 1

    # Get the Tb
    a3 = (1+z)**(-3)
    Tb = 0.3*(OmegaHI)*((Omegam + a3*OmegaL)/0.29)**(-0.5)*((1.+z)/2.5)**0.5
    if params['PKunit']=='mK':
        Tb = Tb/1.e-3

    # Get G
    #def getG(z):
    #   xx = ((1.0/Omegam)-1.0)/(1.0+z)**3
    #   num = 1.0 + 1.175*xx + 0.3046*xx**2 + 0.005335*xx**3
    #   den = 1.0 + 1.875*xx + 1.021 *xx**2 + 0.1530  *xx**3

    #   G = (1.0 + xx)**0.5/(1.0+z)*num/den
    #   return G
    #G = getG(z)
    
    ##print G**2*(1+z)**(-2)
    #D  = grothfunction(z, Omegam, OmegaL)
    #D0 = grothfunction(0, Omegam, OmegaL)
    #print (D/D0)**2

    #PKcamb[1] = PKcamb[1]*(D/D0)**2
    #PKcamb[1] = PKcamb[1]*(G**2*(1.+z)**(-2)*6**2)

    if params['cross']: 
        PKcamb[1] = PKcamb[1]*Tb*b_gbt*b_opt*sqrt(t_opt*t_gbt)
    else: 
        print '\tb=%f'%b_gbt, '\tt=%f'%t_gbt
        PKcamb[1] = PKcamb[1]*(Tb**2)*b_gbt**2*t_gbt**2
        #PKcamb[1] = PKcamb[1]*(Tb**2)

    PKcamb[1] = PKcamb[1]*(PKcamb[0]**3)/2./3.1415926/3.1415926

    return PKcamb
#----------------------------------------------------------------------#

if __name__=="__main__":
    
    map_root = '/Users/ycli/DATA/2df/map_2929.5/'
    map_file = 'real_map_2df.npy'
     
    #map_root = '/Users/ycli/DATA/maps/'
    #map_file = 'secA_15hr_41-80_avg_fdgp_new_clean_map_I_800.npy'

    map_temp = algebra.load(map_root + map_file)
    map_temp = algebra.make_vect(map_temp)

    #x, y, z = getedge(map_temp)
    #print x, y, z

    ibox_kiyo, ibox_kiyo_info = physical_gridding.physical_grid(map_temp, 
                                                                refinement=1, 
                                                                pad=5, 
                                                                order=2)
    ibox_kiyo = algebra.make_vect(ibox_kiyo, ibox_kiyo_info['axes'])
    ibox_kiyo.info = ibox_kiyo_info
    algebra.save(map_root + 'kiyo_fftbox_' + map_file, ibox_kiyo)

    #box_bin, boxunit = get_box_xyz(map_temp, (512,128,64))
    #print box_bin, boxunit
    #ibox, nbox = get_box(box_bin, map_temp, map_temp)

    #np.save(map_root + 'fftbox_' + map_file, ibox)

    #fftbox = BOX((256,128,64), map_temp, map_temp, map_temp, map_temp,)
    #fftbox.mapping_to_xyz()
    #fftbox.convert_3dps_to_2dps()

    #imaps_a = get_mapdict('/mnt/scratch-gl/ycli/cln_result/15hr_ABCD_legendre_modes_0gwj_14conv_new/simmapmode_simmap_beam_000_subreal')

