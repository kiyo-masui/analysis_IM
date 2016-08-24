#Made to calculate zero lag correlation between HIPASS and Parkes maps.
#For calibration.

import numpy as np
from astropy.io import fits
import core.algebra as al
import scipy as sp
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

from kiyopy import parse_ini
from map import beam
import core.algebra as al
import os

#from scripts import parkes_checkscale as pc

map_ra = os.getenv('MAP_RA')
map_dec = os.getenv('MAP_DEC')
map_size = os.getenv('MAP_SIZE')

params_init = {
              'x_beams':[1,2,3,4,5,6,7,8,9,11,12,13],
              'y_beams':[1,2,3,4,5,6,7,8,9,10,12],
              'map_dir': '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/hitconv_sync07/',
              'output_dir': '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/renorm_hitconv_sync07/',
              'noise_dir': '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/hitconv_sync07/',
              'prefix': 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam',
              'map_middle': '_clean_map_bp_div_',
              'noise_middle': '_noise_inv_diag_bp_div_',
              'suffix': '_1316.npy',
              'hipass_file': '/scratch2/p/pen/andersoc/workcopy/HIPASS/CHIPASS_Equ.fits'
              }
prefix = 'pc_'

def spectral_factor(bw=64.0, f=1315.5, index = -0.7):
    #Assumes CHIPASS as survey 1, which has 64 MHz bandwidth 
    #centered at 1394.5 MHz.
    #Our Parkes observations are at 1315.5 MHz center, 64 MHz bandwidth.
    #Calculates ratio of our survey magnitude to CHIPASS, assuming
    #an intensity power law with given index
    i = index + 1.0
    fact = ((f+bw/2)**i - (f - bw/2)**i)/(1426.5**i - 1362.5**i)
    return fact

def degrade_hipass_func(factor, ratio):
    #ratio: largest feq beam is this factor times 14.4 arcmin
    #factor:  multiply beam size by this number
    orig = 14.4/60.
    beam_diff = sp.sqrt((factor*ratio*orig)**2 - orig**2)
    common_resolution = beam.GaussianBeam(beam_diff)
    return common_resolution

def degrade_parkes(factor):
    freq_data = sp.array([1250, 1275, 1300, 1325, 1350, 1430], dtype=float)
    beam_data = sp.array([14.4, 14.4, 14.4, 14.4, 14.4, 14.4])/60. 
    beam_data = beam_data*1420./freq_data
    freq_data *= 1.0e6
    beam_diff = sp.sqrt(max(factor*beam_data)**2-(beam_data)**2)
    common_resolution = beam.GaussianBeam(beam_diff, freq_data)
    return common_resolution

def degrade_noise(noise, common_resolution):
    noise[noise < 1.e-30] = 1.e-30
    noise = 1. / noise
    noise = common_resolution.apply(noise, mode='constant', cval=1.e30)
    noise = common_resolution.apply(noise, mode='constant', cval=1.e30)
    noise = 1. / noise
    noise[noise < 1.e-20] = 0.
    return noise

def add_noise_freq(noise):
    noise[noise < 1.e-30] = 1.e-30
    noise = 1. / noise
    noise = np.sum(noise, axis=0)
    noise = 1. / noise
    noise[noise < 1.e-20] = 0.
    return noise

def w_avg(map, weight):
    weighted = np.multiply(map, weight)
    answer = np.sum(weighted)/np.sum(weight)
    return answer
    
def corr_zero(map1, noise1, map2, noise2):
    #Calculate weighted correlation, zero lag.
    #Noise is actually weight/inverse noise
    ans = np.sum(map1*noise1*map2*noise2)/np.sum(noise1*noise2)
    return ans

#hipass_dir = '/gss01/scratch2/p/pen/andersoc/workcopy/HIPASS/'
#hipass_file = hipass_dir + 'CHIPASS_Equ.fits'
#hdulist = fits.open(hipass_file)
#data = hdulist[0].data

def get_data_from_fits(path):
    hdulist = fits.open(path)
    data = hdulist[0].data
    data = np.fliplr(data)
    data = np.flipud(data)
    return data

def get_map_chunk(data, in_range):
    #data is assumed to be CHIPASS fits data. Dim is [dec, ra]
    #Ra, dec origin of CHIPASS data.
    #Range is ra, dec range desired: [ramin,ramax,decmin,decmax]
    #Must check for negative RA values
    print in_range
    ra_min = in_range[0] % 360
    ra_max = in_range[1] % 360
    #Now, check for wraparound
    #wrap = ra_min>ra_max
    map = np.transpose(data) 
    origin = [0,26]
    range = [ra_min, ra_max, origin[1] - in_range[3], origin[1] - in_range[2]]
    print range
    ra_pix_diff = int(round(in_range[1]*15)) - int(round(in_range[0]*15))
    pixel_range = [int(round(x*15)) for x in range] 
    #Remove the ra=360 values from the map, to wrap over to zero
    map = map[0:-1,:]
    #Now, shift the index so that index 0 is the ra_min pixel
    shift = map.shape[0] - pixel_range[0]
    print shift
    map = np.roll(map, shift, axis=0)
    map = map[0:ra_pix_diff, pixel_range[2]:pixel_range[3]]
    #map = map[pixel_range[0]:pixel_range[1], pixel_range[2]:pixel_range[3]]
    #Reverse ra order and also dec order to match pipeline maps.
    map = np.fliplr(map)
    map = np.flipud(map) 
    # Now, make the map into a proper info_array map.
    dra = -1.0/15.0
    ddec = -dra
    ra_center_ind = ra_pix_diff//2
    dec_center_ind = (pixel_range[3]-pixel_range[2])//2
    map = al.make_vect(map[None,:],axis_names=('freq', 'ra', 'dec'))
    map.set_axis_info('ra', in_range[1] + ra_center_ind*dra, dra)
    map.set_axis_info('dec', in_range[2] + dec_center_ind*ddec, ddec)
    map.set_axis_info('freq', 1394000000, 64000000)
    return map

def map_regrid(map, new_center, ra_space, dec_space, ra_pixels, dec_pixels):
    #Map should be info_array from CHIPASS
    #Standard interpolation is cubic RectBivariateSpline
    s_map = map[0,:,:]
    x = np.arange(map.shape[1])
    y = np.arange(map.shape[2])
    map_func = interpolate.RectBivariateSpline(x, y, s_map)
    print map.info['ra_delta']
    ra_fact = ra_space/map.info['ra_delta']
    print ra_fact
    dec_fact = dec_space/map.info['dec_delta']
    print dec_fact
    ra_new = give_coords(new_center[0], ra_pixels, ra_space)
    ra_orig = give_coords(map.info['ra_centre'], map.shape[1], map.info['ra_delta'])
    #print ra_new
    #print ra_orig
    dec_new = give_coords(new_center[1], dec_pixels, dec_space)
    dec_orig = give_coords(map.info['dec_centre'], map.shape[2], map.info['dec_delta'])
    ra_offset = (-ra_orig[0] + ra_new[0])/map.info['ra_delta']
    #print ra_offset
    dec_offset = (dec_new[0] - dec_orig[0])/map.info['dec_delta']
    new_x = [ra_offset + ra_fact*num for num in xrange(ra_pixels)]
    #print new_x
    new_y = [dec_offset + dec_fact*num for num in xrange(dec_pixels)]
    print new_y[0]
    print new_y[-1]
    new_map = map_func(new_x, new_y)
    #Now, we just need to make the new_map into an info_array map.
    new_map = al.make_vect(new_map[None,:],axis_names=('freq', 'ra', 'dec'))
    new_map.set_axis_info('ra', new_center[0], ra_space)
    new_map.set_axis_info('dec', new_center[1], dec_space)
    new_map.set_axis_info('freq', 1394000000, 64000000)
    print new_map.info
    return new_map

def give_coords(center, size, delta):
    #Returns array of points, with size//2 indexed element equal to center.
    #Array has input size, and the difference between consecutive elements is delta
    start = center - size//2 * delta
    output = [start + delta*x for x in xrange(size)]
    return output

class MapList():
    #Class to hold beam maps, beam weights, and CHIPASS map

    def __init__(self, parameter_file=None, params_dict=None):
        self.params = params_init
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix) 
        params = self.params
        self.maps = {}
        parkes_res = degrade_parkes(1.4)
        #Loading parkes maps and degrading resolution.
        for beam in params['x_beams']:
            map = parkes_res.apply(al.make_vect(al.load(params['map_dir'] + params['prefix'] + str(beam) + params['map_middle'] + 'XX' + params['suffix'])))
            map = np.mean(map, axis=0)
            noise = add_noise_freq(degrade_noise(al.make_vect(al.load(params['noise_dir'] + params['prefix'] + str(beam) + params['noise_middle'] + 'XX' + params['suffix'])), parkes_res))
            map -= w_avg(map, noise)
            self.maps['beam' + str(beam) + 'xx' + 'map'] = map
            self.maps['beam' + str(beam) + 'xx' + 'noise'] = noise
        for beam in params['y_beams']:
            map = parkes_res.apply(al.make_vect(al.load(params['map_dir'] + params['prefix'] + str(beam) + params['map_middle'] + 'YY' + params['suffix'])))
            map = np.mean(map, axis=0)
            noise = add_noise_freq(degrade_noise(al.make_vect(al.load(params['noise_dir'] + params['prefix'] + str(beam) + params['noise_middle'] + 'YY' + params['suffix'])), parkes_res))
            map -= w_avg(map, noise)
            self.maps['beam' + str(beam) + 'yy' + 'map'] = map
            self.maps['beam' + str(beam) + 'yy' + 'noise'] = noise
        #Now, get the appropriate section of the CHIPASS map.
        #Add 1 degree to edges for better interpolation
        temp = self.maps['beam' + str(params['x_beams'][0]) + 'xx' + 'map']
        info = temp.info
        ra_p = temp.shape[0]
        dra = info['ra_delta']
        ddec = info['dec_delta']
        self.ra_center = info['ra_centre']
        self.dec_center = info['dec_centre']
        self.dra = dra
        self.ddec = ddec
        dec_p = temp.shape[1]
        ra = give_coords(info['ra_centre'], ra_p, dra)
        dec = give_coords(info['dec_centre'], dec_p, ddec)
        range = [ra[-1] - 1, ra[0] + 1, dec[0] - 1, dec[-1] + 1]
        print range
        chi_map = get_map_chunk(get_data_from_fits(params['hipass_file']), range)
        chi_map = map_regrid(chi_map, [info['ra_centre'], info['dec_centre']], dra, ddec, ra_p, dec_p)
        #Subtract mean to avoid worst of annoying convolution edge effects.
        chi_map -= np.mean(chi_map)
        #Now, convolve to same resolution as other Parkes maps.  Assumes Parkes fmin = 1283.5 MHz
        hipass_beam = degrade_hipass_func(1.4, 1420./1283.5)
        chi_map = hipass_beam.apply(chi_map)
        self.maps['chipass'] = chi_map[0,:]
        self.calc_chipass_autos()
        self.calc_chipass_cross()
        self.auto_cross_ratio()
        self.tot_noise = np.zeros(self.maps['beam5xxnoise'].shape)
        self.w_map = np.zeros(self.maps['beam5xxmap'].shape)

    def graph_parkes_chipass(self, filename):
        fig, ax = plt.subplots(nrows=3, ncols=1)
        max = np.max(self.maps['chipass'])
        min = np.min(self.maps['chipass'])
        ax[0].imshow(np.transpose(self.maps['chipass']), vmin=min, vmax=max)
        ax[1].imshow(np.transpose(self.auto_cross_ratio['xx']['1']*self.maps['beam1xxmap']),vmin=min,vmax=max)
        ax[2].imshow(np.transpose(self.auto_cross_ratio['xx']['7']*self.maps['beam7xxmap']),vmin=min,vmax=max)
        plt.savefig(filename)
        plt.clf()

    def graph_w_map_chipass(self, filename):
        fig, ax = plt.subplots(nrows=3, ncols=1)
        ax[0].imshow(np.transpose(self.maps['chipass']))
        ax[1].imshow(np.transpose(self.w_map))
        ax[2].imshow(np.transpose(self.tot_noise))
        plt.savefig(filename)
        plt.clf()


    def noise_weight_sum(self, beams = 'x_beams', pol = 'xx'):
        tot_noise = self.tot_noise
        w_map = self.w_map
        for beam in self.params[beams]:
            noise = self.maps['beam' + str(beam) + pol + 'noise']
            tot_noise += noise
            map = self.maps['beam' + str(beam) + pol + 'map']
            w_map += noise*map
        w_map /= tot_noise
        self.w_map = w_map
        self.tot_noise = tot_noise

    def plot_maps(self, filename, beams = 'x_beams', factor = False, type = 'map'):
        #Choose whether to plot noise or map by changing type
        length = len(self.params[beams])
        if beams == 'x_beams':
            pol = 'xx'
        if beams == 'y_beams':
            pol = 'yy'
        fig, ax = plt.subplots(nrows = int(math.ceil(length/3.)), ncols=3, figsize = (15,15))
        largest = 0
        smallest = 0
        if factor:
            acr =  self.auto_cross_ratio
        for beam in self.params[beams]:
            if pol == 'xx' and beam == 11 and type == 'map':
                #This one seems way higher than the others
                largest = largest
            elif factor == False:
                largest = max(largest, np.max(self.maps['beam' + str(beam) + pol +type]))
                smallest = min(smallest, np.min(self.maps['beam' + str(beam) + pol +type]))
            else:
                largest = max(largest, np.max(acr[pol][str(beam)]*self.maps['beam' + str(beam) + pol +type]))
                smallest = min(smallest, np.min(acr[pol][str(beam)]*self.maps['beam' + str(beam) + pol +type]))
        i = 0
        for beam in self.params[beams]:
            if factor:
                mult = acr[pol][str(beam)]
            else:
                mult = 1
            map = np.transpose(mult*self.maps['beam' + str(beam) + pol +type])
            ra_step = self.dra
            dec_step = self.ddec
            ra_range = [self.ra_center - ra_step*map.shape[1]/2, self.ra_center + ra_step*map.shape[1]/2]
            dec_range = [self.dec_center + dec_step*map.shape[0]/2, self.dec_center - dec_step*map.shape[0]/2]
            ax[i//3][i%3].imshow(map, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]], vmin = smallest, vmax = largest)
            i += 1
        #fig.colorbar(ax[i-1])
        plt.savefig(filename)
        plt.clf()


    def save_calibrated_maps(self, index_fix = True):
        #Will save new maps and new noise inv diagonals.
        #If index_fix is true, it will attempt to account for CHIPASS maps
        #being at a higher frequency, assuming a synchrotron power law.
        #It also assumes most power comes from sources that smaller than a beam.
        if index_fix:
            c_fact = spectral_factor()
        else:
            c_fact = 1.0
        params = self.params
        for beam in params['x_beams']:
            map = al.load(params['map_dir'] + params['prefix'] + str(beam) + params['map_middle'] + 'XX' + params['suffix'])
            map *= (c_fact*self.auto_cross_ratio['xx'][str(beam)])/1000.0
            al.save(params['output_dir'] + params['prefix'] + str(beam) + params['map_middle'] + 'XX' + params['suffix'],map)
            noise = al.load(params['noise_dir'] + params['prefix'] + str(beam) + params['noise_middle'] + 'XX' + params['suffix'])
            noise /= ((c_fact*self.auto_cross_ratio['xx'][str(beam)])/1000.0)**2
            al.save(params['output_dir'] + params['prefix'] + str(beam) + params['noise_middle'] + 'XX' + params['suffix'],noise)
        for beam in params['y_beams']:
            map = al.load(params['map_dir'] + params['prefix'] + str(beam) + params['map_middle'] + 'YY' + params['suffix'])
            map *= (c_fact*self.auto_cross_ratio['yy'][str(beam)])/1000.0
            al.save(params['output_dir'] + params['prefix'] + str(beam) + params['map_middle'] + 'YY' + params['suffix'],map)
            noise = al.load(params['noise_dir'] + params['prefix'] + str(beam) + params['noise_middle'] + 'YY' + params['suffix'])
            noise /= ((c_fact*self.auto_cross_ratio['yy'][str(beam)])/1000.0)**2
            al.save(params['output_dir'] + params['prefix'] + str(beam) + params['noise_middle'] + 'YY' + params['suffix'],noise)

    def calc_chipass_autos(self):
        chipass_auto = {}
        chipass_auto['xx']={}
        chipass_auto['yy']={}
        for beam in self.params['x_beams']:
            noise = self.maps['beam' + str(beam) + 'xx' + 'noise']
            chi = self.maps['chipass']
            chi -= w_avg(chi, noise)
            auto = corr_zero(chi, noise, chi, noise)
            chipass_auto['xx']['noise' + str(beam)] = auto
        for beam in self.params['y_beams']:
            noise = self.maps['beam' + str(beam) + 'yy' + 'noise']
            chi = self.maps['chipass']
            chi -= w_avg(chi, noise)
            auto = corr_zero(chi, noise, chi, noise)
            chipass_auto['yy']['noise' + str(beam)] = auto
        self.chipass_auto = chipass_auto

    def calc_chipass_cross(self):
        chipass_cross = {}  
        chipass_cross['xx']={}
        chipass_cross['yy']={}
        for beam in self.params['x_beams']:
            noise = self.maps['beam' + str(beam) + 'xx' + 'noise']
            map = self.maps['beam' + str(beam) + 'xx' + 'map']
            chi = self.maps['chipass']
            chi -= w_avg(chi, noise)
            cross = corr_zero(chi, noise, map, noise)
            chipass_cross['xx']['noise' + str(beam)] = cross
        for beam in self.params['y_beams']:
            noise = self.maps['beam' + str(beam) + 'yy' + 'noise']
            map = self.maps['beam' + str(beam) + 'yy' + 'map']
            chi = self.maps['chipass']
            chi -= w_avg(chi, noise)
            cross = corr_zero(chi, noise, map, noise)
            chipass_cross['yy']['noise' + str(beam)] = cross
        self.chipass_cross = chipass_cross

    def auto_cross_ratio(self):
    #Only call after calculating auto and cross powers with chipass
        auto_cross_ratio = {}
        auto_cross_ratio['xx'] = {}
        auto_cross_ratio['yy'] = {}
        for beam in self.params['x_beams']:
            auto_cross_ratio['xx'][str(beam)] = self.chipass_auto['xx']['noise' + str(beam)]/self.chipass_cross['xx']['noise' + str(beam)]
        for beam in self.params['y_beams']:
            auto_cross_ratio['yy'][str(beam)] = self.chipass_auto['yy']['noise' + str(beam)]/self.chipass_cross['yy']['noise' + str(beam)]
        self.auto_cross_ratio = auto_cross_ratio

if __name__ == '__main__':
    maplist = MapList()
    maplist.save_calibrated_maps()
    #filename = maplist.params['output_dir'] + 'plot_' + map_ra 
    #maplist.graph_parkes_chipass(filename)

