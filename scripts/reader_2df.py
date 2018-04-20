#Script for analyzing 2dF data.

from astropy.io import fits
import os
import glob
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import fitsGBT
import numpy.ma as ma
import datetime
import core.algebra as al
import copy

sgp_dir = '/scratch2/p/pen/andersoc/2df_data/sgp_fits/'
ngp_dir = '/scratch2/p/pen/andersoc/2df_data/ngp_fits/'

def headers(dir):
    headers = []
    for el in os.walk(dir):
        dir_files = glob.glob(el[0] + '/*' + '.fits')
        for file in dir_files:
            hdu = fits.open(file)
            header0 = hdu[0].header
            header1 = hdu[1].header
            headers.append([hdu[0],hdu[1]])
            hdu.close()
    return headers

def make_data_array(header_list):
    #Make numpy array of all the data we want, for quicker access.
    #SEQNUM, in header 0, is the serial number.
    #SRRA, SRDEC are the RA and Dec, in header 0.
    #Z_HELIO, in header 1, is the heliocentric redshift.
    #SRRMAG is the UK-J (R) magnitude, in header 0.
    #SBBJMAG is the UK-J (Bj) magnitude, in header 0.
    #Data columns will be [serial #, ra, dec, redshift, r_mag, b_mag].
    data = np.zeros((len(header_list),6))
    for index in range(len(header_list)):
        header0 = header_list[index][0].header
        header1 = header_list[index][1].header
        data[index, 0] = header0['SEQNUM']
        data[index, 1] = header0['SRRA']
        data[index, 2] = header0['SRDEC']
        data[index, 3] = header1['Z_HELIO']
        data[index, 4] = header0['SRRMAG']
        data[index, 5] = header0['SBBJMAG']
    return data

def color_histogram(data, bins = 60, range = (0.,2.), apply_k_corrections = True):
    #Makes histogram of blue minus red magnitude
    br = data[:,5] - data[:,4]
    if apply_k_corrections:
        y = data[:,3]/(1 + data[:,3])
        k_bj = (-1.63 + 4.53*br)*y + (-4.03 - 2.01*br)*(y**2)
        k_r = (-0.08 + 1.45*br)*y + (-2.88 - 0.48*br)*(y**2)
        br = br - k_bj + k_r
    return np.histogram(br, bins = bins, range = range)

def z_color_histogram(data, br_bins = 60, z_bins = 60, apply_k_corrections = True, range = [[0.05,0.1],[0.5,1.5]]):
    br = data[:,5] - data[:,4]
    if apply_k_corrections:
        y = data[:,3]/(1 + data[:,3])
        k_bj = (-1.63 + 4.53*br)*y + (-4.03 - 2.01*br)*(y**2)
        k_r = (-0.08 + 1.45*br)*y + (-2.88 - 0.48*br)*(y**2)
        br = br - k_bj + k_r
    return np.histogram2d(data[:,3], br, bins = (z_bins, br_bins), range = range)

def save_2dhist_plot(hist, xedges, yedges, fname):
    plt.clf()
    im = plt.imshow(hist, interpolation='None', origin='low',
                extent=[yedges[0], yedges[-1],xedges[0], xedges[-1]], aspect = 'auto')
    plt.savefig(fname)
    plt.clf()

def save_hist_plot(hist_array, fname):
    plt.clf()
    width = 0.7 * (hist_array[1][1] - hist_array[1][0])
    center = (hist_array[1][:-1] + hist_array[1][1:]) / 2
    plt.bar(center, hist_array[0], align='center', width=width)
    plt.savefig(fname)

def open_cat(fname):
    #Returns [header string, data array].
    #Assumes header is two lines.
    f = open(fname)
    a = f.readline()
    a += f.readline()
    f.close()
    data = np.loadtxt(fname)
    ans = [a, data]
    return ans

def save_cat(header, array, fname):
    np.savetxt(fname, array, header=header)

def sub_cat(input_cat, serials):
    #Input catalog is np array, where last column is serial number.
    #Serials is 1d numpy array of serial numbers to be kept.
    #Output is np array in same format.
    bool = np.in1d(input_cat[:,-1], serials)
    output = input_cat[bool]
    return output

def make_shuffled_mocks(dir, num, cat_fname='catalog_2df', shuff_name='catalog_2df'):
    #cat_fname excludes suffix, assumed to be '.out' and prefix, assumed 'real_'
    #Makes n redshift shuffled versions of catalog stored at dir + cat_fname.
    #Outputs are saved in dir with name mock_cat_fname_n.out, where n is zero padded to 4 zeros.
    #Redshift assumed to be the 3rd column.
    cat = open_cat(dir + 'real_' + cat_fname + '.out')
    data = cat[1]
    for n in xrange(num):
        np.random.shuffle(data[:,2:data.shape[1]])
        save_cat(cat[0], data, dir + 'mock_' + cat_fname + '_' + str(n).zfill(4)+ '.out')

def make_red_blue_lists(cat, divider = 1.07, apply_k_corrections = True):
    #Blue minus red magnitude divider defaulted to 1.07.
    br = cat[:,5] - cat[:,4]
    if apply_k_corrections:
        y = cat[:,3]/(1 + cat[:,3])
        k_bj = (-1.63 + 4.53*br)*y + (-4.03 - 2.01*br)*(y**2)
        k_r = (-0.08 + 1.45*br)*y + (-2.88 - 0.48*br)*(y**2)
        br = br - k_bj + k_r
    red_bool = br >= divider

def thin_mock(map, ratio):
    #Map is an algebra array of mock counts.
    #Map is thinned via binomial choice of map[:] elements with probability equal to ratio.
    new = np.random.binomial(map.astype(int), ratio)
    ans = copy.deepcopy(map)
    ans[:] = new
    return ans

def sample_overdensity(over_dense, sel):
    #Each pixel will have an average number of selection*(1+overdensity).
    ans = copy.deepcopy(over_dense)
    ans[:] = np.random.poisson(sel*(1+over_dense))
    return ans

def calc_delta(map, sel):
    #map/sel - 1.
    bool = sel!=0
    ans = copy.deepcopy(map)
    ans[bool] = map[bool]/sel[bool] - 1.
    ans[sel==0] = -1.
    return ans

if __name__=="__main__":
    #Renormalize selection function for red and blue surveys, regular and sep.
    #Save separable selection function as both regular and separable
    ##Then, thin mock catalogs, calculate delta mocks (use regular sel).
    ##Calculate delta_real (regular selection function).
    #Actually, calculate mock catalogs using Nick's CLASS sim realizations.
    #Calculate delta_real and delta_mock (separable sel function).
    base_all = '/scratch2/p/pen/nluciw/parkes/maps/2df/'
    sim_base = '/scratch2/p/pen/nluciw/parkes/simulations/corr_effic/rawsim/'
    #base = '/scratch2/p/pen/andersoc/2df_data/catalog_color_maps/'
    base = '/scratch2/p/pen/andersoc/2df_data/'
    #full = 'full_sel_func/'
    full = 'catalog_color_maps/'
    blue = 'blue_k_corrected/'
    red = 'red_k_corrected/'
    all = 'catalog_maps/'
    fields = ['pks_n1800_2dfmock_full/', 'pks_p1990_2dfmock_full/', 'pks_p1650_2dfmock_full/',  'pks_p3300_2dfmock_full/']
    sim_fields = {'pks_n1800_2dfmock_full/':'ran18/', 'pks_p1990_2dfmock_full/':'ra199/', 'pks_p1650_2dfmock_full/':'ra165/',  'pks_p3300_2dfmock_full/':'ra33/'}
    for field in fields:
        #real_blue = al.load(base + blue + field + 'real_map_2df.npy')
        #al.save(base + full + blue + field + 'real_map_2df.npy', real_blue)
        #real_red = al.load(base + red + field + 'real_map_2df.npy')
        #al.save(base + full + red + field + 'real_map_2df.npy', real_red)
        real_blue = al.load(base + full + blue + field + 'real_map_2df.npy')
        real_red = al.load(base + full + red + field + 'real_map_2df.npy')
        real_all = al.load(base_all + field + 'real_map_2df.npy')
        #al.save(base + all + field + 'real_map_2df.npy', real_all)
        ratio_red = np.sum(real_red)/np.sum(real_all)
        ratio_blue = np.sum(real_blue)/np.sum(real_all)

        #sel = al.load(base_all + field + 'sele_map_2df.npy')
        sel = al.load(base_all + field + 'sele_map_2df_separable.npy')
        red_sel = sel*ratio_red
        blue_sel = sel*ratio_blue
        #al.save(base + full + blue + field + 'sele_map_2df.npy', blue_sel)
        #al.save(base + full + red + field + 'sele_map_2df.npy', red_sel)
        #al.save(base + all + field + 'sele_map_2df.npy', sel)

        sep = al.load(base_all + field + 'sele_map_2df_separable.npy')
        red_sep = sep*ratio_red
        blue_sep = sep*ratio_blue
        #al.save(base + all + field + 'sele_map_2df_separable.npy', sep)
        #al.save(base + full + blue + field + 'sele_map_2df_separable.npy', blue_sep)
        #al.save(base + full + red + field + 'sele_map_2df_separable.npy', red_sep)

        real_delta_blue = calc_delta(real_blue, blue_sel)
        #al.save(base + full + blue + field + 'real_map_2df_delta.npy', real_delta_blue)
        real_delta_red = calc_delta(real_red, red_sel)
        #al.save(base + full + red + field + 'real_map_2df_delta.npy', real_delta_red)
        real_delta_all = calc_delta(real_all, sel)
        #al.save(base + all + field + 'real_map_2df_delta.npy', real_delta_all)       
 
        ##Now cut down mocks and make mock_deltas
        #Now make mocks from Nick's sims and mock_deltas
        #for num in xrange(500):
        #for num in xrange(100):
        for num in xrange(100,200):
            mock = al.load(base_all + field + 'mock_map_2df_' + str(num).zfill(4) + '.npy')
            #sim = al.load(sim_base + sim_fields[field] + 'sim_raw_' + str(num).zfill(3) + '.npy')
            #red_mock = thin_mock(mock, ratio_red)
            #red_mock = sample_overdensity(sim, red_sel)
            #red_mock_delta = calc_delta(red_mock, red_sel)
            #al.save(base + full + red + field + 'mock_map_2df_' + str(num).zfill(4) + '.npy', red_mock)
            #al.save(base + full + red + field + 'mock_map_2df_delta_' + str(num).zfill(3) + '.npy', red_mock_delta)
            #blue_mock = thin_mock(mock, ratio_blue)
            #blue_mock = sample_overdensity(sim, blue_sel)
            #blue_mock_delta = calc_delta(blue_mock, blue_sel)
            #al.save(base + full + blue + field + 'mock_map_2df_' + str(num).zfill(4) + '.npy', blue_mock)
            #al.save(base + full + blue + field + 'mock_map_2df_delta_' + str(num).zfill(3) + '.npy', blue_mock_delta)

            #all_mock = sample_overdensity(sim, sel)
            all_mock = mock
            al.save(base + all + field + 'mock_map_2df_' + str(num).zfill(4) + '.npy', all_mock)
            all_mock_delta = calc_delta(all_mock, sel)
            al.save(base + all + field + 'mock_map_2df_delta_' + str(num).zfill(3) + '.npy', all_mock_delta)
