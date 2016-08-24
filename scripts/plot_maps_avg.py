import core.algebra as al
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import h5py
from itertools import combinations
import glob
import os

def plot_avg(file_list, av_ax, filename):
    count = 0
    rows = float(len(file_list))
    two = float(2)
    rows = int(math.ceil(rows/two))
    fig, ax = plt.subplots(nrows=rows, ncols=2, figsize=(12,18))
    for file in file_list:
        data = al.load(file)
        data = np.transpose(np.mean(data, axis=av_ax))
        print data.shape
        print np.max(data)
        ra_step = data.info['ra_delta']
        dec_step = data.info['dec_delta']
        ra_range = [data.info['ra_centre'] + ra_step*data.shape[1]/2, data.info['ra_centre'] - ra_step*data.shape[1]/2]
        dec_range = [data.info['dec_centre'] + dec_step*data.shape[0]/2, data.info['dec_centre'] - dec_step*data.shape[0]/2]
        data = np.fliplr(data)
        data = np.flipud(data)
        ax[count%rows][count/rows].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]])
        count += 1
        print ra_range
        print dec_range
    plt.savefig(filename)
    plt.clf()

def plot_2df(file_list, noise_list, freq, fname):
    rows = int(len(file_list))
    fig, ax = plt.subplots(nrows=rows, ncols=1)
    fig.subplots_adjust(hspace = 0.02, wspace = 0.02)
    count = 0
    maps = []
    high = 0
    low = 0
    for files in zip(file_list, noise_list):
        print files[0]
        print files[1]
        data = al.load(files[0])
        noise = al.load(files[1])
        data = np.transpose(data[freq,:,:])
        noise = np.transpose(noise[freq,:,:])
        #Subtract noise weighted mean from slice
        data -= np.mean(np.multiply(data,noise))/(np.mean(noise)*1.0)
        maps.append(data)
        high = max(np.max(data),high)
        low = min(np.min(data), low)
        print data.shape
        print np.max(data)
    im = 0
    for data in maps:
        ra_step = data.info['ra_delta']
        dec_step = data.info['dec_delta']
        ra_range = [data.info['ra_centre'] + ra_step*data.shape[1]/2, data.info['ra_centre'] - ra_step*data.shape[1]/2]
        dec_range = [data.info['dec_centre'] + dec_step*data.shape[0]/2, data.info['dec_centre'] - dec_step*data.shape[0]/2]
        data = np.fliplr(data)
        data = np.flipud(data)
        im = ax[count].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]], vmax = high, vmin = low)
        ax[count].set_xlabel('Right Ascension')
        ax[count].set_ylabel('Declination')
        count += 1
        print ra_range
        print dec_range
    fig.subplots_adjust(right=0.7)
    cbar_ax = fig.add_axes([0.75, 0.3, 0.05, 0.4])
    cbar = fig.colorbar(im, cbar_ax)
    cbar.set_label('Temperature fluctuations in mK', rotation = 270)
    #fig.tight_layout()
    plt.savefig(fname)
    plt.clf()

def save_2df_map_avg(file_list, noise_list, path_out, fname):
    print len(file_list)
    print len(noise_list)
    count = 0
    map_sum = 0
    noise_sum = 0
    for files in zip(file_list, noise_list):
        print files[0]
        print files[1]
        data = al.load(files[0])
        noise = al.load(files[1])
        data = np.multiply(data,noise)
        if count == 0:
            map_sum = data
            noise_sum = noise
        else:
            map_sum += data
            noise_sum += noise
        count += 1
    map_sum = np.divide(map_sum, noise_sum)
    al.save(path_out + fname + 'map' + '.npy', map_sum)
    al.save(path_out + fname + 'weights' + '.npy', noise_sum)



if __name__ == '__main__':
    base_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/conv_fixed/renormalized/'
    noise_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/'
    output_dir = '/home/p/pen/andersoc/plots/'
    #Plot beam 1 xx maps for 2df1, at freq 32.  
    combine_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/conv_fixed/renormalized/beams_combined/'
    cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/beams_removed/flux_cal/proper_conv/XX/'

#for el in os.walk(base_dir):
#    print el[0]
#noise_files = glob.glob(base_dir + '*' + 'ra165dec0' + '*noise_inv*' + '.npy')
#map_files = glob.glob(base_dir + '*' + 'ra165dec0' + '*clean*' + '.npy')
#noise_files = glob.glob(base_dir + '*' + 'ra33decn30' + '*noise_inv*' + '.npy')
    #map_files = glob.glob(base_dir + '*' + 'ran18decn30' + '*beam*' + '*clean*'  + '*XX*' + '.npy')
    #map_files = glob.glob(combine_dir + 'ran18decn30*' + '*map*' + '.npy')
    #map_files = glob.glob(base_dir + '*' + 'ran18decn30' + '*beam*' + '*clean*'  + '*XX*' + '.npy')
    #noise_files = glob.glob(noise_dir + '*' + 'ran18decn30' + '*beam*' + '*noise_inv_diag*'  + '*XX*' + '.npy')
    #noise_files = glob.glob(combine_dir + 'ran18decn30*' + '*weights*' + '.npy')

    #map_files.append(glob.glob(combine_dir + 'ra33decn30*' + '*map*' + '.npy')[0])
    #noise_files.append(glob.glob(combine_dir + 'ra33decn30*' + '*weights*' + '.npy')[0])

    map_files = [cleaned_dir + 'ran18/' + 'combined_clean_map_25modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_map_25modes.npy']
    noise_files = [cleaned_dir + 'ran18/' + 'combined_clean_weight_25modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_weight_25modes.npy']

    #noise_files = []
    #for file in map_files:
    #    noise_files.append(glob.glob(noise_dir + '*' + 'ran18decn30*' + file[file.index('beam'):file.index('clean')] + '*noise_inv_diag*'  + '*XX*' + '.npy')[0])
    plot_2df(map_files, noise_files, 32, output_dir + '2df1_xx_f32_allbeamsbut10_calibrated25modes')
    #save_2df_map_avg(map_files, noise_files, combine_dir, 'ran18decn30_XX_')
#noise_files.sort()
#map_files.sort()

#print map_files
#print len(map_files)

#plot_dir = '/home/p/pen/andersoc/plots/'

#plot_avg(noise_files, 0, plot_dir + 'noise_diags') 
#plot_avg(map_files, 0, plot_dir + 'freq_avg_maps_ra33decn30')
