import core.algebra as al
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
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

base_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps/'

#for el in os.walk(base_dir):
#    print el[0]
#noise_files = glob.glob(base_dir + '*' + 'ra165dec0' + '*noise_inv*' + '.npy')
#map_files = glob.glob(base_dir + '*' + 'ra165dec0' + '*clean*' + '.npy')
noise_files = glob.glob(base_dir + '*' + 'ra33decn30' + '*noise_inv*' + '.npy')
map_files = glob.glob(base_dir + '*' + 'ra33decn30' + '*clean*' + '.npy')
noise_files.sort()
map_files.sort()

print map_files
print len(map_files)

plot_dir = '/home/p/pen/andersoc/plots/'

#plot_avg(noise_files, 0, plot_dir + 'noise_diags') 
plot_avg(map_files, 0, plot_dir + 'freq_avg_maps_ra33decn30')
