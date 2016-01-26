import core.algebra as al
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import interpolate
import scipy as sp

beams= range(1,14)

dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps/'

first = 'parkes_parallel_thread_ra33decn30_p08_568by106_beam'

last = '_clean_map_XX_1316.npy'

#last = '_noise_inv_diag_XX_1316.npy'

small_maps =[]

for beam in beams:
    map = al.load(dir+first+str(beam)+last)
    #map = map[:,310:323,32:45]
    map = np.average(map, 0)
    map = np.transpose(map)
    map = np.fliplr(map)
    map = np.flipud(map)
    small_maps.append(map)

def plot(maps, filename):
    fig, ax = plt.subplots(nrows=7, ncols=2, figsize=(12,18))
    count = 0
    for file in maps:
        ax[count%7][count//7].imshow(file)
        count += 1
    plt.savefig(filename)
    plt.clf()

save_dir = '/home/p/pen/andersoc/plots/'

#plot(small_maps, save_dir + 'parkes_point_source_weights')

hipass_dir = '/gss01/scratch2/p/pen/andersoc/workcopy/HIPASS/'

hipass_file = hipass_dir + 'CHIPASS_Equ.fits'

hdulist = fits.open(hipass_file)

data = hdulist[0].data

dra = 1.0/15.0

ddec = -dra

def plot_single(map, filename):
    fig = plt.imshow(map)
    plt.savefig(filename)
    plt.clf()

def plot_pair(maps, filename):
    fig, ax = plt.subplots(nrows=len(maps), ncols=1)
    count = 0
    for map in maps:
        ax[count].imshow(map)
        count += 1
    plt.savefig(filename)
    plt.clf()

all = np.array(data)
all = np.fliplr(all)
all = np.flipud(all)

#Order is ra, dec
origin = [0,26]
#ra start, finish, followed by dec start, finish
#range = [147.4,182.6,origin[1] - 4.94,origin[1] + 5.94]
range = [10,55,origin[1] + 26 ,origin[1] + 34]

pixel_range = [int(round(x*15)) for x in range]

print pixel_range

print all.shape


map = np.array(all[pixel_range[2]:pixel_range[3],pixel_range[0]:pixel_range[1]])

print map.shape

#test interpolation
x = np.arange(range[2], range[3], 1.0/15)
y = np.arange(range[0], range[1], 1.0/15)
map_func = interpolate.RectBivariateSpline(x, y, map)

new_x = np.arange(range[2], range[3], 2.0/25)
new_y  = np.arange(range[0], range[1], (2.0/25)/sp.cos(sp.pi/6))

interp_map = map_func(new_x, new_y)

map_pairs = [map, interp_map]

plot_pair(map_pairs, hipass_dir + 'redgridding_decn30ra30')

#plot_single(map, hipass_dir + 'decn26n34ra10to55')


#plot_single(all, hipass_dir + 'entire_map')

#plot_single(all[0:1400,0:2800], hipass_dir + 'test2')

#test = np.zeros((2,2))
#test[0,1]=5
#test[1,1]=20

#plot_single(test, hipass_dir + 'test1')
