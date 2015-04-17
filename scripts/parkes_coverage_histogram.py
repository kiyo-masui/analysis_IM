import fnmatch
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *

def find_pattern(pattern,root_dir):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

#matches = find_pattern("*west2*.sdfits","/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/")
#matches = open('sorted_datalist_2012_ycf.txt', 'r').read().splitlines()
#matches = open('final_datalist.txt', 'r').read().splitlines()
#matches = open('datalist_2008_center.txt', 'r').read().splitlines()
matches = []
data_dir = '/scratch/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS/'
for el in os.walk(data_dir):
    print el[0]
    dir_files = glob.glob(el[0] + '/*' + '2df1' + '*.fits')
    for file in dir_files:
        matches.append(file)

print len(matches)
print matches

def fix(list):
    wh_wrap = (list>180)
    list[wh_wrap] = list[wh_wrap] - 360
    return list

def ra_max(list):
    ra_max_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['RA'])
        ra_max_file = np.max(ra_vec)
        ra_max_list.append(ra_max_file)
    max_ra = max(ra_max_list)
    return max_ra

def ra_min(list):
    ra_min_list = []
    for file in list:
        #print file
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['RA'])
        ra_min_file = np.min(ra_vec)
        ra_min_list.append(ra_min_file)
    min_ra = min(ra_min_list)
    return min_ra

def dec_max(list):
    dec_max_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        dec_vec = hdu_data.data['DEC']
        dec_max_file = np.max(dec_vec)
        dec_max_list.append(dec_max_file)
    max_dec = max(dec_max_list)
    return max_dec

def dec_min(list):
    dec_min_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        dec_vec = hdu_data.data['DEC']
        dec_min_file = np.min(dec_vec)
        dec_min_list.append(dec_min_file)
    min_dec = min(dec_min_list)
    return min_dec

def saveplot(hist, extent, j):
    currentplot=plt.imshow(hist, vmax=45, vmin=0, extent=extent, aspect=1, interpolation='nearest')
    plt.colorbar()
    #plt.savefig('/cita/h/home-2/anderson/anderson/parkes_analysis_IM/parkes_2012_movie_yc/' + '{0:03}'.format(j), bbox_inches=0)
    plt.savefig('/scratch/p/pen/andersoc/plots/' + '{0:03}'.format(j), bbox_inches=0)
    plt.close()

#hitmap = np.zeros((50,200))
hitmap = np.zeros((1000,100))
ran=ra_min(matches)
rax=ra_max(matches)
decn=dec_min(matches)
decx=dec_max(matches)

print ran
print rax
print decn
print decx
#ran=38.5
#rax=41.5
#decn=-31.5
#decx=-28.5

#ran=24.5
#rax=27.5
#decn=-31.5
#decx=-28.5


i=1
#saveplot(hitmap, [-20,20,-24,-34],1)
for file in matches:
    #print file
    hdulist = pyfits.open(file)
    hdu_data = hdulist[1]
    ra_vec = fix(hdu_data.data['RA'])
    dec_vec = hdu_data.data['DEC']
    (skycov, dec_edge, ra_edge) = np.histogram2d(dec_vec, ra_vec,
                                                 range=[[-35,-25],[-40,60]],                                                  bins=hitmap.shape)
    hitmap += skycov
    extent=[ra_edge[0], ra_edge[-1], dec_edge[0], dec_edge[-1]]
    #extent=[22,32,-26.5,-33.5]
    #saveplot(hitmap[17:183,(3*125):(3*144)], extent, i)
    i += 1

'''for j in np.arange(0,200):
    if j+1>=67 and j+1<=183:
        hitmap[18,j]=1000
        hitmap[184,j]=1000
    if j+1>=18 and j+1<= 184:
        hitmap[j,67]=1000
        hitmap[j,86]=1000
        hitmap[j,105]=1000
        hitmap[j,125]=1000
        hitmap[j,144]=1000
        hitmap[j,164]=1000
        hitmap[j,183]=1000'''
saveplot(hitmap, extent, 1)

#extent=[ra_edge[0], ra_edge[-1], dec_edge[0], dec_edge[-1]]
#plt.imshow(hitmap, vmax=50, vmin=0, extent=extent, interpolation='nearest')
#plt.imshow(hitmap, extent=extent, interpolation='nearest')
#plt.colorbar()
#plt.show()
#print hitmap
