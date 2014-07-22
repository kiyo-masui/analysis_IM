import fnmatch
import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def find_pattern(pattern,root_dir):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

matches = find_pattern("*east*.sdfits","/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/")
#a=np.arange(len(matches))
#b=str(a.tolist()).strip()
#print b
#print "\n".join(matches)

def fix(list):
    wh_wrap = (list>180)
    list[wh_wrap] = list[wh_wrap] - 360
    return list

def ra_max(list):
    ra_max_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['CRVAL3'])
        ra_max_file = np.max(ra_vec)
        ra_max_list.append(ra_max_file)
    max_ra = max(ra_max_list)
    return max_ra

def ra_min(list):
    ra_min_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['CRVAL3'])
        ra_min_file = np.min(ra_vec)
        ra_min_list.append(ra_min_file)
    min_ra = min(ra_min_list)
    return min_ra

def dec_max(list):
    dec_max_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        dec_vec = hdu_data.data['CRVAL4']
        dec_max_file = np.max(dec_vec)
        dec_max_list.append(dec_max_file)
    max_dec = max(dec_max_list)
    return max_dec

def dec_min(list):
    dec_min_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        dec_vec = hdu_data.data['CRVAL4']
        dec_min_file = np.min(dec_vec)
        dec_min_list.append(dec_min_file)
    min_dec = min(dec_min_list)
    return min_dec

hitmap = np.zeros((50,200))
weight_map = hitmap
ran=ra_min(matches)
rax=ra_max(matches)
decn=dec_min(matches)
decx=dec_max(matches)
i=1
for file in matches:
    hdulist = pyfits.open(file)
    hdu_data = hdulist[1]
    ra_vec = fix(hdu_data.data['CRVAL3'])
    dec_vec = hdu_data.data['CRVAL4']
    #map_data = np.mean(hdu_data.data['data'], axis=1)
    map_data = hdu_data.data['data'][:,300]
    #print map_data.shape
    (skycov, dec_edge, ra_edge) = np.histogram2d(dec_vec, ra_vec,                                                                                 range=[[decn,decx],[ran,rax]],
                                                 bins=hitmap.shape, weights=None)
    (weighted, dec_edge, ra_edge) = np.histogram2d(dec_vec, ra_vec,
                                range=[[decn,decx],[ran,rax]],
                                bins=weight_map.shape,
                                normed=False,                                                                                    weights=map_data)
    hitmap = hitmap + skycov
    weight_map = weight_map + weighted
    #print hitmap
    #print skycov
    #print i
    #i=i+1
hitrap = (hitmap < 1)
hitmap[hitrap] = hitmap[hitrap]+1
map=weight_map/hitmap
extent=[ra_edge[0], ra_edge[-1], dec_edge[-1], dec_edge[0]]
plt.imshow(map, extent=extent, interpolation='nearest')
#print weight_map
#print hitmap
#print weight_map.shape
plt.colorbar()
plt.show()
