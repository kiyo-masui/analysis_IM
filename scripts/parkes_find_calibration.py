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

#matches = find_pattern("**.sdfits","/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/")
matches = open('sorted_datalist_2012_yc.txt', 'r').read().splitlines()

def fix(list):
    wh_wrap = (list>180)
    list[wh_wrap] = list[wh_wrap] - 360
    return list

def ra_diff_min(list, beam):
    ra_diff_list = []
    i=1;
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        beam_mask = hdu_data.data['beam'] == beam
        ra_vec = fix(hdu_data.data['CRVAL3'])[beam_mask]
        ra_diff_file = np.max(ra_vec) - np.min(ra_vec)
        ra_diff_list.append(ra_diff_file)
    min_diff = min(ra_diff_list)
    indices = ra_diff_list == min_diff
    for file in indices:
        if file == 1:
            print i
        i = i + 1
    return {'mindiff':min_diff, 'indices':indices}

def dec_diff_min(list, beam):
    ra_diff_list = []
    i=1;
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        beam_mask = hdu_data.data['beam'] == beam
        ra_vec = hdu_data.data['CRVAL4'][beam_mask]
        ra_diff_file = np.max(ra_vec) - np.min(ra_vec)
        ra_diff_list.append(ra_diff_file)
    min_diff = min(ra_diff_list)
    #print min_diff.argmin()
    indices = ra_diff_list == min_diff
    for file in indices:
        if file == 1:
            print i
        i = i + 1
    return {'mindiff':min_diff, 'indices':indices}

def dec_diff_max(list, beam):
    ra_diff_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        beam_mask = hdu_data.data['beam'] == beam
        ra_vec = hdu_data.data['CRVAL4'][beam_mask]
        ra_diff_file = np.max(ra_vec) - np.min(ra_vec)
        ra_diff_list.append(ra_diff_file)
    max_diff = max(ra_diff_list)
    #print max_diff.argmax()
    indices = ra_diff_list == max_diff
    return {'mindiff':max_diff, 'indices':indices}


def ra_diff_max(list, beam):
    ra_diff_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        beam_mask = hdu_data.data['beam'] == beam
        ra_vec = fix(hdu_data.data['CRVAL3'])[beam_mask]
        ra_diff_file = np.max(ra_vec) - np.min(ra_vec)
        ra_diff_list.append(ra_diff_file)
    max_diff = max(ra_diff_list)
    indices = ra_diff_list == max_diff
    return {'mindiff':max_diff, 'indices':indices}

def ra_big(list):
    ra_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['CRVAL3'])
        ra =ra_vec[10]
        ra_list.append(ra)
        i=1
    for file in ra_list:
        if file >= 65:
            print i
        i = i + 1
    return ra_list

ra_big(matches)
#print ra_diff_max(matches, 1)['mindiff']
#print dec_diff_min(matches, 1)['mindiff']
#print dec_diff_max(matches, 1)['mindiff']
