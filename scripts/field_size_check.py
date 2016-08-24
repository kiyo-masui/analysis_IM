import os
import glob
from astropy.io import fits
import numpy as np

file_middles = []


data_dir = '/scratch/p/pen/andersoc/second_parkes_pipe/rebinned/'
dir_len = len(data_dir)
for el in os.walk(data_dir):
    dir_files = glob.glob(el[0] + '/*' + '2df1' + '*.fits')
    for file in dir_files:
        file_middles.append(file[dir_len:])

ra = 'CRVAL3'
dec = 'CRVAL4'
def max(file, field):
    hdulist = fits.open(file)
    hdu_data = hdulist[1]
    ra_vec = hdu_data.data[field]
    max_ra = np.max(ra_vec)
    hdulist.close()
    return max_ra

def min(file, field):
    hdulist = fits.open(file)
    hdu_data = hdulist[1]
    vec = hdu_data.data[field]
    min = np.min(vec)
    hdulist.close()
    return min

ra_max = 0
ra_min = 360
dec_min = 90
dec_max = -90
for string in file_middles:
    ra_max_temp = max(data_dir + string, ra)
    ra_min_temp = min(data_dir + string, ra)
    dec_max_temp = max(data_dir + string, dec)
    dec_min_temp = min(data_dir + string, dec)
    if ra_max_temp > ra_max:
        ra_max = ra_max_temp
    if dec_max_temp > dec_max:
        dec_max = dec_max_temp
    if ra_min_temp < ra_min:
        ra_min = ra_min_temp
    if dec_min_temp < dec_min:
        dec_min = dec_min_temp

print ra_max
print ra_min
print dec_max
print dec_min
