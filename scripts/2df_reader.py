#Script for analyzing 2dF data.

from astropy.io import fits
import os
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import fitsGBT
import numpy.ma as ma
import datetime

sgp_dir = '/scratch2/p/pen/andersoc/2df_data/sgp_fits/'
ngp_dir = '/scratch2/p/pen/andersoc/2df_data/ngp_fits/'

def headers(dir):
    headers = []
    for el in os.walk(dir):
    dir_files = glob.glob(el[0] + '/*' + '.fits')
    for file in dir_files:
        hdu = fits.open(file.split('.')[0][dir_len:])
        headers.append(hdu[0].header)
    return headers
