"""Script to run profile rebin_freq module."""

import os
import cProfile

import numpy.ma as ma

from core import fitsGBT
import time_stream.flag_data as flag_data

data_file_name = (os.getenv('GBT10B_DATA') +
                  '06_wigglez1hr_azel_180-187.raw.acs.fits')
out_file_name = 'flag_data.prof'

Reader = fitsGBT.Reader(data_file_name)
Data = Reader.read(5, 3)
#Data.data = ma.copy(Data.data, 'F')

cProfile.run('flag_data.apply_cuts(Data)', out_file_name)

