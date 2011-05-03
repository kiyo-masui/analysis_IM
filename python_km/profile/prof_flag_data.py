"""Script to run profile rebin_freq module."""

import os
import cProfile

import numpy.ma as ma

from core import fitsGBT
import time_stream.flag_data as flag_data

data_file_name = (os.getenv('GBT10B_OUT') + 'guppi_data/' +
                  '44_wigglez15hrst_ralongmap_297-304.fits')
out_file_name = 'flag_data.prof'

Reader = fitsGBT.Reader(data_file_name)
Data = Reader.read(5, 0)
#Data.data = ma.copy(Data.data, 'F')

cProfile.run('flag_data.apply_cuts(Data, 5, 5)', out_file_name)

