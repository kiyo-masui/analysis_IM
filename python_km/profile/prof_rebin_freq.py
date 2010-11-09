"""Script to run profile rebin_freq module."""

import os
import cProfile

import numpy.ma as ma

from core import fitsGBT
from time_stream import rebin_freq

data_file_name = (os.getenv('GBT10B_DATA') +
                  '06_wigglez1hr_azel_180-187.raw.acs.fits')
out_file_name = 'rebin_freq.prof'

Reader = fitsGBT.Reader(data_file_name)
Data = Reader.read(5, 3)
#Data.data = ma.copy(Data.data, 'F')

cProfile.run('rebin_freq.rebin(Data, 2.0)', out_file_name)

