"""Script to run profile rebin_freq module."""

import os
import cProfile

import numpy.ma as ma

from core import fitsGBT
from time_stream import hanning

data_file_name = (os.getenv('GBT10B_DATA') +
                  '06_wigglez1hr_azel_180-187.raw.acs.fits')
out_file_name = 'hanning.prof'

Reader = fitsGBT.Reader(data_file_name)
Data = Reader.read(5, 3)
#Data.data = ma.copy(Data.data, 'F')

cProfile.run('hanning.hanning_smooth(Data)', out_file_name)

