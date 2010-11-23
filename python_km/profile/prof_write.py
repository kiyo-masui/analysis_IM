"""Script to run profile rebin_freq module."""

import os
import cProfile

import numpy.ma as ma

from core import fitsGBT

data_file_name = (os.getenv('GBT10B_DATA') +
                  '06_wigglez1hr_azel_180-187.raw.acs.fits')
out_file_name = 'write.prof'

Reader = fitsGBT.Reader(data_file_name)
Blocks = Reader.read()

Writer = fitsGBT.Writer()
#Writer.add_data(Blocks)

cProfile.run("Writer.add_data(Blocks)", out_file_name)

#Writer.write('write_prof_temp.fits')
#os.remove('write_prof_temp.fits')

