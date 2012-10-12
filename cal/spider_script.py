import os

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import optimize
import ephem

from core import fitsGBT, dir_data
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import rebin_time, combine_cal

data_root = os.getenv('GBT_DATA') + 'GBT12A_418/'
end = '.fits'

# These files we will use to calibrate.
# Slow scans.
cal_files = ['22_3C147_track_' + str(ii) for ii in range(27, 35)]
# 3C147
cal_coords = ephem.Equatorial("05:42:36.155", "+49:51:07.28", 
                              epoch=ephem.B1950)


# Convert cal to a Body object.
cal_source = ephem.FixedBody()
cal_source._ra = cal_coords.ra
cal_source._dec = cal_coords.dec
cal_source._epoch = cal_coords.epoch


cal_Blocks = []
for fname in cal_files:
    # Read.
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    Data = Reader.read(0,0)
    cal_Blocks.append(Data)

for Data in cal_Blocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    rebin_freq.rebin(Data, 16, True, True)
    #rebin_time.rebin(Data, 4)

# Now calculate the pointing relative to the source.
for Data in cal_Blocks[:1]:
    az = Data.field['CRVAL2']
    el = Data.field['CRVAL3']
    date_time = Data.field['DATE-OBS']
    source_az = np.empty_like(az)
    source_el = np.empty_like(el)
    for ii in range(Data.dims[0]):







