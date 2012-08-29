import os
import sys
import numpy as np
import ephem
from core import fitsGBT
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import combine_cal
# rebin_time?

def load_moonscan(filename):
    cal_coords = ephem.Equatorial("05:42:36.155", "+49:51:07.28",
                                  epoch=ephem.B1950)

    # Convert cal to a Body object.
    cal_source = ephem.FixedBody()
    cal_source._ra = cal_coords.ra
    cal_source._dec = cal_coords.dec
    cal_source._epoch = cal_coords.epoch

    Reader = fitsGBT.Reader(filename)
    Data = Reader.read(0,0)

    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    rebin_freq.rebin(Data, 16, True, True)
    #rebin_time.rebin(Data, 4)

    az = Data.field['CRVAL2']
    el = Data.field['CRVAL3']
    date_time = Data.field['DATE-OBS']
    source_az = np.empty_like(az)
    source_el = np.empty_like(el)
    #for ii in range(Data.dims[0])

    print az

datapath = "/mnt/raid-project/gmrt/kiyo/data/guppi_data/GBT12A_418/"
filename = datapath + "22_MOON_track_10.fits"
load_moonscan(filename)
