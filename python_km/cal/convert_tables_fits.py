"""A script to make the cal tables kevin sent me to be in a format something
like I will eventually use.
"""

import os

import scipy as sp
import pyfits

from core import fitsGBT, data_block
from time_stream import stitch_windows_crude

# A calibration file from which I will copy headers and the like.
data_file_name = (os.getenv('GBT10B_DATA') + 
                  '/05_3c286_onoff_161-164.raw.acs.fits')

# Where I saved Kevin's Tables.
kevin_table_dir = os.getenv('GBT10B_OUT') + '/kevin_cal/'
out_file = kevin_table_dir + 'cal_21.fits'
input_root = kevin_table_dir + '21_3c48_113-116_noisecal_'
Writer = fitsGBT.Writer()

Reader = fitsGBT.Reader(data_file_name)
Blocks = Reader.read(0, ())

required_fields = ('CRVAL1', 'CDELT1', 'CRPIX1', 'SCAN')
CalBlocks = ()

# Loop over IFs and read in the appropriate cal data.
for Data in Blocks :
    CalT = data_block.DataBlock()

    for field in required_fields :
        CalT.set_field(field, Data.field[field], Data.field_axes[field],
                       Data.field_formats[field])

    CalT.set_field('CAL', ['R'], ('cal',), '1A')
    CalT.set_field('CRVAL4', [-5, -6], ('pol',), '1I')
    CalT.set_field('LST', [3546.2], ('time',), '1D')

    freq_centre = int(round(Data.field['CRVAL1']/1.0e6)) # MHz

    indata_xx = sp.loadtxt(input_root + str(freq_centre) + '_xx.txt')
    indata_yy = sp.loadtxt(input_root + str(freq_centre) + '_yy.txt')

    cal_data = sp.concatenate((indata_xx[:,[1]], indata_yy[:,[1]]), axis=1)
    cal_data = sp.swapaxes(cal_data, 0, 1)
    CalT.set_data(cal_data[sp.newaxis,:,sp.newaxis,:])
    CalT.verify()

    CalBlocks = CalBlocks + (CalT,)

BigData = stitch_windows_crude.stitch(CalBlocks)

Writer.add_data(BigData)
Writer.write(out_file)
