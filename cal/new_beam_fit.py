import os

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import optimize
import ephem
from numpy import ma
import matplotlib.animation as animation

from core import fitsGBT, dir_data
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import rebin_time, combine_cal
from map import pol_beam
from utils import misc
import cal.source
from cal import beam_fit

#np.seterr(all='warn')

# GBT_DATA = '/mnt/raid-project/gmrt/kiyo/data/guppi_data/'
data_root = os.getenv('GBT_DATA') + 'GBT12A_418/'
end = '.fits'

#source = '3C147'
#source = '3C295'
source = '3C286'

# These files we will use to calibrate.
# Slow scans.
#cal_files = ['22_3C147_track_' + str(ii) for ii in range(27, 35)]
#cal_files = ['22_3C295_track_' + str(ii) for ii in range(59, 67)]
cal_files = ['21_3C286_track_' + str(ii) for ii in range(18, 26)]

# The following two loops is a standard set of things we do to our raw data
# when it comes from the telescope.  Our data format is roughly SDfits.
# Read and preprocess the Data.
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
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    #combine_cal.combine(Data, (0.5, 0.5), False, True)
    combine_cal.combine(Data, (0., 1.), False, True)
    #rebin_time.rebin(Data, 4)

Data.calc_freq()

# Put all the data into a same format for the fits.
BeamData = beam_fit.FormattedData(cal_Blocks)

# Source object.  This just calculates the ephemeris of the source compared to
# where the telescope is pointing.
S = cal.source.Source(source)

# Do a preliminary fit to just the XX and YY polarizations.  This is a
# non-linear fit to the Gaussian and gets things like the centriod and the
# Gaussian width.  All fits are channel-by-channel (independantly).
center_offset, width, amps, Tsys = beam_fit.fit_simple_gaussian(BeamData, S)

# Basis basis functions to be used in the fit.
HermiteBasis = pol_beam.HermiteBasis(Data.freq, center_offset, width)

# Perform the fit.
beam_params, scan_params, model_data = beam_fit.linear_fit(BeamData, HermiteBasis,
                                                           S, 3, 2)
# Make a beam object from the basis funtions and the fit parameters (basis
# coefficients).
Beam = pol_beam.LinearBeam(HermiteBasis, beam_params)

# Some plots.
beam_map = Beam.get_full_beam(100, 1)

plt.figure()
this_data, this_weight = BeamData.get_data_weight_chan(35)
plt.plot(this_data[0,:])
plt.plot(model_data[35,0,:])

pol_beam.plot_beam_map(beam_map[128,...], color_map=0.5, side=1.,
                       normalize='max03', rotate='XXYYtoIQ')
cbar = plt.colorbar()
#cbar.set_label(r"Square root intensity, $ \sgn(I) \sqrt(|I|) $")
cbar.set_label(r"Square root intensity, (${\rm{sgn}}(I)\,\sqrt{|I|}$)")
plt.xlabel(r"IQUV, Azimuth (degrees)")
plt.ylabel(r"Elevation (degrees)")



# Other usefull plots, if you have two beam fits you want to compare.
pol_beam.compare_beam_maps(beam_map_295[35], beam_map_147[35], 1)
