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
from cal import flux_diff_gain_gen_beamfit


data_root = 'path/to/data/'
end = '.fits'

source = '3C286'

beam_cal_files = [file_middles]

gain_cal_files = [file_middles]

beam_cal_Blocks = []
for fname in beam_cal_files:
    # Read.
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    Data = Reader.read(0,0)
    cal_Blocks.append(Data)

for Data in beam_cal_Blocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    #combine_cal.combine(Data, (0.5, 0.5), False, True)
    combine_cal.combine(Data, (0., 1.), False, True)
    #rebin_time.rebin(Data, 4)

gain_cal_OnBlocks = []
gain_cal_OffBlocks = []
for fname in gain_cal_files:
    # Read.
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    OnData = Reader.read(0,0)
    OffData = Reader.read(1,0)
    gain_cal_OnBlocks.append(OnData)
    gain_cal_OffBlocks.append(OffData)

for Data in gain_cal_OnBlocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    #combine_cal.combine(Data, (0.5, 0.5), False, True)
    combine_cal.combine(Data, (0., 1.), False, True)
    #rebin_time.rebin(Data, 4)

for Data in gain_cal_OffBlocks:
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

BeamCalData = beam_fit.FormattedData(beam_cal_Blocks)
#GainCalOnData = beam_fit.FormattedData(gain_cal_OnBlocks)
#GainCalOffData = beam_fit.FormattedData(gain_cal_OffBlocks)

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
diff_gain_cal_gen_beamfit.calcGain(gain_cal_OnBlocks,gain_cal_OffBlocks,len(gain_cal_files),len(Data.freq),source,width)

# Make a beam object from the basis funtions and the fit parameters (basis
# coefficients).
Beam = pol_beam.LinearBeam(HermiteBasis, beam_params)

# Some plots.
beam_map = Beam.get_full_beam(100, 1)

#plt.figure()
#this_data, this_weight = BeamData.get_data_weight_chan(35)
#plt.plot(this_data[0,:])
#plt.plot(model_data[35,0,:])

#pol_beam.plot_beam_map(beam_map[35,...])






