import os

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import optimize
import ephem
from numpy import ma

from core import fitsGBT, dir_data
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import rebin_time, combine_cal
from map import pol_beam
from utils import misc

np.seterr(all='raise')

data_root = os.getenv('GBT_DATA') + 'GBT12A_418/'
end = '.fits'

# Basic parameters.
n_beam_comp = 3 # number of hermite polynomials to Fit to each polarization.


# These files we will use to calibrate.
# Slow scans.
cal_files = ['22_3C147_track_' + str(ii) for ii in range(27, 35)]


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
    rebin_freq.rebin(Data, 16, True, True)
    #rebin_time.rebin(Data, 4)

n_blocks = len(cal_Blocks)
first_Data = cal_Blocks[0]
first_Data.calc_freq()
freq = first_Data.freq.copy()
n_f = len(freq)
# Initialize all parameters.
first_Data.calc_pointing()
source = [np.mean(first_Data.ra), np.mean(first_Data.dec)]
width = 0.3
# XX and YY peaks.
time_mid_ind = first_Data.dims[0] // 2
XX_peak = np.mean(first_Data.data[time_mid_ind,0,0,:])
YY_peak = np.mean(first_Data.data[time_mid_ind,3,0,:])
beam = np.zeros((n_beam_comp, n_beam_comp, 4))
beam[0,0,0] = XX_peak
beam[0,0,3] = XX_peak
# Scan means and slopes.
scan_params = np.zeros((n_blocks, 4, 2, 2))
# System temperature is about 10 * T_cal.
scan_params[:,0,:,0] = 10
scan_params[:,3,:,0] = 10
# Pack them up and initialize space for to store fits.
init_pars = pack_parameters(source, width, beam, scan_params)
zero_pars = pack_parameters(source, width, 
                            np.zeros((n_beam_comp, n_beam_comp, 4)),
                            scan_params)
all_pars = np.empty((n_f, len(init_pars)), dtype=np.float64)


for ii in range(n_f):
    print "Processing channel ", ii
    this_chan_data, this_n_time = preprocess_blocks(cal_Blocks, ii)
    get_residuals(init_pars)
    fit_pars, ier = optimize.leastsq(get_residuals, init_pars, ftol=1e-4, 
                                     xtol=1e-4)
    all_pars[ii] = fit_pars






# ==== Utilities ====

# Parameters are:
#    source location              2
#    width                        1
#    coefficients                 4 * n_modes * (n_modes + 1) / 2
#    mean/slope for each channel  2 * 4 * 2 * n_blocks

def plot_from_params(pars):
    source, width, beam, scans = unpack_parameters(pars)
    Beam = pol_beam.SimpleBeam()
    Beam.set_width(width)
    Beam.set_coefficients(beam)
    n_side = 200
    side = 2 * width
    beam_map = Beam.get_full(n_side, side)
    axis = np.arange(n_side, dtype=float) / n_side * side
    for ii in range(4):
        plt.figure()
        plt.imshow(beam_map[0,ii,:,:])
        #plt.axis('equal')
        plt.colorbar()

def pack_parameters(source, width, beam, scans):
    # Check inputs.
    if len(source) != 2:
        raise ValueError()
    width = float(width)
    if beam.shape != (n_beam_comp, n_beam_comp, 4):
        raise ValueError()
    if scans.shape != (n_blocks, 4, 2, 2):
        raise ValueError()
    # Get the number and size of each class of parameters.
    pars_beam_start = 3
    n_pars_beam =  n_beam_comp * (n_beam_comp + 1) / 2 * 4
    pars_scans_start = pars_beam_start + n_pars_beam
    n_pars_scans = n_blocks * 4 * 2 * 2
    n_pars = 2 + 1 + n_pars_beam + n_pars_scans
    pars = np.empty(n_pars, dtype=np.float64)
    # Pack up the parameters.
    pars[0:2] = source
    pars[2] = width
    # Pack the beam compoentes.
    start = 3
    for ii in range(n_beam_comp):
        pars[start:start + (n_beam_comp - ii) * 4] = \
                beam[ii,-(n_beam_comp - ii):,:].flat[:]
        start += (n_beam_comp - ii) * 4
    # Pack the scan nuance parameters.
    pars[pars_scans_start:] = scans.flat
    return pars

def unpack_parameters(pars):
    # Get the number and size of each class of parameters.
    pars_beam_start = 3
    n_pars_beam =  n_beam_comp * (n_beam_comp + 1) / 2 * 4
    pars_scans_start = pars_beam_start + n_pars_beam
    n_pars_scans = n_blocks * 4 * 2 * 2
    n_pars = 2 + 1 + n_pars_beam + n_pars_scans
    if pars.shape != (n_pars,):
        raise ValueError()
    # Pack up the parameters.
    source = pars[0:2]
    width = pars[2]
    # Pack the beam compoentes.
    beam = np.zeros((n_beam_comp, n_beam_comp, 4))
    start = 3
    for ii in range(n_beam_comp):
        beam[ii,-(n_beam_comp - ii):,:].flat[:] = \
                pars[start:start + (n_beam_comp - ii) * 4]
        start += (n_beam_comp - ii) * 4
    beam.shape = (1, n_beam_comp, n_beam_comp, 4)
    # Pack the scan nuance parameters.
    scans = np.empty((n_blocks, 4, 2, 2))
    scans.flat[:] = pars[pars_scans_start:]
    return source, width, beam, scans
 
def get_residuals(pars):
    source, width, beam, scans = unpack_parameters(pars)
    # Make the beam object.
    Beam = pol_beam.SimpleBeam()
    Beam.set_width(width)
    Beam.set_coefficients(beam)
    # Initialize the residuals array
    residuals = np.empty((4, 2, this_n_time), dtype=np.float64)
    ra_source = source[0]
    dec_source = source[1]
    # Loop through the data calculate the residuals.
    time_ind = 0
    for jj in range(n_blocks):
        block_data = this_chan_data[jj]
        data = block_data['data']
        # Subtract off the mean.
        data = data - scans[jj,:,:,0,None]
        # Subtract off the slope.
        time = block_data['time']
        n_t = len(time)
        s = time - np.mean(time)
        s /= np.sum(s**2)
        data -= scans[jj,:,:,1,None] * s
        # Calculate the model.
        ra = block_data['ra']
        dec = block_data['dec']
        ra_factor = np.cos(dec_source * np.pi / 180)
        model = Beam.get_slice((ra_source - ra) * ra_factor,
                               dec_source - dec)
        data -= model[0,:,None,:]
        # Noise weight.
        data *= block_data['weight']
        residuals[:,:,time_ind:time_ind + n_t] = data
        time_ind += n_t
    return residuals.flat[:]

def preprocess_blocks(Blocks, ind):
    n = len(Blocks)
    data_list = []
    n_time = 0
    for Data in Blocks:
        n_time += Data.dims[0]
        # Get the data for this channel.
        data = np.rollaxis(Data.data[:,:,:,ind].filled(0), 0, 3)
        # Make the noise weights.
        mask = np.rollaxis(ma.getmaskarray(Data.data[:,:,:,ind]), 0, 3)
        weight = data.copy()
        weight[[1,2],:,:] = np.sqrt(weight[0,:,:] * weight[3,:,:])
        weight = 1/weight
        weight[mask] = 0
        # Get the pointing.
        Data.calc_pointing()
        ra = Data.ra
        dec = Data.dec
        # Get the time array.
        Data.calc_time()
        time = Data.time
        # Store it all and add to list.
        this_data = {}
        this_data['data'] = data
        this_data['weight'] = weight
        this_data['ra'] = ra
        this_data['dec'] = dec
        this_data['time'] = time
        data_list.append(this_data)
    return data_list, n_time



# ==== Unused ====

# Now calculate the pointing relative to the source.
# Convert cal to a Body object.
cal_source = ephem.FixedBody()
cal_source._ra = cal_coords.ra
cal_source._dec = cal_coords.dec
cal_source._epoch = cal_coords.epoch

for Data in cal_Blocks[:1]:
    az = Data.field['CRVAL2']
    el = Data.field['CRVAL3']
    date_time = Data.field['DATE-OBS']
    source_az = np.empty_like(az)
    source_el = np.empty_like(el)
    for ii in range(Data.dims[0]):


# 3C147
cal_coords = ephem.Equatorial("05:42:36.155", "+49:51:07.28", 
                              epoch=ephem.B1950)







