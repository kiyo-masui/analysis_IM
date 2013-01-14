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

np.seterr(all='warn')

data_root = os.getenv('GBT_DATA') + 'GBT12A_418/'
end = '.fits'

# Basic parameters.
g_n_beam_comp = 4 # number of hermite polynomials to Fit to each polarization.

# Polinomials to take out of each scan as noise. 1 = constant, 2 = linear, etc.
g_n_scan_comp = 2 

#g_source = '3C147'
g_source = '3C295'

# These files we will use to calibrate.
# Slow scans.
#cal_files = ['22_3C147_track_' + str(ii) for ii in range(27, 35)]
cal_files = ['22_3C295_track_' + str(ii) for ii in range(59, 66)]


# Read and preprocess the Data.
cal_Blocks = []
g_n_time = 0
for fname in cal_files:
    # Read.
    fpath = data_root + fname + end
    Reader = fitsGBT.Reader(fpath)
    Data = Reader.read(0,0)
    cal_Blocks.append(Data)
    g_n_time += Data.dims[0]

for Data in cal_Blocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    #rebin_freq.rebin(Data, 16, True, True)
    rebin_freq.rebin(Data, 16, True, True)
    combine_cal.combine(Data, (0.5, 0.5), False, True)
    #rebin_time.rebin(Data, 4)

# Set some globals.
g_Blocks = cal_Blocks
g_n_blocks = len(cal_Blocks)
# Here we get an initial guess for the important parameters.
first_Data = cal_Blocks[0]
g_n_cal = first_Data.dims[2]
first_Data.calc_freq()
freq = first_Data.freq.copy()
n_f = len(freq)
# Initialize all parameters.
first_Data.calc_pointing()
#source = [np.mean(first_Data.ra), np.mean(first_Data.dec)]
source = [0, 0]
width = 0.3
# XX and YY system Temperatures.
XX_T = np.mean(first_Data.data[0,0,0,:])
YY_T = np.mean(first_Data.data[0,3,0,:])
# XX and YY peaks.
time_mid_ind = first_Data.dims[0] // 2
XX_peak = np.mean(first_Data.data[time_mid_ind,0,0,:]) - XX_T
YY_peak = np.mean(first_Data.data[time_mid_ind,3,0,:]) - YY_T
beam = np.zeros((g_n_beam_comp, g_n_beam_comp, 4))
beam[0,0,0] = XX_peak
beam[0,0,3] = YY_peak
# Scan means and slopes.
scan_params = np.zeros((g_n_blocks, 4, g_n_cal, g_n_scan_comp))
scan_params[:,0,:,0] = XX_T
scan_params[:,3,:,0] = YY_T
# Pack them up and initialize space for to store fits.
init_pars = pack_parameters(source, width, beam, scan_params)
zero_pars = pack_parameters(source, width, 
                            np.zeros((g_n_beam_comp, g_n_beam_comp, 4)),
                            scan_params)
all_pars = np.empty((n_f, len(init_pars)), dtype=np.float64)


ii = 30
g_chan_ind = ii
g_chan_data = preprocess_blocks(cal_Blocks, ii)
#get_residuals(init_pars)
fit_pars, ier = optimize.leastsq(get_residuals, init_pars, ftol=1e-4, 
                                 xtol=1e-4)

d, m, w = get_chan_fit_data(fit_pars, cal_Blocks, ii)
plt.figure()
plt.plot(d[0,0,:])
plt.plot(m[0,0,:])


for ii in range(n_f):
    print "Processing channel ", ii
    g_chan_ind = ii
    g_chan_data = preprocess_blocks(cal_Blocks, ii)
    #get_residuals(init_pars)
    fit_pars, ier = optimize.leastsq(get_residuals, init_pars, ftol=1e-4, 
                                     xtol=1e-4)
    all_pars[ii] = fit_pars


# Get the beam array.
n_side = 200
size = 2  # degrees
full_beam = np.empty((n_f, 4, n_side, n_side), dtype=np.float64)
for ii in range(n_f):
    full_beam[ii,...] = beam_map_from_params(size, n_side, all_pars[ii])

np.save("fit_pars_" + g_source + '.npy', all_pars)
np.save("full_beam" + g_source + '.npy', full_beam)
np.save("freq" + g_source + '.npy', freq)




# ==== Utilities ====

# Parameters are:
#    source location              2
#    width                        1
#    coefficients                 4 * n_modes * (n_modes + 1) / 2
#    mean/slope for each channel  n_scan_comp * 4 * 2 * n_blocks

def get_residuals(pars):
    data, model, weights = get_chan_fit_data(pars, g_Blocks, g_chan_ind)
    residuals = (data - model) * weights
    return residuals.flat[:]

def get_chan_fit_data(pars, Blocks, chan_ind):
    chan_data = preprocess_blocks(Blocks, chan_ind)
    source, width, beam, scans = unpack_parameters(pars)
    # Make the beam object.
    Beam = pol_beam.SimpleBeam()
    Beam.set_width(width)
    Beam.set_coefficients(beam)
    # Initialize the residuals array.
    data = np.empty((4, g_n_cal, g_n_time), dtype=np.float64)
    model =  np.zeros((4, g_n_cal, g_n_time), dtype=np.float64)
    weight =  np.empty((4, g_n_cal, g_n_time), dtype=np.float64)
    ra_source = source[0]
    dec_source = source[1]
    # Loop through the data calculate the residuals.
    time_ind = 0
    for jj in range(g_n_blocks):
        block_data = chan_data[jj]
        n_t = block_data['data'].shape[-1]
        data[:,:,time_ind:time_ind + n_t] = block_data['data']
        weight[:,:,time_ind:time_ind + n_t] = block_data['weight']
        # Add mean, slope and other scan components.
        time = block_data['time']
        polys = misc.ortho_poly(time, g_n_scan_comp)
        # Renormalize the 0th polynomial such that the amplitude is interpreted
        # as the system temperature.
        polys[0,:] = 1
        for kk in range(g_n_scan_comp):
            model[:,:,time_ind:time_ind + n_t] += (scans[jj,:,:,kk,None]
                                                   * polys[kk,:])
        # Calculate the beam model.
        ra = block_data['ra']
        dec = block_data['dec']
        # XXX
        # ra_factor = np.cos(dec_source * np.pi / 180)
        ra_factor = 1
        beam_model = Beam.get_slice((ra_source - ra) * ra_factor,
                               dec_source - dec)
        model[:,:,time_ind:time_ind + n_t] += beam_model[0,:,None,:]
        # Noise weight.
        time_ind += n_t
    return data, model, weight

def beam_map_from_params(size, n_side, pars):
    source, width, beam, scans = unpack_parameters(pars)
    Beam = pol_beam.SimpleBeam()
    Beam.set_width(width)
    Beam.set_coefficients(beam)
    return Beam.get_full(n_side, size)[0]


def plot_from_params(pars):
    source, width, beam, scans = unpack_parameters(pars)
    Beam = pol_beam.SimpleBeam()
    Beam.set_width(width)
    Beam.set_coefficients(beam)
    n_side = 200
    side = 3 * width
    beam_map = Beam.get_full(n_side, side)
    axis = np.arange(n_side, dtype=float) / n_side * side
    for ii in range(4):
        plt.figure()
        plt.imshow(beam_map[0,ii,:,:])
        #plt.axis('equal')
        plt.colorbar()


def beam_movie(beam, fname):    
    # Get dimensions.
    nx = beam.shape[2]
    ny = beam.shape[3]
    ns = beam.shape[1]
    nf = beam.shape[0]
    # Normalize the data.
    # Should I interpolate to the middle?
    norm_x = beam[:,0,nx//2,ny//2].copy()
    norm_y = beam[:,3,nx//2,ny//2].copy()
    beam = beam.copy()
    beam[:,0,:,:] /= norm_x[:,None,None]
    beam[:,3,:,:] /= norm_y[:,None,None]
    beam[:,[1,2],:,:] /= np.sqrt(norm_x * norm_y)[:,None,None,None]
    # Compress colors for more contrast.
    sign = beam < 0
    beam = abs(beam) ** (1./3)
    beam[sign] *= -1
    # Reshape the data, to display all 4 polarizations along the x axis.
    beam.shape = (nf, ns * nx, ny)
    # Set up plotting axes.
    fig = plt.figure(figsize=(5 * ns, 5))
    ax = plt.axes(xlim=(0, ns * nx), ylim=(0, ny))
    ax.set_aspect('equal')
    ax.imshow(beam[50,...].T, vmin=-1, vmax=1)
    # Define the plotting function for every frame
    def animate(ii):
        img = ax.imshow(beam[ii,...].T, vmin=-1, vmax=1)
    # Animate it.
    anim = animation.FuncAnimation(fig, animate, frames=nf, interval=100,
                                   blit=True)
    anim.save('basic_animation.mp4', fps=10)





def pack_parameters(source, width, beam, scans):
    # Check inputs.
    if len(source) != 2:
        raise ValueError('Source')
    width = float(width)
    if beam.shape == (1, g_n_beam_comp, g_n_beam_comp, 4):
        beam = beam.view()
        beam.shape = (g_n_beam_comp, g_n_beam_comp, 4)
    elif beam.shape != (g_n_beam_comp, g_n_beam_comp, 4):
        raise ValueError('beam')
    if scans.shape != (g_n_blocks, 4, g_n_cal, g_n_scan_comp):
        raise ValueError('scans')
    # Get the number and size of each class of parameters.
    pars_beam_start = 3
    n_pars_beam =  g_n_beam_comp * (g_n_beam_comp + 1) / 2 * 4
    pars_scans_start = pars_beam_start + n_pars_beam
    n_pars_scans = g_n_blocks * 4 * g_n_cal * g_n_scan_comp
    n_pars = 2 + 1 + n_pars_beam + n_pars_scans
    pars = np.empty(n_pars, dtype=np.float64)
    # Pack up the parameters.
    pars[0:2] = source
    pars[2] = width
    # Pack the beam compoentes.
    start = 3
    for ii in range(g_n_beam_comp):
        pars[start:start + (g_n_beam_comp - ii) * 4] = \
                beam[ii,:(g_n_beam_comp - ii),:].flat[:]
        start += (g_n_beam_comp - ii) * 4
    # Pack the scan nuance parameters.
    pars[pars_scans_start:] = scans.flat
    return pars

def unpack_parameters(pars):
    # Get the number and size of each class of parameters.
    pars_beam_start = 3
    n_pars_beam =  g_n_beam_comp * (g_n_beam_comp + 1) / 2 * 4
    pars_scans_start = pars_beam_start + n_pars_beam
    n_pars_scans = g_n_blocks * 4 * g_n_cal * g_n_scan_comp
    n_pars = 2 + 1 + n_pars_beam + n_pars_scans
    if pars.shape != (n_pars,):
        raise ValueError()
    # Pack up the parameters.
    source = pars[0:2]
    width = pars[2]
    # Pack the beam compoentes.
    beam = np.zeros((g_n_beam_comp, g_n_beam_comp, 4))
    start = 3
    for ii in range(g_n_beam_comp):
        beam[ii,:(g_n_beam_comp - ii),:].flat[:] = \
                pars[start:start + (g_n_beam_comp - ii) * 4]
        start += (g_n_beam_comp - ii) * 4
    beam.shape = (1, g_n_beam_comp, g_n_beam_comp, 4)
    # Pack the scan nuance parameters.
    scans = np.empty((g_n_blocks, 4, g_n_cal, g_n_scan_comp))
    scans.flat[:] = pars[pars_scans_start:]
    return source, width, beam, scans
 
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
        # In sky fram.
        #Data.calc_pointing()
        #ra = Data.ra
        #dec = Data.dec
        # Telescope frame.
        az = Data.field['CRVAL2']
        el = Data.field['CRVAL3']
        UT = Data.field['DATE-OBS']
        Source = cal.source.Source(g_source)
        az_s, el_s = Source.azelGBT(UT)
        ra_factor = np.cos(np.mean(el_s) * np.pi / 180)
        rel_az = (az - az_s) * ra_factor
        rel_el = el - el_s
        # Get the time array.
        Data.calc_time()
        time = Data.time
        # Store it all and add to list.
        this_data = {}
        this_data['data'] = data
        this_data['weight'] = weight
        #this_data['ra'] = ra
        #this_data['dec'] = dec
        this_data['ra'] = rel_az
        this_data['dec'] = rel_el
        this_data['time'] = time
        data_list.append(this_data)
    return data_list




# == Unused ==

# 3C147
cal_coords = ephem.Equatorial("05:42:36.155", "+49:51:07.28", 
                              epoch=ephem.B1950)







