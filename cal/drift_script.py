import os

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import optimize

from core import fitsGBT
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import rebin_time, combine_cal

fnames = ['02_3C295_track_104.fits', '02_3C295_track_48.fits']
#fnames = ['02_3C295_track_104.fits'] #, '02_3C295_track_48.fits']
root = os.getenv('GBT_DATA') + 'GBT12A_418/'

Blocks = []

for fname in fnames:
    # Read.
    fpath = root + fname
    Reader = fitsGBT.Reader(fpath)
    Data = Reader.read(0,0)
    Blocks.append(Data)

for Data in Blocks:
    # Preprocess.
    rotate_pol.rotate(Data, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(Data, True, False, False, False, True)
    flag_data.flag_data(Data, 5, 0.1, 40)
    rebin_freq.rebin(Data, 8, True, True)
    rebin_time.rebin(Data, 4)
    #rotate_pol.rotate(Data, (1, 2, 3, 4))

def model(n_time, centre, width, amp_xx, amp_yy, amp_xy, amp_yx,
          off_xx, off_yy, off_yx, off_xy, slope_xx, slope_yy, slope_xy,
          slope_yx, gain_xx, gain_yy, quad_xx, quad_yy):
    
    # Preliminaries.
    time = np.arange(n_time, dtype=float) - centre
    out = np.empty((4, 2, n_time), dtype=float)
    # Generate a unit gaussian.
    gauss = np.exp(- time**2 / (2 * width**2))
    # Generate the four time series.
    out[0,:,:] = amp_xx * gauss + off_xx + slope_xx * time
    out[3,:,:] = amp_yy * gauss + off_yy + slope_yy * time
    out[1,:,:] = amp_xy * gauss + off_xy + slope_xy * time
    out[2,:,:] = amp_yx * gauss + off_yx + slope_yx * time
    # Add the noise cal to XX, YY and XY only.
    out[[0,3,1],0,:] += 1.0
    # Adjust by the nonlinear gain.
    out[0,:,:] = (1.0 + gain_xx) * out[0,:,:] + quad_xx * out[0,:,:]**2
    out[3,:,:] = (1.0 + gain_yy) * out[3,:,:] + quad_yy * out[3,:,:]**2
    # Crude model for cross polarization non-linear gain.
    cross_gain = np.sqrt((1 + gain_xx + quad_xx * out[0,:,:])
                         * (1 + gain_yy + quad_yy * out[3,:,:]))
    out[1,:,:] = out[1,:,:] * cross_gain
    out[2,:,:] = out[2,:,:] * cross_gain
    return out

def chi_2(params, data):
    nt = data.shape[2]
    model_data = model(nt, params[0], params[1], params[2], params[3],
                       params[4], params[5], params[6], params[7],
                       params[8], params[9], params[10], params[11],
                       params[12], params[13], params[14], params[15],
                       params[16], params[17])
    mask = ma.getmaskarray(data)
    data = data.filled(0)
    model_data[mask] = 0
    noise = model_data.copy()
    noise[[1,2],:,:] = np.sqrt(noise[0,:,:] * noise[3,:,:])
    noise[mask] = 1000000
    return ((data - model_data)/noise).flat

initial_params = np.array([500., 100., 15., 15., 0., 0., 10., 10., 0., 0.,
                           0., 0., 0., 0., 0., 0., 0., 0.])

for Data in Blocks:
    nt = Data.data.shape[0]
    for ii in [30, 100]:
        this_data = Data.data[:,:,:,ii]
        this_data = np.rollaxis(this_data, 0, 3)
        this_chi_2 = lambda x: chi_2(x, this_data)
        fit, err = optimize.leastsq(this_chi_2, initial_params)
        this_model = model(*(nt,) + tuple(fit))
        print fit
        plt.figure()
        plt.plot(this_data[0,1,:])
        plt.plot(this_data[3,1,:])
        plt.plot(this_model[0,1,:])
        plt.plot(this_model[3,1,:])
        plt.figure()
        plt.plot(this_data[1,1,:])
        plt.plot(this_data[2,1,:])
        plt.plot(this_model[1,1,:])
        plt.plot(this_model[2,1,:])
        plt.figure()
        plt.plot(this_data[0,0,:] - this_data[0,1,:])
        plt.plot(this_data[3,0,:] - this_data[3,1,:])
        plt.plot(this_data[1,0,:] - this_data[1,1,:])
        plt.plot(this_model[0,0,:] - this_model[0,1,:])
        plt.plot(this_model[3,0,:] - this_model[3,1,:])
        plt.plot(this_model[1,0,:] - this_model[1,1,:])


