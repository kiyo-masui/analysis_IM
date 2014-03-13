import os
import math

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg

from core import fitsGBT, dir_data
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import rebin_time, combine_cal, reflag
from noise import noise_power

data_root = os.getenv('GBT_DATA') + 'GBT11B_055/'
end = '.fits'

middle = '05_NCP_track_120'
#middle = '12_NCP_track_126'
#middle = '13_NCP_track_135' # Why all NANs?
#middle = '16_NCP_track_208'
#middle = '17_NCP_track_132'

#### Functions ####

def mean_diag(power_mat):
    s = power_mat.shape
    n_nu = s[-1]
    power_reshaped = np.reshape(power_mat, s[:-2] + (n_nu**2,))
    out = np.mean(power_reshaped[...,::n_nu + 1], -1)
    return out

def log_bin(power, ratio=2.):
    """Only impormented on the 0th axis."""
    axis = 0
    n = power.shape[axis]
    if ratio <= 1:
        raise ValueError("`ratio` must be greater than unity.")
    edges = [0, float(ratio)]
    while edges[-1] < n - 1:
        edges = edges + [ratio * edges[-1]]
    s = power.shape[1:]
    new_n = len(edges) - 1
    brackets = []
    for ii in range(new_n):
        lower_ind = int(math.floor(edges[ii]))
        upper_ind = int(math.floor(edges[ii+1]))
        if lower_ind==upper_ind:
            continue
        brackets.append((lower_ind, upper_ind))
    new_n = len(brackets)
    out = np.zeros((new_n,) + s, dtype=power.dtype)
    for ii in range(new_n):
        lower_ind = brackets[ii][0]
        upper_ind = brackets[ii][1]
        out[ii,...] = np.mean(power[lower_ind:upper_ind,...], 0)
    return out


#### Analysis and plots ####

file_name = data_root + middle + end

# Initial read.
Reader = fitsGBT.Reader(file_name)
Data = Reader.read(0,0)

# Preprocess.
rotate_pol.rotate(Data, (-5, -7, -8, -6))
cal_scale.scale_by_cal(Data, True, False, False, False, True)
data_preflag = Data.data.copy()
flag_data.flag_data(Data, 5, 0.1, 40)
data_postflag = Data.data.copy()
Data.calc_freq()
nu_full = Data.freq
rebin_freq.rebin(Data, 16, True, True)
#rotate_pol.rotate(Data, (1,))
combine_cal.combine(Data, (0., 1.), False, True)
data_preflag2 = Data.data.copy()
reflag.flag(Data, Data, thres=3., modes_subtract=10,
            filter_type='gaussian/edge')
data_postflag2 = Data.data.copy()
Data.calc_freq()
nu = Data.freq
Data.calc_time()
time = Data.time

# Calculate the power spectrum.
DECON = False  # XXX.
tmp = noise_power.full_power_mat((Data,), n_time=-1, window='hanning',
                                 deconvolve=DECON, subtract_slope=True,
                                 normalize=not DECON, split_scans=False)
power_mat, window_function, dt, channel_means = tmp
# Calculate the thermal expectation level.
fudge = 2. #  To deal with only using cal-on, for instants.
d_nu = abs(Data.field['CDELT1'])
thermal_expectation = channel_means * np.sqrt(fudge / d_nu / dt) 
power_mat = power_mat / (thermal_expectation[None,:,:,None,:] 
                         * thermal_expectation[None,:,:,:,None])

# Some dimensions.
n_t = power_mat.shape[0]
n_nu = power_mat.shape[-1]
n_pol = power_mat.shape[1]
n_cal = power_mat.shape[2]
freq = noise_power.ps_freq_axis(dt, n_t)
n_f = len(freq)
# Prune to meaninful modes.
power_mat = noise_power.prune_power(power_mat, 0)
window_function = noise_power.prune_power(window_function, 0)

# Calculate the eigenmodes and subtract them out.
n_modes = 10
power_modes_removed = np.empty(power_mat.shape[:3] + (n_modes + 1,), 
                               dtype=np.float64)
all_modes = np.zeros(power_mat.shape[1:3] + (n_modes + 1,n_nu),
                     dtype=np.complex128)
for ii in range(n_pol):
    for jj in range(n_cal):
        reduced_power = power_mat[:,ii,jj,:,:].copy().real
        mat = np.mean(reduced_power[:n_f//4,:,:], 0)
        e, v = linalg.eigh(mat)
        power_modes_removed[:,ii,jj,0] = mean_diag(reduced_power).real
        for kk in range(1,n_modes + 1):
            mode = v[:,-kk]
            all_modes[ii,jj,kk,:] = mode
            # Calculate the components to remove.
            tmp_amp1 = np.sum(reduced_power * mode, -1)
            tmp_amp2 = np.sum(reduced_power * mode[:,None].conj(), -2)
            tmp_amp3 = np.sum(tmp_amp2 * mode, -1)
            # Remove them from the matrix.
            reduced_power -= tmp_amp1[:,:,None] * mode.conj()
            reduced_power -= tmp_amp2[:,None,:] * mode[:,None]
            reduced_power += (tmp_amp3[:,None,None] * mode[:,None] *
                              mode.conj())
            power_modes_removed[:,ii,jj,kk] = mean_diag(reduced_power).real

plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
for ii in range(7):
    plt.plot(log_bin(freq, 1.1), log_bin(power_modes_removed[:,0,0,ii], 1.1))
plt.xlim([5e-3, 4])
plt.ylim([1, 100])
plt.ylabel('noise power')
plt.xlabel('frequency (Hz)')

plt.figure()
base_line = 0
for ii in range(1, 8):
    plt.plot(nu/1e6, np.zeros(n_nu) + base_line, 'k:')
    plt.plot(nu/1e6, all_modes[0,0,ii,:].real + base_line)
    #plt.plot(nu/1e6, all_modes[0,0,ii,:].imag + base_line)
    base_line += 0.3
plt.ylabel('mode amplitude (normalized)')
plt.xlabel('spectral channel frequency (MHz)')

# Rebin preflagged data.
n_rebin = 16
pre = data_preflag[:,0,0,:]
pre.shape = (pre.shape[0], pre.shape[-1] // n_rebin, n_rebin)
pre = np.mean(pre, -1)
pre = (pre/np.mean(pre, 0)).filled(1) - 1
#pre = (data_preflag[:,0,0,:]/np.mean(data_preflag[:,0,0,:], 0)).filled(1) - 1
pre = pre[:,::-1]
# Rebin postflagged data.
#post = data_postflag[:,0,0,:]
#post.shape = (post.shape[0], post.shape[-1] // n_rebin, n_rebin)
#post = np.mean(post, -1)
post = data_postflag2[:,0,0,:]
post = (post/np.mean(post, 0)).filled(1) - 1
#post = (data_postflag[:,0,0,:]/np.mean(data_postflag[:,0,0,:], 0)).filled(1) - 1
post = post[:,::-1]
t = time - time[0]
#aspect = abs((nu[-1]/1e6 - nu[0]/1e6) / (t[-1] - t[0]))
aspect = 1
extent = [nu[-1]/1e6, nu[0]/1e6, t[-1], t[0]]
plt.figure()
scale_lim = 0.01
plt.imshow(pre, vmin=-scale_lim, vmax=scale_lim, extent=extent, aspect=aspect)
#plt.colorbar()
xlab = plt.xlabel('spectral channel frequency (MHz)')
ylab = plt.ylabel('time (s)')
plt.savefig('flag_before.eps', bbox_extra_artists=(xlab, ylab),
            bbox_inches='tight')
plt.figure()
plt.imshow(post, vmin=-scale_lim, vmax=scale_lim, extent=extent, aspect=aspect)
plt.colorbar()
xlab = plt.xlabel('spectral channel frequency (MHz)')
ylab = plt.ylabel('time (s)')
plt.savefig('flag_after.eps', bbox_extra_artists=(xlab, ylab),
            bbox_inches='tight')



plt.show()
