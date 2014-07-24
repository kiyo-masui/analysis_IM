import sys

import numpy as np
#import matplotlib.pyplot as plt

from core import fitsGBT



def bandpass_I_from_onoff(filename):
    reader = fitsGBT.Reader(f)
    scans = reader.read()
    scans[0].calc_freq()
    freq = scans[0].freq
    src = scans[0].field['OBJECT']
    on_spectrum = I_in_cal_units(scans[0])
    off_spectrum = I_in_cal_units(scans[1])
    src_spectrum = on_spectrum - off_spectrum

    src_nominal = calibrator_spectrum(src, freq)

    cal_spectrum = src_nominal / src_spectrum

    out = np.empty(len(freq),
            dtype=[('freq', np.float64), ('cal_T', np.float64)])
    out['freq'] = freq
    out['cal_T'] = cal_spectrum

    #plt.figure()
    #plt.plot(freq, on_spectrum)
    #plt.plot(freq, off_spectrum)

    #plt.figure()
    #plt.plot(freq, cal_spectrum)

    return out


def I_in_cal_units(scan):
    cal_on = np.mean(scan.data[:,0,0,:], 0).filled(0)
    cal_off = np.mean(scan.data[:,0,1,:], 0).filled(0)
    bad_inds = cal_on == cal_off
    cal_on[bad_inds] = 1.
    cal_off[bad_inds] = 0
    power = (cal_on + cal_off) / (cal_on - cal_off)
    power[bad_inds] = float('nan')
    return power


def calibrator_spectrum(src, freq):
    gain = 2    # K/Jy
    # Data from arXiv:1211.1300v1.
    if src == '3C48':
        coeff = [1.3324, -0.7690, -0.1950, 0.059]
    elif src == '3C295':
        coeff = [1.4866, -0.7871, -0.3440, 0.0749]

    l_freq_ghz = np.log10(freq / 1e9)    # Convert to GHz.
    poly = 1
    spec = 0
    for A in coeff:
        spec += A * poly
        poly *= l_freq_ghz
    spec = 10.**spec
    #plt.figure()
    #plt.plot(freq, spec)
    return spec




if __name__ == "__main__":
    files = sys.argv[1:]

    bp = []
    for f in files:
        bp.append(bandpass_I_from_onoff(f))

    n = len(bp)
    cal_T_n = np.empty((n, len(bp[0])), dtype=np.float64)
    freq = bp[0]['freq']
    for ii in range(n):
        cal_T_n[ii,:] = bp[ii]['cal_T']
        if not np.allclose(bp[ii]['freq'], freq):
            raise ValueError("Frequency axis")

    bad_chans = np.any(np.isnan(cal_T_n), 0)
    cal_T_n[:,bad_chans] = 1
    med_spec = np.median(cal_T_n, 0)
    std_spec = np.std(cal_T_n, 0)
    norm_std = std_spec / med_spec
    bad_chans = np.logical_or(bad_chans, norm_std > 10 * np.median(norm_std))
    med_spec[bad_chans] = float('nan')

    #plt.figure()
    #plt.plot(freq, med_spec)

    #plt.show()

    out = bp[0]
    out['cal_T'] = med_spec
    np.save("cal_spectrum_I.npy", out)

