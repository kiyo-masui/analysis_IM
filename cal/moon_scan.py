import os
import sys
import numpy as np
import scipy as sp
import ephem
from core import fitsGBT
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import combine_cal
import pickle
import utils.misc as utils
import matplotlib.pyplot as plt
"""This is some quick code analyze moon scan data"""

def load_moonscan(filename):
    cal_coords = ephem.Equatorial("05:42:36.155", "+49:51:07.28",
                                  epoch=ephem.B1950)

    # convert cal to a Body object.
    cal_source = ephem.FixedBody()
    cal_source._ra = cal_coords.ra
    cal_source._dec = cal_coords.dec
    cal_source._epoch = cal_coords.epoch

    Reader = fitsGBT.Reader(filename)
    moon_dataobj = Reader.read(0,0)

    rotate_pol.rotate(moon_dataobj, (-5, -7, -8, -6))
    cal_scale.scale_by_cal(moon_dataobj, True, False, False, False, True)
    flag_data.flag_data(moon_dataobj, 5, 0.1, 40)
    rebin_freq.rebin(moon_dataobj, 16, True, True)
    #rebin_time.rebin(moon_dataobj, 4)

    return moon_dataobj


def plot_obs(time_data):
    fig = plt.figure(1, figsize = (20,2))
    ax = fig.add_subplot(111)

    cax = ax.imshow(time_data, interpolation='nearest')
    cbar = fig.colorbar(cax)

    plt.show()


def process_moonscan(moon_dataobj, outfile, avg_width=30):
    (ndata, npol, ntoggle, nfreq) = moon_dataobj.dims

    moon_dataobj.calc_freq()
    moon_dataobj.calc_pointing()
    moon_dataobj.calc_time()
    freq = moon_dataobj.freq
    ra = moon_dataobj.ra
    dec = moon_dataobj.dec
    date_time = moon_dataobj.field['DATE-OBS']
    time = moon_dataobj.time
    az = moon_dataobj.field['CRVAL2']
    el = moon_dataobj.field['CRVAL3']

    print moon_dataobj.field.keys()
    print moon_dataobj.field['CAL']

    pols = list(moon_dataobj.field['CRVAL4'])
    #print moon_dataobj.field['CRVAL4']
    pol_names = {}
    for pol_idx in range(npol):
        pol_names[utils.polint2str(pols[pol_idx])] = pol_idx

    print pol_names

    # find the peak in XX, YY
    maxlist = np.argmax(moon_dataobj.data[:, 0, 0, :], axis=0)
    maxlist = maxlist[np.where(maxlist > 0)]
    max_xx = np.mean(maxlist)
    print "max_xx", max_xx, np.std(maxlist)

    maxlist = np.argmax(moon_dataobj.data[:, 3, 0, :], axis=0)
    maxlist = maxlist[np.where(maxlist > 0)]
    max_yy = np.mean(maxlist)
    print "max_yy", max_yy, np.std(maxlist)

    time_max = int((max_xx + max_yy) / 2.)
    print time_max

    # find the max/min in YX (well-defined because of dipole leakage)
    maxlist = np.argmax(moon_dataobj.data[:, 2, 0, :], axis=0)
    maxlist = maxlist[np.where(maxlist > 0)]
    max_yx = int(np.mean(maxlist))
    print "max_yx", max_yx, np.std(maxlist)

    minlist = np.argmin(moon_dataobj.data[:, 2, 0, :], axis=0)
    minlist = minlist[np.where(minlist > 0)]
    min_yx = int(np.mean(minlist))
    print "min_yx", min_yx, np.std(minlist)

    before = np.mean(moon_dataobj.data[0:(2 * avg_width), :, 0, :], axis=0)
    during = np.mean(moon_dataobj.data[(time_max - avg_width): \
                                       (time_max + avg_width), :, 0, :], axis=0)
    during_u = np.mean(moon_dataobj.data[(max_yx - avg_width): \
                                       (max_yx + avg_width), :, 0, :], axis=0)
    during_d = np.mean(moon_dataobj.data[(min_yx - avg_width): \
                                       (min_yx + avg_width), :, 0, :], axis=0)
    after = np.mean(moon_dataobj.data[-(2 * avg_width):-1, :, 0, :], axis=0)

    offsrc = (before + after) / 2.
    during -= offsrc
    during_u -= offsrc
    during_d -= offsrc

    outfp = open(outfile, "w")

    for freq_idx in range(nfreq):
        compilation = [freq_idx]
        #compilation.extend(before[:, freq_idx])
        compilation.extend(during[:, freq_idx])
        compilation.extend(during_u[:, freq_idx])
        compilation.extend(during_d[:, freq_idx])
        #compilation.extend(after[:, freq_idx])
        outfp.write(("%.5g " * 13) % tuple(compilation) + "\n")

    outfp.close()

    #for pol_idx in range(npol):
    #    time_data = moon_dataobj.data[:, pol_idx, 0, :]
    #    #plot_obs(time_data.transpose())
    #    maxlist = np.argmax(time_data, axis=0)
    #    maxlist = maxlist[np.where(maxlist > 0)]
    #    print pol_idx, "max", np.mean(maxlist), np.std(maxlist)

    #    minlist = np.argmin(time_data, axis=0)
    #    minlist = minlist[np.where(minlist > 0)]
    #    print pol_idx, "min", np.mean(minlist), np.std(minlist)

    #for ii in range(moon_dataobj.dims[0])


def process_moon(filename, outfile):
    moon_data = load_moonscan(filename)
    outpkl = open('data.pkl', 'wb')
    pickle.dump(moon_data, outpkl, -1)
    outpkl.close()

    inpkl = open('data.pkl', 'rb')
    moon_data = pickle.load(inpkl)
    inpkl.close()

    process_moonscan(moon_data, outfile)


def find_obsvec(coherency, alternate=True):
    obsvec = np.zeros((coherency.shape[0], 4))

    obsvec[:, 0] = coherency[:, 0, 0]
    obsvec[:, 1] = 0.5 * (coherency[:, 1, 0] + coherency[:, 0, 1])
    obsvec[:, 2] = 0.5 * (coherency[:, 1, 0] - coherency[:, 0, 1]) * -1.j
    obsvec[:, 3] = coherency[:, 1, 1]

    if alternate:
        coh_xy = 0.5 * (obsvec[:, 1] + obsvec[:, 2])
        coh_yx = 0.5 * (obsvec[:, 1] - obsvec[:, 2])
        obsvec[:, 1] = coh_xy
        obsvec[:, 2] = coh_yx

    return obsvec


def find_coherency(filename, alternate=False):
    fp = open(filename, "r")
    data = np.genfromtxt(fp, names=["freq", "XX", "XY", "YX", "YY", \
                                              "XXu", "XYu", "YXu", "YYu", \
                                              "XXd", "XYd", "YXd", "YYd"])
    fp.close()

    coherency = np.zeros((data["freq"].shape[0], 2, 2), dtype=complex)

    stokes_i = 0.5 * (data["XX"] + data["YY"])
    stokes_q = 0.5 * (data["XX"] - data["YY"])
    if alternate:
        stokes_u = data["XY"] + data["YX"]
        stokes_v = data["XY"] - data["YX"]
    else:
        stokes_u = data["XY"]
        stokes_v = data["YX"]

    coherency[:, 0, 0] = stokes_i + stokes_q
    coherency[:, 1, 0] = stokes_u + stokes_v * 1.j
    coherency[:, 0, 1] = stokes_u - stokes_v * 1.j
    coherency[:, 1, 1] = stokes_i - stokes_q

    obsvec = find_obsvec(coherency, alternate=alternate)
    print obsvec[:, 0] - data["XX"]
    print obsvec[:, 1] - data["XY"]
    print obsvec[:, 2] - data["YX"]
    print obsvec[:, 3] - data["YY"]

    return coherency


def verify_moon(filename1, filename2, outfile, normalize=True):
    coherency1 = find_coherency(filename1)
    coherency2 = find_coherency(filename2)
    coherency1out = np.zeros_like(coherency1)
    coherency2out = np.zeros_like(coherency2)

    jones = np.zeros_like(coherency1)
    jones_inv = np.zeros_like(coherency1)

    num_freq = coherency1.shape[0]
    for idx in range(num_freq):
        print coherency1[idx, :, :]
        try:
            #jones[idx, :, :] = np.linalg.cholesky(coherency1[idx, :, :])
            jones[idx, :, :] = sp.linalg.sqrtm(coherency1[idx, :, :])
        except (np.linalg.LinAlgError, ValueError):
            coherency2out[idx, :, :] = np.zeros((2,2))
            print "not positive definite"
            continue

        print jones[idx, :, :]
        print np.linalg.det(jones[idx, :, :])

        jones[idx, :, :] /= np.sqrt(np.linalg.det(jones[idx, :, :]))

        jones_inv[idx, :, :] = np.linalg.inv(jones[idx, :, :])
        print jones_inv[idx, :, :]
        print np.dot(jones_inv[idx, :, :], np.dot(coherency1[idx, :, :],
                     jones_inv[idx, :, :].T.conj()))
        coherency2out[idx, :, :] = np.dot(jones_inv[idx, :, :],
                                          np.dot(coherency2[idx, :, :],
                                          jones_inv[idx, :, :].T.conj()))
        coherency1out[idx, :, :] = np.dot(jones_inv[idx, :, :],
                                          np.dot(coherency1[idx, :, :],
                                          jones_inv[idx, :, :].T.conj()))
        print coherency2out[idx, :, :]
        print "-"*80

    obsvec_ref = find_obsvec(coherency1)
    obsvec_in = find_obsvec(coherency2)
    obsvec_out = find_obsvec(coherency2out)
    fout = open(outfile, "w")
    for idx in range(num_freq):
        obs = np.append(obsvec_ref[idx, :], [obsvec_in[idx, :],
                         obsvec_out[idx, :]])
        fout.write(("%10.15g " * 12) % tuple(obs) + "\n")

    fout.close()


def batch_moon(generate_fluxes=False):
    datapath = "/mnt/raid-project/gmrt/kiyo/data/guppi_data/GBT12A_418/"
    filename1 = datapath + "22_MOON_track_54.fits"
    filename2 = datapath + "22_MOON_track_10.fits"
    fluxfile1 = "moon_track_54.txt"
    fluxfile2 = "moon_track_10.txt"

    if generate_fluxes:
        process_moon(filename1, fluxfile1)
        process_moon(filename2, fluxfile2)

    verify_moon(fluxfile1, fluxfile2, "recovered.dat")

batch_moon()
