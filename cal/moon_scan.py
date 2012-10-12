import os
import sys
import numpy as np
import scipy as sp
import ephem
from core import fitsGBT
from time_stream import rotate_pol, cal_scale, flag_data, rebin_freq
from time_stream import combine_cal, moon_rotation
import pickle
import utils.misc as utils
import matplotlib.pyplot as plt
import h5py
"""This is some quick code analyze moon scan data"""

def load_moonscan(filename, rotate_moon=True):
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
    cal_scale.scale_by_cal(moon_dataobj, scale_t_ave=True, scale_f_ave=False,
                           sub_med=False, scale_f_ave_mod=False, rotate=True)
    flag_data.flag_data(moon_dataobj, 5, 0.1, 40)
    rebin_freq.rebin(moon_dataobj, 16, True, True)
    #rebin_time.rebin(moon_dataobj, 4)

    if rotate_moon:
        moon_rotation.rotate_pol_moon(moon_dataobj)

    fgc_mueler_file = '/mnt/raid-project/gmrt/tcv/diff_gain_params/GBT12A_418/22_diff_gain_calc.txt'
    fgc_RM_file = ' '
    fgc_R_to_sky = True
    fgc_DP_correct = False  # this is already handled in scale_by_cal's rotate
    fgc_RM_correct = False

    from time_stream import flux_diff_gain_cal as fdg
    m_total = fdg.flux_dg(fgc_mueler_file)
    fdg.calibrate_pol(moon_dataobj, m_total, fgc_RM_file,
                      fgc_R_to_sky, fgc_DP_correct, fgc_RM_correct)

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


def find_coherency(pol_vec, alternate=False):
    coherency = np.zeros((2,2), dtype=complex)

    stokes_i = 0.5 * (pol_vec[0] + pol_vec[3])
    stokes_q = 0.5 * (pol_vec[0] - pol_vec[3])
    stokes_u = pol_vec[1]
    stokes_v = pol_vec[2]

    if alternate:
        stokes_u = pol_vec[1] + pol_vec[2]
        stokes_v = pol_vec[1] - pol_vec[2]

    coherency[0, 0] = stokes_i + stokes_q
    coherency[1, 0] = stokes_u + stokes_v * 1.j
    coherency[0, 1] = stokes_u - stokes_v * 1.j
    coherency[1, 1] = stokes_i - stokes_q

    return coherency


def find_polvec(coherency, alternate=False):
    #pol_vecout = np.zeros(4, dtype=complex)
    pol_vecout = np.zeros(4)

    pol_vecout[0] = coherency[0, 0]
    pol_vecout[1] = 0.5 * (coherency[1, 0] + coherency[0, 1])
    pol_vecout[2] = 0.5 * (coherency[1, 0] - coherency[0, 1]) * -1.j
    pol_vecout[3] = coherency[1, 1]
    # alternate form
    if alternate:
        coh_xy = 0.5 * (pol_vecout[1] + pol_vecout[2])
        coh_yx = 0.5 * (pol_vecout[1] - pol_vecout[2])
        pol_vecout[1] = coh_xy
        pol_vecout[2] = coh_yx

    return pol_vecout


def rotate_single(pol_vec, rotation):
    r"""transform an observed polarization vector using a Jones matrix"""
    coherency = find_coherency(pol_vec)
    coherency = np.dot(rotation, np.dot(coherency, rotation.T.conj()))

    return find_polvec(coherency)


def rotate_mueller(pol_vec, rotation):
    r"""transform an observed polarization vector using an analogous mueller
    matrix"""

    return np.dot(pol_vec, rotation)


def find_rotation(pol_vec):
    r"""Factorize polarization response of observations of an unpolarized
    source"""
    coherency = find_coherency(pol_vec)
    jones = np.zeros_like(coherency)
    jones_inv = np.zeros_like(coherency)

    try:
        #jones = np.linalg.cholesky(coherency)
        jones = sp.linalg.sqrtm(coherency)
    except (np.linalg.LinAlgError, ValueError):
        jones = np.zeros((2,2))
        print "not positive definite"

    jones /= np.sqrt(np.linalg.det(jones))

    jones_inv = np.linalg.inv(jones)
    print coherency, jones, np.linalg.det(jones), jones_inv

    return jones_inv


def convert_to_matrix(rotation):
    r"""convert a jones matrix to a 4x4 rotation on the observed pols"""
    mueller = np.zeros((4,4))
    mueller[0, :] = rotate_single(np.array([1., 0., 0., 0.]), rotation)
    mueller[1, :] = rotate_single(np.array([0., 1., 0., 0.]), rotation)
    mueller[2, :] = rotate_single(np.array([0., 0., 1., 0.]), rotation)
    mueller[3, :] = rotate_single(np.array([0., 0., 0., 1.]), rotation)

    print mueller
    return mueller


def open_moon_obs(filename, alternate=False):
    fp = open(filename, "r")
    data = np.genfromtxt(fp, names=["freq", "XX", "XY", "YX", "YY", \
                                              "XXu", "XYu", "YXu", "YYu", \
                                              "XXd", "XYd", "YXd", "YYd"])
    print data.dtype.names
    fp.close()

    print data['freq']
    return data


def verify_moon(filename1, filename2, outfile, rot_file, normalize=True):
    observation1 = open_moon_obs(filename1)
    observation2 = open_moon_obs(filename2)
    summaryfile = h5py.File(rot_file, "w")
    rotation_array = np.zeros((4, 4, 256))

    num_freq = observation1['freq'].shape[0]
    fout = open(outfile, "w")
    for idx in range(num_freq):
        obsvec_ref = np.array([observation1['XX'][idx],
                               observation1['XY'][idx],
                               observation1['YX'][idx],
                               observation1['YY'][idx]])

        obsvec_in = np.array([observation2['XX'][idx],
                               observation2['XY'][idx],
                               observation2['YX'][idx],
                               observation2['YY'][idx]])

        print "-"*80
        jones_inv = find_rotation(obsvec_ref)
        mueller = convert_to_matrix(jones_inv)
        rotation_array[:, :, idx] = mueller
        obsvec_out1 = rotate_single(obsvec_in, jones_inv)
        obsvec_out2 = rotate_mueller(obsvec_in, mueller)
        print "this", obsvec_out2 - obsvec_out1
        print "-"*80

        print obsvec_ref, obsvec_in
        obs = np.append(obsvec_ref, [obsvec_in, obsvec_out2])
        fout.write(("%10.15g " * 12) % tuple(obs) + "\n")

    summaryfile['moon_scan1'] = rotation_array
    summaryfile.close()
    fout.close()


def batch_moon(generate_fluxes=False):
    datapath = "/mnt/raid-project/gmrt/kiyo/data/guppi_data/GBT12A_418/"
    filename1 = datapath + "22_MOON_track_54.fits"
    filename2 = datapath + "22_MOON_track_10.fits"
    fluxfile1 = "moon_track_54_rot.txt"
    fluxfile2 = "moon_track_10_rot.txt"

    if generate_fluxes:
        process_moon(filename1, fluxfile1)
        process_moon(filename2, fluxfile2)

    file_rot = "/mnt/raid-project/gmrt/eswitzer/GBT/calibration/"
    file_rot += "moon_rotation_hyb.hd5"
    verify_moon(fluxfile1, fluxfile2, "recovered_new.dat", file_rot)


batch_moon()
