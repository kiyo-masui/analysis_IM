'''
wigglez_xcorr: cross correlate GBT with WiggleZ data
'''
# TODO: for some reason, PYTHONPATH is clobbered on the compute nodes; it has
# the local directory as ::, but this never seems to make it to the sys.path,
# so add it explicitly. Should contact admin.?
import sys, site
site.addsitedir('./')

import numpy as np
import scipy as sp
from core import algebra
#import matplotlib
#matplotlib.use('Agg')
from correlate import freq_slices as fs
from kiyopy import parse_ini
import shelve
from map import beam
# TODO: this seemed to be necessary on the Sunnyvale compute nodes because it
# was clobbing the python path?
import sys, site
site.addsitedir('/home/eswitzer/local/lib/')
site.addsitedir('/home/eswitzer/local/lib/python2.6/site-packages/')

params_default = {
      'optical_root': '/cita/h/home-2/eswitzer/data/binned_wigglez/',
      'radio_root': '/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/',
      'radio_data_file': 'sec_A_15hr_41-69_cleaned_clean_map_I.npy',
      'radio_noiseinv_file': 'sec_A_15hr_41-69_cleaned_noise_inv_I.npy',
      'optical_data_file': 'reg15data.npy',
      'optical_selection_file': 'reg15separable.npy',
      'freq': (),
      'lags': (),
      'output_shelve_file': 'test.shelve',
      'convolve': False,
      'subtract_mean': True,
      'speedup': False
      }
prefix = 'fs_'


def wigglez_correlation(init_filename):
    """Perform the cross-correlation between the WiggleZ dataset and GBT
    """
    np.set_printoptions(threshold=np.nan)

    print "starting wigglez correlation run with file: " + init_filename
    params = parse_ini.parse(init_filename, params_default,
                             prefix=prefix, feedback=10)

    radio_file = params['radio_root'] + params['radio_data_file']
    noiseinv_file = params['radio_root'] + params['radio_noiseinv_file']
    optical_file = params['optical_root'] + params['optical_data_file']
    optical_selection_file = params['optical_root'] + \
                             params['optical_selection_file']

    map_radio = algebra.make_vect(algebra.load(radio_file))
    #noiseinv_radio = np.abs(algebra.make_vect(algebra.load(noiseinv_file)))
    noiseinv_radio = algebra.make_vect(algebra.load(noiseinv_file))
    map_opt = algebra.make_vect(algebra.load(optical_file))
    map_nbar = algebra.make_vect(algebra.load(optical_selection_file))

    algebra.compressed_array_summary(map_opt, "opt map as loaded")
    algebra.compressed_array_summary(map_nbar, "nbar map as loaded")
    algebra.compressed_array_summary(map_radio, "radio map as loaded")
    algebra.compressed_array_summary(noiseinv_radio, "radio N^-1 as loaded")

    noiseinv_radio[noiseinv_radio < 1.e-20] = 0.
    nan_array = np.isnan(noiseinv_radio)
    noiseinv_radio[nan_array] = 0.
    inf_array = np.isinf(noiseinv_radio)
    noiseinv_radio[inf_array] = 0.

    if params['convolve']:
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                              0.281176247549, 0.270856788455, 0.26745856078,
                              0.258910010848, 0.249188429031])
        # also consider the beam model where everything is degraded to the
        # lowest resolution; e.g. in a MapPair
        beam_data_degrade = sp.zeros(8) + max(beam_data)
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                          dtype=float)
        freq_data *= 1.e6
        psf = beam.GaussianBeam(beam_data_degrade, freq_data)
        # Convolve the optical data to the lowest radio resolution
        print "convolving the map"
        map_opt = psf.apply(map_opt)
        algebra.compressed_array_summary(map_opt, "opt after convolution")

        # how should the covariance be convolved?
        # nbar is N^-1; convolve B N B^T?,
        # or just convolve nbar? use?:
        #map_nbar[map_nbar<1.e-30] = 1.e-30
        #map_nbar = 1./map_nbar
        #map_nbar = psf.apply(map_nbar, cval=1.e30)
        #map_nbar = 1./map_nbar
        #map_nbar[map_nbar<1.e-20] = 0
        print "convolving the covariance/selection"
        map_nbar = psf.apply(map_nbar)
        algebra.compressed_array_summary(map_nbar, "nbar after convolution")

    # convert to delta-overdensity
    map_opt = map_opt / map_nbar - 1.
    #algebra.compressed_array_summary(map_opt, "opt after conversion to delta")

    # set the NaNs and infs to zero in data and weights
    nan_array = np.isnan(map_opt)
    map_opt[nan_array] = 0.
    map_nbar[nan_array] = 0.
    inf_array = np.isinf(map_opt)
    map_opt[inf_array] = 0.
    map_nbar[inf_array] = 0.

    freqlist = params['freq']  # full: range(map_radio.shape[0])

    algebra.compressed_array_summary(map_opt[freqlist, :, :],
                    "opt map as entering the correlation function")
    algebra.compressed_array_summary(map_nbar[freqlist, :, :],
                    "nbar as entering the correlation function")
    algebra.compressed_array_summary(map_radio[freqlist, :, :],
                    "radio map as entering the correlation function")
    algebra.compressed_array_summary(noiseinv_radio[freqlist, :, :],
                    "radio map N^-1 as entering the correlation function")

    cross_pair = fs.MapPair(map_opt, map_radio, map_nbar, noiseinv_radio,
                           freqlist)

    if params['subtract_mean']:
        cross_pair.subtract_weighted_mean()

    (corr, counts) = cross_pair.correlate(params['lags'],
                            speedup=params['speedup'])

    corr_shelve = shelve.open(params['output_shelve_file'])
    corr_shelve["corr"] = corr
    corr_shelve["counts"] = counts
    corr_shelve["freq_axis"] = map_radio.get_axis('freq')
    corr_shelve["params"] = params
    corr_shelve.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        initfilename = str(sys.argv[1])
        wigglez_correlation(initfilename)
    else:
        print 'Maximum one argument, a parameter file name.'
