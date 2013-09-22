'''
radio_corr: cross radio volumes; this is specifically for simulations where the
treatment is different than the real radio autocorr in freq slices
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
      'radio_noiseroot1': '/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/',
      'radio_noiseroot2': '/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/',
      'radio_root1': '/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/',
      'radio_root2': '/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/',
      'radio_data_file1': 'sec_A_15hr_41-69_cleaned_clean_map_I.npy',
      'radio_noiseinv_file1': 'sec_A_15hr_41-69_cleaned_noise_inv_I.npy',
      'radio_data_file2': 'sec_A_15hr_41-69_cleaned_clean_map_I.npy',
      'radio_noiseinv_file2': 'sec_A_15hr_41-69_cleaned_noise_inv_I.npy',
      'freq': (),
      'lags': (),
      'output_shelve_file': 'test.shelve',
      'convolve': False,
      'subtract_mean': True,
      'speedup': False
      }
prefix = 'fs_'


def radio_correlation(init_filename):
    """Perform the cross-correlation between two mock GBT cubes
    """
    np.set_printoptions(threshold=np.nan)

    print "starting wigglez correlation run with file: " + init_filename
    params = parse_ini.parse(init_filename, params_default,
                             prefix=prefix, feedback=10)

    radio_file1 = params['radio_root1'] + params['radio_data_file1']
    noiseinv_file1 = params['radio_noiseroot1'] + params['radio_noiseinv_file1']
    radio_file2 = params['radio_root2'] + params['radio_data_file2']
    noiseinv_file2 = params['radio_noiseroot2'] + params['radio_noiseinv_file2']

    map_radio1 = algebra.make_vect(algebra.load(radio_file1))
    noiseinv_radio1 = algebra.make_vect(algebra.load(noiseinv_file1))
    map_radio2 = algebra.make_vect(algebra.load(radio_file2))
    noiseinv_radio2 = algebra.make_vect(algebra.load(noiseinv_file2))

    algebra.compressed_array_summary(map_radio1, "radio 1 map as loaded")
    algebra.compressed_array_summary(noiseinv_radio1, "radio 1 N^-1 as loaded")
    algebra.compressed_array_summary(map_radio2, "radio 2 map as loaded")
    algebra.compressed_array_summary(noiseinv_radio2, "radio 2 N^-1 as loaded")

    noiseinv_radio1[noiseinv_radio1 < 1.e-20] = 0.
    nan_array = np.isnan(noiseinv_radio1)
    noiseinv_radio1[nan_array] = 0.
    inf_array = np.isinf(noiseinv_radio1)
    noiseinv_radio1[inf_array] = 0.

    noiseinv_radio2[noiseinv_radio2 < 1.e-20] = 0.
    nan_array = np.isnan(noiseinv_radio2)
    noiseinv_radio2[nan_array] = 0.
    inf_array = np.isinf(noiseinv_radio2)
    noiseinv_radio2[inf_array] = 0.

    freqlist = params['freq']  # full: range(map_radio.shape[0])

    algebra.compressed_array_summary(map_radio1[freqlist, :, :],
                    "radio map 1 as entering the correlation function")
    algebra.compressed_array_summary(noiseinv_radio1[freqlist, :, :],
                    "radio map 1 N^-1 as entering the correlation function")
    algebra.compressed_array_summary(map_radio2[freqlist, :, :],
                    "radio map 2 as entering the correlation function")
    algebra.compressed_array_summary(noiseinv_radio2[freqlist, :, :],
                    "radio map 2 N^-1 as entering the correlation function")

    cross_pair = fs.MapPair(map_radio1, map_radio2,
                            noiseinv_radio1, noiseinv_radio2,
                            freqlist)

    if params['subtract_mean']:
        cross_pair.subtract_weighted_mean()

    if params['convolve']:
        print "WARNING: you are degrading the beams"
        cross_pair.degrade_resolution()

    (corr, counts) = cross_pair.correlate(params['lags'],
                            speedup=params['speedup'])

    corr_shelve = shelve.open(params['output_shelve_file'])
    corr_shelve["corr"] = corr
    corr_shelve["counts"] = counts
    corr_shelve["freq_axis"] = map_radio1.get_axis('freq')
    corr_shelve["params"] = params
    corr_shelve.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        initfilename = str(sys.argv[1])
        radio_correlation(initfilename)
    else:
        print 'Maximum one argument, a parameter file name.'
