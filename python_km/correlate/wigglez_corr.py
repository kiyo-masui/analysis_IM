'''
wigglez_corr: WiggleZ x WiggleZ 
'''
# TODO: for some reason, PYTHONPATH is clobbered on the compute nodes; it has
# the local directory as ::, but this never seems to make it to the sys.path,
# so add it explicitly. Should contact admin.?
import sys, site
site.addsitedir('./')

import numpy as np
import scipy as sp
from core import algebra
import operator
import matplotlib
matplotlib.use('Agg')
from correlate import freq_slices as fs
import correlate.correlation_plots as cp
from kiyopy import parse_ini
import shelve
from map import beam
# TODO: this seemed to be necessary on the Sunnyvale compute nodes because it
# was clobbing the python path?
import sys, site
site.addsitedir('/home/eswitzer/local/lib/')
site.addsitedir('/home/eswitzer/local/lib/python2.6/site-packages/')

# TODO: delete /mnt/raid-project/gmrt/eswitzer/wiggleZ/maps/sec_A_15hr_41-73_clean_map_I.npy

params_default = {
      'optical_root1': '/cita/h/home-2/eswitzer/data/binned_wigglez/',
      'optical_root2': '/cita/h/home-2/eswitzer/data/binned_wigglez/',
      'optical_data_file1': 'reg15data.npy',
      'optical_selection_file1': 'reg15separable.npy',
      'optical_data_file2': 'reg15data.npy',
      'optical_selection_file2': 'reg15separable.npy',
      'freq': (),
      'lags': (),
      'output_shelve_file': 'test.shelve',
      'convolve': False,
      'subtract_mean': True,
      'speedup': False
      }
prefix = 'fs_'


def array_summary(array, testname, axes, meetall=False, identify_entries=True):
    """helper function for summarizing arrays
    meetall: prints those entries for which all values in the slice meet the
    criteria (normal behavior is print all entries where _any_ value in the
    slice meets the criteria
    identify_entries: prints entries meeting the criteria
    """
    total_matching = array.sum()
    if total_matching != 0:
        print testname + "s:"
        match_count = np.apply_over_axes(np.sum, array, axes)
        print match_count.flatten()
        if identify_entries:
            if meetall:
                arrayshape = array.shape
                subarray_size = reduce(operator.mul,
                                       [arrayshape[i] for i in axes])
                print "with all " + testname + "s: " + \
                      repr(np.where(match_count.flatten() == subarray_size))
            else:
                print "has " + testname + ": " + \
                      repr(np.where(match_count.flatten() != 0))

            print "total " + testname + "s: " + repr(total_matching)

        print "-" * 80
    else:
        print "There are no " + testname + " entries"


# TODO: move this somewhere for utilities
def compressed_array_summary(array, name, axes=[1, 2], extras=False):
    """print various summaries of arrays compressed along specified axes"""

    print "-" * 80
    print "array property summary for " + name + ":"
    array_summary(np.isnan(array), "nan", axes)
    array_summary(np.isinf(array), "inf", axes)
    array_summary((array == 0.), "zero", axes, meetall=True)
    array_summary((array < 0.), "negative", axes, identify_entries=False)

    if extras:
        sum_nu = np.apply_over_axes(np.sum, array, axes)
        min_nu = np.apply_over_axes(np.min, array, axes)
        max_nu = np.apply_over_axes(np.max, array, axes)
        print sum_nu.flatten()
        print min_nu.flatten()
        print max_nu.flatten()

    print ""


def wigglez_correlation(init_filename):
    """Perform the cross-correlation between the WiggleZ dataset and GBT
    """
    np.set_printoptions(threshold=np.nan)

    print "starting wigglez correlation run with file: " + init_filename
    params = parse_ini.parse(init_filename, params_default,
                             prefix=prefix, feedback=10)

    optical_file1 = params['optical_root1'] + params['optical_data_file1']
    optical_selection_file1 = params['optical_root1'] + \
                             params['optical_selection_file1']
    optical_file2 = params['optical_root2'] + params['optical_data_file2']
    optical_selection_file2 = params['optical_root2'] + \
                             params['optical_selection_file2']

    map_opt1 = algebra.make_vect(algebra.load(optical_file1))
    map_nbar1 = algebra.make_vect(algebra.load(optical_selection_file1))
    map_opt2 = algebra.make_vect(algebra.load(optical_file2))
    map_nbar2 = algebra.make_vect(algebra.load(optical_selection_file2))

    compressed_array_summary(map_opt1, "opt map 1 as loaded")
    compressed_array_summary(map_nbar1, "nbar map 1 as loaded")
    compressed_array_summary(map_opt2, "opt map 2 as loaded")
    compressed_array_summary(map_nbar2, "nbar map 2 as loaded")


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
        print "convolving the first map"
        map_opt1 = psf.apply(map_opt1)
        compressed_array_summary(map_opt1, "opt 1 after convolution")
        print "convolving the second map"
        map_opt2 = psf.apply(map_opt2)
        compressed_array_summary(map_opt2, "opt 2 after convolution")

        # how should the covariance be convolved?
        # nbar is N^-1; convolve B N B^T?,
        # or just convolve nbar? use?:
        #map_nbar[map_nbar<1.e-30] = 1.e-30
        #map_nbar = 1./map_nbar
        #map_nbar = psf.apply(map_nbar, cval=1.e30)
        #map_nbar = 1./map_nbar
        #map_nbar[map_nbar<1.e-20] = 0
        print "convolving the first covariance/selection"
        map_nbar1 = psf.apply(map_nbar1)
        compressed_array_summary(map_nbar1, "nbar 1 after convolution")
        print "convolving the second covariance/selection"
        map_nbar2 = psf.apply(map_nbar2)
        compressed_array_summary(map_nbar2, "nbar 2 after convolution")

    # convert to delta-overdensity
    map_opt1 = map_opt1 / map_nbar1 - 1.
    map_opt2 = map_opt2 / map_nbar2 - 1.
    #compressed_array_summary(map_opt1, "opt 1 after conversion to delta")
    #compressed_array_summary(map_opt2, "opt 2 after conversion to delta")

    # set the NaNs and infs to zero in data and weights
    nan_array = np.isnan(map_opt1)
    map_opt1[nan_array] = 0.
    map_nbar1[nan_array] = 0.
    inf_array = np.isinf(map_opt1)
    map_opt1[inf_array] = 0.
    map_nbar1[inf_array] = 0.

    nan_array = np.isnan(map_opt2)
    map_opt2[nan_array] = 0.
    map_nbar2[nan_array] = 0.
    inf_array = np.isinf(map_opt2)
    map_opt2[inf_array] = 0.
    map_nbar2[inf_array] = 0.

    freqlist = params['freq']  # full: range(map_radio.shape[0])

    compressed_array_summary(map_opt1[freqlist, :, :],
                    "opt map 1 as entering the correlation function")
    compressed_array_summary(map_nbar1[freqlist, :, :],
                    "nbar 1 as entering the correlation function")
    compressed_array_summary(map_opt2[freqlist, :, :],
                    "opt map 2 as entering the correlation function")
    compressed_array_summary(map_nbar2[freqlist, :, :],
                    "nbar 2 N^-1 as entering the correlation function")

    cross_pair = fs.MapPair(map_opt1, map_opt2, map_nbar1, map_nbar2,
                           freqlist)

    if params['subtract_mean']:
        cross_pair.subtract_weighted_mean()

    (corr, counts) = cross_pair.correlate(params['lags'],
                            speedup=params['speedup'])

    corr_shelve = shelve.open(params['output_shelve_file'])
    corr_shelve["corr"] = corr
    corr_shelve["counts"] = counts
    corr_shelve["freq_axis"] = map_opt1.get_axis('freq')
    corr_shelve["params"] = params
    corr_shelve.close()


# For running this module from the command line
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        initfilename = str(sys.argv[1])
        wigglez_correlation(initfilename)
    else:
        print 'Maximum one argument, a parameter file name.'
