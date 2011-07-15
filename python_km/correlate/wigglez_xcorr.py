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
#import sys, site
#site.addsitedir('/home/eswitzer/local/lib/')
#site.addsitedir('/home/eswitzer/local/lib/python2.6/site-packages/')

# TODO: delete /mnt/raid-project/gmrt/eswitzer/wiggleZ/maps/sec_A_15hr_41-73_clean_map_I.npy

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

    radio_file = params['radio_root'] + params['radio_data_file']
    noiseinv_file = params['radio_root'] + params['radio_noiseinv_file']
    optical_file = params['optical_root'] + params['optical_data_file']
    optical_selection_file = params['optical_root'] + \
                             params['optical_selection_file']

    map_radio = algebra.make_vect(algebra.load(radio_file))
    # TODO: why are the N^-1 matrices negative?
    noiseinv_radio = np.abs(algebra.make_vect(algebra.load(noiseinv_file)))
    map_opt = algebra.make_vect(algebra.load(optical_file))
    map_nbar = algebra.make_vect(algebra.load(optical_selection_file))

    compressed_array_summary(map_opt, "opt map as loaded")
    compressed_array_summary(map_nbar, "nbar map as loaded")
    compressed_array_summary(map_radio, "radio map as loaded")
    compressed_array_summary(noiseinv_radio, "radio N^-1 as loaded")

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
        compressed_array_summary(map_opt, "opt after convolution")

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
        compressed_array_summary(map_nbar, "nbar after convolution")

    # convert to delta-overdensity
    map_opt = map_opt / map_nbar - 1
    #compressed_array_summary(map_opt, "opt after conversion to delta")

    # set the NaNs and infs to zero in data and weights
    nan_array = np.isnan(map_opt)
    map_opt[nan_array] = 0.
    map_nbar[nan_array] = 0.
    inf_array = np.isinf(map_opt)
    map_opt[inf_array] = 0.
    map_nbar[inf_array] = 0.

    freqlist = params['freq']  # full: range(map_radio.shape[0])

    compressed_array_summary(map_opt[freqlist, :, :],
                    "opt map as entering the correlation function")
    compressed_array_summary(map_nbar[freqlist, :, :],
                    "nbar as entering the correlation function")
    compressed_array_summary(map_radio[freqlist, :, :],
                    "radio map as entering the correlation function")
    compressed_array_summary(noiseinv_radio[freqlist, :, :],
                    "radio map N^-1 as entering the correlation function")

    cross_pair = fs.MapPair(map_opt, map_radio, map_nbar, noiseinv_radio,
                           freqlist)

    if params['subtract_mean']:
        cross_pair.subtract_weighted_mean()

    (corr, counts) = cross_pair.correlate(params['lags'],
                            speedup=params['speedup'])

    corr_shelve = shelve.open(params['output_shelve_file'])
    corr_shelve["corr"] = corr
    corr_shelve["freq_axis"] = map_radio.get_axis('freq')
    corr_shelve["params"] = params
    corr_shelve.close()


# TODO: move this into the plot correlation class and make the plotting class
# more of a correlation object data class
def plot_crosscorr(filename, outputname, collapse=False, printcorr=False,
                   title=None, coloraxis=None):
    """wrap the plot correlation class which reads correlation object shelve
    files"""
    print title + " to " + outputname + " from " + filename
    corr_shelve = shelve.open(filename)
    corr = corr_shelve["corr"]
    run_params = corr_shelve["params"]

    print "run parameters \n" + "-" * 80
    for key in run_params:
        print key + ": " + repr(run_params[key])

    if printcorr:
        np.set_printoptions(threshold=np.nan)
        print corr_shelve["corr"]

    corr[np.isnan(corr)] = 0.
    corr[np.isinf(corr)] = 0.

    plot_object = cp.CorrelationPlot(run_params["lags"], corr,
                                     run_params["freq"],
                                     corr_shelve["freq_axis"])

    if collapse:
        plot_object.plot_collapsed(outputname, norms=False,
                                   lag_inds=range(len(run_params["lags"])),
                                   cross_power=True, title=title)

    else:
        plot_object.plot_contour(outputname, norms=False,
                                 lag_inds=range(len(run_params["lags"])),
                                 cross_power=True, title=title,
                                 coloraxis=coloraxis)

    corr_shelve.close()


# TODO: make plots in a more automated way
def make_plots():
    """Driver to make a lot of plots (not production code)
    """
    coloraxis_a = np.linspace(-0.03, 0.03, 100, endpoint=True)
    #coloraxis_a = np.around(coloraxis_a * 1000.) / 1000.
    print "using color A axis " + repr(coloraxis_a)
    coloraxis_b = np.linspace(-0.2, 0.2, 100, endpoint=True)
    #coloraxis_b = np.around(coloraxis_b * 1000.) / 1000.
    print "using color B axis " + repr(coloraxis_b)

    plot_crosscorr("radio_opt_xcorr_mapA.shelve",
                   "opt_x_radio_mapA_noconv_sparse.png", collapse=False,
                   title="GBT_A x WiggleZ, not PSF convolved, sparse selection",
                   coloraxis=coloraxis_a)
    plot_crosscorr("radio_opt_xcorr_mapA.shelve",
                   "opt_x_radio_mapA_noconv_sparse_coll.png", collapse=True,
                   title="GBT_A x WiggleZ, not PSF convolved, sparse selection")
    plot_crosscorr("radio_opt_xcorr_mapB.shelve",
                   "opt_x_radio_mapB_noconv_sparse.png", collapse=False,
                   title="GBT_B x WiggleZ, not PSF convolved, sparse selection",
                   coloraxis=coloraxis_b)
    plot_crosscorr("radio_opt_xcorr_mapB.shelve",
                   "opt_x_radio_mapB_noconv_sparse_coll.png", collapse=True,
                   title="GBT_B x WiggleZ, not PSF convolved, sparse selection")

    plot_crosscorr("opt_x_radio_mapA_noconv_sep.shelve",
                   "opt_x_radio_mapA_noconv_sep.png", collapse=False,
                   title="GBT_A x WiggleZ, not PSF convolved, \
                          separable selection",
                   coloraxis=coloraxis_a)
    plot_crosscorr("opt_x_radio_mapA_noconv_sep.shelve",
                   "opt_x_radio_mapA_noconv_sep_coll.png", collapse=True,
                   title="GBT_A x WiggleZ, not PSF convolved, \
                          separable selection")
    plot_crosscorr("opt_x_radio_mapB_noconv_sep.shelve",
                   "opt_x_radio_mapB_noconv_sep.png", collapse=False,
                   title="GBT_B x WiggleZ, not PSF convolved, \
                          separable selection",
                   coloraxis=coloraxis_b)
    plot_crosscorr("opt_x_radio_mapB_noconv_sep.shelve",
                   "opt_x_radio_mapB_noconv_sep_coll.png", collapse=True,
                   title="GBT_B x WiggleZ, not PSF convolved, \
                          separable selection")

    plot_crosscorr("opt_x_radio_mapAnull_noconv_sep.shelve",
                   "opt_x_radio_mapAnull_noconv_sep.png", collapse=False,
                   title="GBT_A x WiggleZ_random, not PSF convolved, \
                          separable selection",
                   coloraxis=coloraxis_a)
    plot_crosscorr("opt_x_radio_mapAnull_noconv_sep.shelve",
                   "opt_x_radio_mapAnull_noconv_sep_coll.png", collapse=True,
                   title="GBT_A x WiggleZ_random, not PSF convolved, \
                          separable selection")
    plot_crosscorr("opt_x_radio_mapBnull_noconv_sep.shelve",
                   "opt_x_radio_mapBnull_noconv_sep.png", collapse=False,
                   title="GBT_B x WiggleZ_random, not PSF convolved, \
                          separable selection",
                   coloraxis=coloraxis_b)
    plot_crosscorr("opt_x_radio_mapBnull_noconv_sep.shelve",
                   "opt_x_radio_mapBnull_noconv_sep_coll.png", collapse=True,
                   title="GBT_B x WiggleZ_random, not PSF convolved, \
                          separable selection")

    plot_crosscorr("opt_x_radio_mapA_conv_sep.shelve",
                   "opt_x_radio_mapA_conv_sep.png", collapse=False,
                   title="GBT_A x WiggleZ, PSF convolve, separable selection",
                   coloraxis=coloraxis_a)
    plot_crosscorr("opt_x_radio_mapA_conv_sep.shelve",
                   "opt_x_radio_mapA_conv_sep_coll.png", collapse=True,
                   title="GBT_A x WiggleZ, PSF convolve, separable selection")
    plot_crosscorr("opt_x_radio_mapB_conv_sep.shelve",
                   "opt_x_radio_mapB_conv_sep.png", collapse=False,
                   title="GBT_B x WiggleZ, PSF convolve, separable selection",
                   coloraxis=coloraxis_b)
    plot_crosscorr("opt_x_radio_mapB_conv_sep.shelve",
                   "opt_x_radio_mapB_conv_sep_coll.png", collapse=True,
                   title="GBT_B x WiggleZ, PSF convolve, separable selection")

# For running this module from the command line
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        initfilename = str(sys.argv[1])
        wigglez_correlation(initfilename)
    else:
        #print 'Maximum one argument, a parameter file name.'
        print "no arguments given, just doing some user-specified nonsense \
               instead"
        coloraxis_a = np.linspace(-0.2, 0.2, 100, endpoint=True)
        plot_crosscorr("data/opt_x_radio_mapA_noconv_sep.shelve",
                       "signal_xcorr.png", title="opt_x_radio_mapA_noconv", collapse=True,
                       coloraxis=coloraxis_a )
