import numpy as np
import scipy as sp
from core import algebra 
import operator
import matplotlib
matplotlib.use('Agg')
from correlate import freq_slices as fs
from map import optical_catalog as oc 
import correlation_plots as cp
from kiyopy import parse_ini
import shelve
from map import beam
# TODO: this seemed to be necessary on the Sunnyvale compute nodes because it
# was clobbing the python path?
#import sys, site
#site.addsitedir('/home/eswitzer/local/lib/')
#site.addsitedir('/home/eswitzer/local/lib/python2.6/site-packages/')

# TODO: delete /mnt/raid-project/gmrt/eswitzer/wiggleZ/maps/sec_A_15hr_41-73_clean_map_I.npy
params_init = {
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
               'threading': False
               }
prefix = 'fs_'

# TODO: move this somewhere for utilities
def compressed_array_summary(array, name, axes=[1,2], extras=False):
    """print various summaries of arrays compressed along specified axes""" 

    print "-"*80
    print "array property summary for "+name+":"
    arrayshape = array.shape 
    subarray_size = reduce(operator.mul,  [ arrayshape[i] for i in axes ])
    nanarray = np.isnan(array)
    infarray = np.isinf(array)
    zeroarray = (array == 0.)
    negarray = (array < 0.)

    total_nan = nanarray.sum()
    if total_nan != 0:
        print "nans:"
        nan_count = np.apply_over_axes(np.sum, nanarray, axes)
        print nan_count.flatten()
        print "has nan: "+repr(np.where( nan_count.flatten() != 0))
        print "total nans: "+repr(total_nan)
        print "-"*80
    else:
        print "There are no nans"

    total_inf = infarray.sum()
    if total_inf != 0:
        print "infs:"
        inf_count = np.apply_over_axes(np.sum, infarray, axes)
        print inf_count.flatten()
        print "has inf: "+repr(np.where( inf_count.flatten() != 0))
        print "total infs: "+repr(total_inf)
    else:
        print "There are no infs"

    total_zero = zeroarray.sum()
    if total_zero != 0:
        print "zeros:"
        zero_count = np.apply_over_axes(np.sum, zeroarray, axes)
        print zero_count.flatten()
        print "all zero: "+repr(np.where(zero_count.flatten() == subarray_size)) 
        #print "total zeros: "+repr(total_zero)
    else:
        print "There are no zeros"

    total_neg = negarray.sum()
    if total_neg != 0:
        print "negatives:"
        neg_count = np.apply_over_axes(np.sum, negarray, axes)
        print neg_count.flatten()
        #print np.where( neg_count.flatten() != 0)
        #print "total negatives: "+repr(total_zero)
    else:
        print "There are no negative entries"

    if extras:
        sum_nu = np.apply_over_axes(np.sum, array, axes)
        min_nu = np.apply_over_axes(np.min, array, axes)
        max_nu = np.apply_over_axes(np.max, array, axes)
        print sum_nu.flatten()
        print min_nu.flatten()
        print max_nu.flatten()

    print ""

def wigglez_correlation(init_filename):
    np.set_printoptions(threshold=np.nan)
    
    print "starting wigglez correlation run with file: "+init_filename
    params = parse_ini.parse(init_filename, params_init, prefix=prefix, feedback=10)

    radio_file = params['radio_root']+params['radio_data_file']
    noiseinv_file = params['radio_root']+params['radio_noiseinv_file']
    optical_file = params['optical_root']+params['optical_data_file']
    optical_selection_file = params['optical_root']+params['optical_selection_file']
    
    Map_radio = algebra.make_vect(algebra.load(radio_file))
    # TODO: why are the N^-1 matrices negative?
    Noiseinv_radio = np.abs(algebra.make_vect(algebra.load(noiseinv_file)))
    Map_opt = algebra.make_vect(algebra.load(optical_file))
    Map_nbar = algebra.make_vect(algebra.load(optical_selection_file))

    compressed_array_summary(Map_opt, "opt map as loaded") 
    compressed_array_summary(Map_nbar, "nbar map as loaded") 
    compressed_array_summary(Map_radio, "radio map as loaded") 
    compressed_array_summary(Noiseinv_radio, "radio N^-1 as loaded") 

    Noiseinv_radio[Noiseinv_radio<1.e-20] = 0.
    nan_array = np.isnan(Noiseinv_radio)
    Noiseinv_radio[nan_array] = 0.
    inf_array = np.isinf(Noiseinv_radio)
    Noiseinv_radio[inf_array] = 0.

    if params['convolve']:
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                              0.281176247549, 0.270856788455, 0.26745856078,
                              0.258910010848, 0.249188429031])
        # also consider the beam model where everything is degraded to the lowest
        # resolution; e.g. in a MapPair
        beam_data_degrade = sp.zeros(8)+max(beam_data)  
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                          dtype=float)
        freq_data *= 1.e6
        b = beam.GaussianBeam(beam_data_degrade, freq_data)
        # Convolve the optical data to the lowest radio resolution
        print "convolving the map"
        Map_opt = b.apply(Map_opt) 
        compressed_array_summary(Map_opt, "opt after convolution")  

        # how should the covariance be convolved? nbar is N^-1; convolve B N B^T?,
        # or just convolve nbar? use?:
        #Map_nbar[Map_nbar<1.e-30] = 1.e-30
        #Map_nbar = 1./Map_nbar
        #Map_nbar = b.apply(Map_nbar, cval=1.e30)
        #Map_nbar = 1./Map_nbar
        #Map_nbar[Map_nbar<1.e-20] = 0
        print "convolving the covariance/selection"
        Map_nbar = b.apply(Map_nbar)
        compressed_array_summary(Map_nbar, "nbar after convolution")

    # convert to delta-overdensity
    Map_opt = Map_opt/Map_nbar - 1
    #compressed_array_summary(Map_opt, "opt after conversion to delta") 

    # set the NaNs and infs to zero in data and weights
    nan_array = np.isnan(Map_opt)
    Map_opt[nan_array] = 0.
    Map_nbar[nan_array] = 0.
    inf_array = np.isinf(Map_opt)
    Map_opt[inf_array] = 0.
    Map_nbar[inf_array] = 0.
 
    freqlist = params['freq'] # full: range(Map_radio.shape[0])

    compressed_array_summary(Map_opt[freqlist,:,:], "opt map as entering the correlation function")
    compressed_array_summary(Map_nbar[freqlist,:,:], "nbar as entering the correlation function")
    compressed_array_summary(Map_radio[freqlist,:,:], "radio map as entering the correlation function")
    compressed_array_summary(Noiseinv_radio[freqlist,:,:], "radio map N^-1 as entering the correlation function")

    Cross_pair = fs.MapPair(Map_opt, Map_radio, Map_nbar, Noiseinv_radio,
                           freqlist)

    if params['subtract_mean']:
        Cross_pair.subtract_weighted_mean()

    corr = Cross_pair.correlate(params['lags'], threading=params['threading'])

    corr_shelve = shelve.open(params['output_shelve_file'])
    corr_shelve["corr"] = corr
    corr_shelve["freq_axis"] = Map_radio.get_axis('freq')
    corr_shelve["params"] = params
    corr_shelve.close()


def plot_crosscorr(filename, outputname, collapse=False, printcorr=False,
                   title=None, coloraxis=None):
    print title+" to "+outputname+" from "+filename
    corr_shelve = shelve.open(filename)
    corr = corr_shelve["corr"]
    run_params = corr_shelve["params"]

    print "run parameters \n"+"-"*80
    for key in run_params:
        print key + ": "+ repr(run_params[key])
    
    if printcorr:
        np.set_printoptions(threshold=np.nan)
        print corr_shelve["corr"]

    corr[np.isnan(corr)] = 0.
    corr[np.isinf(corr)] = 0.

    plot_object = cp.CorrelationPlot(run_params["lags"], corr,
                                     run_params["freq"], corr_shelve["freq_axis"])

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


def make_plots():
    coloraxisA = np.linspace(-0.03, 0.03, 100, endpoint=True)
    #coloraxisA = np.around(coloraxisA*1000.)/1000.
    print "using color A axis "+repr(coloraxisA)
    coloraxisB = np.linspace(-0.2, 0.2, 100, endpoint=True)
    #coloraxisB = np.around(coloraxisB*1000.)/1000.
    print "using color B axis "+repr(coloraxisB)

    plot_crosscorr("radio_opt_xcorr_mapA.shelve",
                   "opt_x_radio_mapA_noconv_sparse.png", collapse=False,
                   title="GBT_A x WiggleZ, not PSF convolved, sparse selection", 
                   coloraxis = coloraxisA)
    plot_crosscorr("radio_opt_xcorr_mapA.shelve",
                   "opt_x_radio_mapA_noconv_sparse_coll.png", collapse=True,
                   title="GBT_A x WiggleZ, not PSF convolved, sparse selection")
    plot_crosscorr("radio_opt_xcorr_mapB.shelve",
                   "opt_x_radio_mapB_noconv_sparse.png", collapse=False,
                   title="GBT_B x WiggleZ, not PSF convolved, sparse selection",
                   coloraxis = coloraxisB)
    plot_crosscorr("radio_opt_xcorr_mapB.shelve",
                   "opt_x_radio_mapB_noconv_sparse_coll.png", collapse=True,
                   title="GBT_B x WiggleZ, not PSF convolved, sparse selection")

    plot_crosscorr("opt_x_radio_mapA_noconv_sep.shelve",
                   "opt_x_radio_mapA_noconv_sep.png", collapse=False,
                   title="GBT_A x WiggleZ, not PSF convolved, separable selection",
                   coloraxis = coloraxisA)
    plot_crosscorr("opt_x_radio_mapA_noconv_sep.shelve",
                   "opt_x_radio_mapA_noconv_sep_coll.png", collapse=True,
                   title="GBT_A x WiggleZ, not PSF convolved, separable selection")
    plot_crosscorr("opt_x_radio_mapB_noconv_sep.shelve",
                   "opt_x_radio_mapB_noconv_sep.png", collapse=False,
                   title="GBT_B x WiggleZ, not PSF convolved, separable selection",
                   coloraxis = coloraxisB)
    plot_crosscorr("opt_x_radio_mapB_noconv_sep.shelve",
                   "opt_x_radio_mapB_noconv_sep_coll.png", collapse=True,
                   title="GBT_B x WiggleZ, not PSF convolved, separable selection")

    plot_crosscorr("opt_x_radio_mapAnull_noconv_sep.shelve",
                   "opt_x_radio_mapAnull_noconv_sep.png", collapse=False,
                   title="GBT_A x WiggleZ_random, not PSF convolved, separable selection",
                   coloraxis = coloraxisA)
    plot_crosscorr("opt_x_radio_mapAnull_noconv_sep.shelve",
                   "opt_x_radio_mapAnull_noconv_sep_coll.png", collapse=True,
                   title="GBT_A x WiggleZ_random, not PSF convolved, separable selection")
    plot_crosscorr("opt_x_radio_mapBnull_noconv_sep.shelve",
                   "opt_x_radio_mapBnull_noconv_sep.png", collapse=False,
                   title="GBT_B x WiggleZ_random, not PSF convolved, separable selection",
                   coloraxis = coloraxisB)
    plot_crosscorr("opt_x_radio_mapBnull_noconv_sep.shelve",
                   "opt_x_radio_mapBnull_noconv_sep_coll.png", collapse=True,
                   title="GBT_B x WiggleZ_random, not PSF convolved, separable selection")

    plot_crosscorr("opt_x_radio_mapA_conv_sep.shelve",
                   "opt_x_radio_mapA_conv_sep.png", collapse=False,
                   title="GBT_A x WiggleZ, PSF convolve, separable selection",
                   coloraxis = coloraxisA)
    plot_crosscorr("opt_x_radio_mapA_conv_sep.shelve",
                   "opt_x_radio_mapA_conv_sep_coll.png", collapse=True,
                   title="GBT_A x WiggleZ, PSF convolve, separable selection")
    plot_crosscorr("opt_x_radio_mapB_conv_sep.shelve",
                   "opt_x_radio_mapB_conv_sep.png", collapse=False,
                   title="GBT_B x WiggleZ, PSF convolve, separable selection",
                   coloraxis = coloraxisB)
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
        print 'Maximum one argument, a parameter file name.'

    #coloraxisA = np.linspace(-2., 2., 100, endpoint=True)
    #plot_crosscorr("opt_x_radio_mapA_noconv_sep.shelve",
    #               "stupid.png", title="ok", collapse=False, coloraxis=coloraxisA )

