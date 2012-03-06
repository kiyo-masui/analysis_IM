import sys
from plotting import summary_plot as splt
from simulations import corr21cm
# to unpickle freq_slices class instances, they needs the class definitions
# imported at the top level: (TODO: swtich to more transparent databases)
from correlate.freq_slices import NewSlices
from correlate.freq_slices import MapPair
import multiprocessing
import shelve
import utils.misc as utils
import numpy as np

# note that notes are purely human-readable and the keys do not mean anything
rootdir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/"

xcorr_real_notes = {
    "runname": "xcorr_real",
    "machine": "scinet",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "this has the wigglez xcorr measurement with the new corr function",
    "notes2": "it also has the 100 random optical xcorr"
    }
xcorr_real_param = {
    "path": rootdir + "xcorr_newcorr",
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "notes": xcorr_real_notes
    }

xcorr_sim_notes = {
    "runname": "xcorr_sim",
    "machine": "scinet",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "cross-correlation of (real map (15 modes removed) + sim) x sim",
    "notes2": "this uses radio N^-1 for sim1, nbar for sim2; simulate xcorr"
    }
xcorr_sim_param = {
    "path": rootdir + "xcorr_newcorr_sim_variants",
    "randsim:list": {"prefix": "radioplussim_x_optsim_",
                  "suffix": "_est_noconvopt",
                  "indices": range(1,101),
                  "indexfmt": "%d",
                  "id_prefix": "rand"},
    "notes": xcorr_sim_notes
    }
#    "path": rootdir + "xcorr_newcorr",
#                  "suffix": "",

xcorr_loss_notes = {
    "runname": "xcorr_loss",
    "machine": "prawn",
    "selection_function": "non-separable, 200,000 catalogs",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "mode loss test relevant for cross-correlation",
    }
xcorr_loss_param = {
    "path": rootdir + "xcorr_modeloss",
    "modedata:list": {"prefix": "opt_x_radio_combined_mode",
                  "suffix": "_noconv_est",
                  "indices": range(0, 55, 5),
                  "indexfmt": "%d",
                  "id_prefix": "modedata"},
    "notes": xcorr_loss_notes
    }

xcorr_loss_sim_notes = {
    "runname": "xcorr_loss_sim",
    "machine": "prawn",
    "selection_function": "non-separable, 200,000 catalogs",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "mode loss test relevant for cross-correlation",
    "notes2": "sim map with N modes removed x sim map with zero modes removed",
    "notes3": "use radio noise as weight on mode-removed side",
    "notes4": "use nbar weight on signal side"
    }
xcorr_loss_sim_param = {
    "path": rootdir + "xcorr_modeloss_sim",
    "modedata:list": {"prefix": "simmodeloss_x_optsim_",
                  "suffix": "_est",
                  "indices": range(0, 55, 5),
                  "indexfmt": "%d",
                  "id_prefix": "modedata"},
    "notes": xcorr_loss_sim_notes
    }

xcorr_variants_notes = {
    "runname": "xcorr_variants",
    "machine": "prawn",
    "speedup": "on",
    "notes": "cases with non-separable selection: convolving, taking out means"
    }
xcorr_variants_param = {
    "path": rootdir + "xcorr_newcorr_variants",
    "signal:file": "opt_x_radio_combined_noconv_est",
    "signalconv:file": "opt_x_radio_combined_conv_est",
    "signalnomean:file": "opt_x_radio_combined_noconv_nomean_est",
    "notes": xcorr_variants_notes
    }

autocorr_sim_notes = {
    "runname": "autocorr_sim",
    "machine": "scinet",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "cross-correlation of (real map (15 modes removed) + sim) x sim",
    "notes2": "this uses radio N^-1 for sim1, nbar for sim2; simulate xcorr"
    }
autocorr_sim_param = {
    "path": rootdir + "radio_x_radio_newcorr",
    "randsim:list": {"prefix": "radiosim_x_radiosim_",
                  "suffix": "",
                  "indices": range(1,101),
                  "indexfmt": "%d",
                  "id_prefix": "rand"},
    "notes": autocorr_sim_notes
    }

def make_xcorr_plotdata():
    # find the mean brightness used in simulations
    corrobj = corr21cm.Corr21cm()
    T_b_sim = corrobj.T_b(1420./800.-1)
    print T_b_sim

    #splt.process_batch_correlations(xcorr_real_param,
    #                multiplier=1., cross_power=True)
    ##splt.process_batch_correlations(xcorr_sim_param, cross_power=True,
    ##                                multiplier=1./(T_b_sim/1.e3))
    #splt.process_batch_correlations(xcorr_sim_param, cross_power=True,
    #                                multiplier=1./(T_b_sim))
    #splt.process_batch_correlations(xcorr_loss_sim_param,
    #                multiplier=1./T_b_sim*1.e-3, cross_power=True)
    #splt.process_batch_correlations(xcorr_variants_param,
    #                multiplier=1., cross_power=True)

    #splt.plot_batch_correlations(xcorr_real_param,
    #                        dir_prefix="plots/xcorr_real/",
    #                        color_range=[-0.04, 0.04], cross_power=True)
    #splt.plot_batch_correlations(xcorr_sim_param,
    #                        dir_prefix="plots/xcorr_sim/",
    #                        color_range=[-0.04, 0.04], cross_power=True)
    #splt.plot_batch_correlations(xcorr_loss_param,
    #                        dir_prefix="plots/xcorr_loss/",
    #                        color_range=[-0.04, 0.04], cross_power=True)
    #splt.plot_batch_correlations(xcorr_loss_sim_param,
    #                        dir_prefix="plots/xcorr_loss_sim/",
    #                        color_range=[-0.04, 0.04], cross_power=True)
    #splt.plot_batch_correlations(xcorr_variants_param,
    #                        dir_prefix="plots/xcorr_variants/",
    #                        color_range=[-0.04, 0.04], cross_power=True)

    # find the real xcorr signal with errors
    xcorr_data = splt.batch_correlations_statistics(xcorr_real_param,
                                       randtoken="rand",
                                       include_signal=True)
    # find the simulated xcorr signal with errors
    xcorr_sim_data = splt.batch_correlations_statistics(xcorr_sim_param,
                                       randtoken="rand",
                                       include_signal=False)

    # find the compensation function from simulations
    compmode = splt.batch_compensation_function(xcorr_loss_sim_param)

    lagl = xcorr_data[0]
    lagc = xcorr_data[1]
    lagr = xcorr_data[2]
    xcorr_null = xcorr_data[3]
    xcorr_signal = xcorr_data[6]
    xcorr_cov = xcorr_data[5]
    xcorr_err = xcorr_data[4]
    xcorr_sim = xcorr_sim_data[3]
    xcorr_sim_err = xcorr_sim_data[4]
    xcorr_sim_cov = xcorr_sim_data[5]
    # 0=0, 5=1, 10=2, 15=3, etc.
    compensation = compmode[3 ,:]

    (amp, amp_err) = utils.ampfit(xcorr_signal, xcorr_cov + xcorr_sim_cov, xcorr_sim)
    print amp, amp_err

    # now convert from the 1e-3 Omega_HI in simulations to 0.5e-3
    #xcorr_sim /= 2.
    #xcorr_sim_err /= 2.
    xcorr_sim *= amp
    xcorr_sim_err *= amp

    for correntry in zip(lagl, lagc, lagr, xcorr_signal, xcorr_null, xcorr_err,
                      xcorr_sim, xcorr_sim_err, compensation):
        print "%5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g" % correntry


def wrap_pickle_loader(runitem):
    (run_id, run_file) = runitem
    print run_file
    return splt.make_autocorr(run_file, identifier=run_id)


def make_corr_plotdata(filename):
    master = shelve.open(filename)

    autocorr_root = "/mnt/raid-project/gmrt/kiyo/wiggleZ/scinet_runs/"
    autocorr_list = [("autocorr"+repr(mode_num)+"_mode", autocorr_root +
                      "auto_corr_"+repr(mode_num)+"_modes/New_Slices_object.pkl")
                      for mode_num in range(0, 55, 5)]
    autocorr_sim_list = [("sim"+repr(sim_num)+"_mode"+repr(mode_num),
                         autocorr_root + "sim_" + repr(sim_num) +
                         "_corr_" + repr(mode_num) + "_modes/New_Slices_object.pkl")
                         for mode_num in range(0, 55, 5)
                         for sim_num in range(1,10)]

    #print autocorr_list
    #print autocorr_sim_list

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    results = pool.map(wrap_pickle_loader, autocorr_sim_list)

    for item in results:
        master[item[0]] = item[1]

    results = pool.map(wrap_pickle_loader, autocorr_list)

    for item in results:
        master[item[0]] = item[1]

    master.close()


def make_corr_plots(filename):
    corrobj = corr21cm.Corr21cm()
    T_b_sim = corrobj.T_b(1420./800.-1)
    print T_b_sim
    master = shelve.open(filename, "r")
    nlag = 15
    nmode = 11
    nsim = 9

    # load the autocorr simulations --------------------------------------------
    #splt.process_batch_correlations(autocorr_sim_param, cross_power=False,
    #                                multiplier=1.)

    autocorr_sim_data = splt.batch_correlations_statistics(autocorr_sim_param,
                                       randtoken="rand",
                                       include_signal=False)

    lagl = autocorr_sim_data[0]
    lagc = autocorr_sim_data[1]
    lagr = autocorr_sim_data[2]
    autocorr_sim = autocorr_sim_data[3]
    autocorr_sim_err = autocorr_sim_data[4]
    autocorr_sim_cov = autocorr_sim_data[5]

    # now convert from the 1e-3 Omega_HI in simulations to 0.5e-3
    #autocorr_sim /= 2.
    #autocorr_sim_err /= 2.
    # calibrate to xcorr: 0.56126646473 0.229801269703
    autocorr_sim_low = autocorr_sim*(0.56-0.23)
    autocorr_sim_high = autocorr_sim*(0.56+0.23)
    autocorr_sim_center = autocorr_sim*0.56

    for correntry in zip(lagl, lagc, lagr,
                    autocorr_sim_low, autocorr_sim_center, autocorr_sim_high):
        print "%5.3g %5.3g %5.3g %5.3g %5.3g %5.3g" % correntry

    sys.exit()

    # treat the autocorr and loss simulations----------------------------------
    lags = np.zeros((3, nlag))
    autocorr = np.zeros((nmode, nlag))
    autocorr_err = np.zeros((nmode, nlag))
    autocorr_sim = np.zeros((nsim, nmode, nlag))
    autocorr_sim_avg = np.zeros((nmode, nlag))
    compensated_autocorr = np.zeros((nmode, nlag))
    compensated_autocorr_err = np.zeros((nmode, nlag))
    mode_compensation = np.zeros((nmode, nlag))
    for mode_index in range(0, nmode):
        mode_num = mode_index*5
        autocorr_name = "autocorr" + repr(mode_num) + "_mode"
        entry = master[autocorr_name]
        autocorr[mode_index, :] = entry["corr1D"]
        autocorr_err[mode_index, :] = entry["corr1D_std"]
        lags[0, :] = entry["x_axis"][0]
        lags[1, :] = entry["x_axis"][1]
        lags[2, :] = entry["x_axis"][2]
        for sim_index in range(0, nsim):
            sim_name = "sim" + repr(sim_index+1) + "_mode" + repr(mode_num)
            entry = master[sim_name]
            autocorr_sim[sim_index, mode_index, :] = entry["corr1D"]
        autocorr_sim_avg = np.mean(autocorr_sim, axis=0)/1000.

    zero_modes = autocorr_sim_avg[0,:]
    mode_compensation = autocorr_sim_avg/zero_modes[None, :]
    compensated_autocorr = autocorr/mode_compensation
    compensated_autocorr_err = autocorr_err/mode_compensation
    print mode_compensation

    for lagind in range(0, nlag):
        lag_bounds = splt.fancy_vector(lags[:,lagind], "%5.2g")
        modes = splt.fancy_vector(autocorr[:,lagind], "%5.2g")
        modes_err = splt.fancy_vector(autocorr_err[:,lagind], "%5.2g")
        print lag_bounds + modes + modes_err

    for lagind in range(0, nlag):
        lag_bounds = splt.fancy_vector(lags[:,lagind], "%5.2g")
        modes = splt.fancy_vector(autocorr_sim_avg[:,lagind], "%5.2g")
        print lag_bounds + modes

    for lagind in range(0, nlag):
        lag_bounds = splt.fancy_vector(lags[:,lagind], "%5.2g")
        modes = splt.fancy_vector(compensated_autocorr[:,lagind], "%5.2g")
        modes_err = splt.fancy_vector(compensated_autocorr_err[:,lagind], "%5.2g")
        print lag_bounds + modes + modes_err

    master.close()

if __name__ == '__main__':
    #make_xcorr_plotdata()
    #make_corr_plotdata("autocorr_data.shelve")
    make_corr_plots("autocorr_data.shelve")
