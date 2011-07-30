from plotting import summary_plot as splt
from simulations import corr21cm

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

run23_notes = {
    "runname": "run23",
    "machine": "prawn",
    "speedup": "on",
    "notes": "autocorrelation with repaired correlation function; for proposal"
    }
batch23_param = {
    "path": rootdir + "new_autocorr_test",
    "pair:list": {"prefix": "radio_autocorr_15_modes_",
                  "suffix": "",
                  "indices": range(0,6),
                  "indexfmt": "%d",
                  "id_prefix": "pair"},
    "notes": run23_notes
    }

def make_xcorr_plotdata():
    # find the mean brightness used in simulations
    corrobj = corr21cm.Corr21cm()
    T_b_sim = corrobj.T_b(1420./800.-1)
    print T_b_sim

    #splt.process_batch_correlations(xcorr_real_param,
    #                multiplier=1., cross_power=True)
    #splt.process_batch_correlations(xcorr_sim_param, cross_power=True,
    #                                multiplier=1./(T_b_sim/1.e3))
    splt.process_batch_correlations(xcorr_sim_param, cross_power=True,
                                    multiplier=1./(T_b_sim))
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
    xcorr_signal = xcorr_data[5]
    xcorr_err = xcorr_data[4]
    xcorr_sim = xcorr_sim_data[3]
    xcorr_sim_err = xcorr_sim_data[4]
    # 0=0, 5=1, 10=2, 15=3, etc.
    compensation = compmode[3 ,:]

    # now convert from the 1e-3 Omega_HI in simulations to 0.5e-3
    xcorr_sim /= 2.
    xcorr_sim_err /= 2.

    for correntry in zip(lagl, lagc, lagr, xcorr_signal, xcorr_null, xcorr_err,
                      xcorr_sim, xcorr_sim_err, compensation):
        print "%5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g" % correntry

def make_corr_plotdata():
    print "write me"

if __name__ == '__main__':
    make_xcorr_plotdata()
