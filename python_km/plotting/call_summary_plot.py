from plotting import summary_plot as splt
from simulations import corr21cm

# note that notes are purely human-readable and the keys do not mean anything
rootdir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/"

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

xcorr_sim_notes = {
    "runname": "xcorr_sim",
    "machine": "scinet",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "cross-correlation of (real map (15 modes removed) + sim) x sim",
    "notes2": "this uses radio N^-1 for sim1, nbar for sim2; simulate xcorr"
    }
xcorr_sim_param = {
    "path": rootdir + "xcorr_newcorr",
    "randsim:list": {"prefix": "radioplussim_x_optsim_",
                  "suffix": "",
                  "indices": range(1,101),
                  "indexfmt": "%d",
                  "id_prefix": "rand"},
    "notes": xcorr_sim_notes
    }


if __name__ == '__main__':
    # find the mean brightness used in simulations
    corrobj = corr21cm.Corr21cm()
    T_b_sim = corrobj.T_b(1420./800.-1)
    print T_b_sim

    #splt.process_batch_correlations(xcorr_loss_param,
    #                multiplier=1., cross_power=True)
    #splt.process_batch_correlations(xcorr_loss_sim_param,
    #                multiplier=1./T_b_sim*1.e-3, cross_power=True)
    splt.process_batch_correlations(xcorr_sim_param, cross_power=True,
                                    multiplier=1./(T_b_sim/1.e3))

    #splt.plot_batch_correlations(xcorr_loss_param,
    #                        dir_prefix="plots/xcorr_loss/",
    #                        color_range=[-0.04, 0.04], cross_power=True)
    #splt.plot_batch_correlations(xcorr_loss_sim_param,
    #                        dir_prefix="plots/xcorr_loss_sim/",
    #                        color_range=[-0.04, 0.04], cross_power=True)
    splt.plot_batch_correlations(xcorr_sim_param,
                            dir_prefix="plots/xcorr_sim/",
                            color_range=[-0.04, 0.04], cross_power=True)

    splt.batch_correlations_statistics(xcorr_sim_param, randtoken="rand",
                                       include_signal=False)

    #splt.batch_compensation_function(xcorr_loss_sim_param)

