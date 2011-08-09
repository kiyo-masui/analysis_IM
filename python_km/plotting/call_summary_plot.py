from plotting import summary_plot as splt
from simulations import corr21cm

# note that notes are purely human-readable and the keys do not mean anything
rootdir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/"

# /mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/simmodeloss_x_optsim/simmodeloss_x_optsim_0.shelve
xcorr_loss_notes = {
    "runname": "xcorr_loss",
    "machine": "prawn",
    "selection_function": "non-separable, 200,000 catalogs",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "mode loss test relevant for cross-correlation",
    "notes2": "sim map with N modes removed x sim map with zero modes removed",
    "notes3": "use radio noise as weight on mode-removed side",
    "notes4": "use nbar weight on signal side"
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

if __name__ == '__main__':
    # find the mean brightness used in simulations
    corrobj = corr21cm.Corr21cm()
    T_b_sim = corrobj.T_b(1420./800.-1)
    print T_b_sim

    #splt.process_batch_correlations(xcorr_loss_param,
    #                multiplier=1./T_b_sim*1.e-3, cross_power=True)
    splt.process_batch_correlations(xcorr_loss_param,
                    multiplier=1./T_b_sim, cross_power=True)

    splt.plot_batch_correlations(xcorr_loss_param,
                            dir_prefix="plots/xcorr_loss/",
                            color_range=[-0.04, 0.04], cross_power=True)

    splt.batch_compensation_function(xcorr_loss_param)

