from plotting import summary_plot as splt

# note that notes are purely human-readable and the keys do not mean anything
rootdir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/"
run1_notes = {
    "runname": "run1",
    "machine": "sunnyvale",
    "selection_function": "not separable, 1000 catalogs",
    "radio_map": "A cleaned with D",
    "threading": "off",
    "error_1": "job did not complete because of 48-hour limit",
    "error_2": "filename has _sep even though separability was not assumed",
    "error_3": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    "todo": "file filenames"
    }
batch1_param = {
    "path": rootdir + "data_run1",
    "rand:list": {"prefix": "opt_x_radio_mapArand",
                  "suffix": "_noconv_sep",
                  "indices": [0, 1, 2, 3, 4, 6, 7, 15, 16, 17, \
                              18, 19, 20, 21, 42, 43, 44, 45, 59, \
                              60, 61, 62, 63, 64, 97, 98, 99],
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_mapA_noconv_sep",
    "notes": run1_notes
    }

run2_notes = {
    "runname": "run1",
    "machine": "sunnyvale",
    "selection_function": "not separable, 1000 catalogs",
    "radio_map": "A cleaned with D",
    "threading": "on",
    "note": "should be same as run1, but all 100 rand. catalogs finished",
    "error": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    }
batch2_param = {
    "path": rootdir + "data_run2",
    "rand:list": {"prefix": "opt_x_radio_mapArand",
                  "suffix": "_noconv",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_mapA_noconv",
    "notes": run2_notes
    }

run3_notes = {
    "runname": "run3",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "A cleaned with D",
    "threading": "on",
    "error": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    }
batch3_param = {
    "path": rootdir + "data_run3",
    "rand:list": {"prefix": "opt_x_radio_mapArand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_mapA_noconv_sep",
    "notes": run3_notes
    }

run4_notes = {
    "runname": "run4",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD",
    "threading": "on",
    "error": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    }
batch4_param = {
    "path": rootdir + "data_run4",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "notes": run4_notes
    }

run5_notes = {
    "runname": "run5",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "sec A of Kiyo's old 2-way split (signal seen)",
    "threading": "on",
    "note": "this run also includes xcorr with selection function",
    }
truncated = range(100)
truncated.remove(2)
batch5_param = {
    "path": rootdir + "data_run5",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": truncated,
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "selxcorr:file": "optsel_x_radio_combined_noconv_sep",
    "notes": run5_notes
    }

run6_notes = {
    "runname": "run6",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "threading": "on",
    }
batch6_param = {
    "path": rootdir + "data_run6",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "selxcorr:file": "optsel_x_radio_combined_noconv_sep",
    "notes": run6_notes
    }

run7_notes = {
    "runname": "run7",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "notes": "should match run6 exactly, but with new correlate()"
    }
batch7_param = {
    "path": rootdir + "data_run7",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    #"selxcorr:file": "optsel_x_radio_combined_noconv_sep",
    "notes": run7_notes
    }

run8_notes = {
    "runname": "mode_xcorr_run1",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 0-25 modes removed",
    "speedup": "on"
    }
batch8_param = {
    "path": rootdir + "mode_xcorr_run1",
    "mode:list": {"prefix": "opt_x_radio_combined_mode",
                  "suffix": "_noconv_sep",
                  "indices": range(26),
                  "indexfmt": "%d",
                  "id_prefix": "mode"},
    "nullmode:list": {"prefix": "opt_x_radio_combined_nullmode",
                  "suffix": "_noconv_sep",
                  "indices": range(26),
                  "indexfmt": "%d",
                  "id_prefix": "nullmode"},
    "notes": run8_notes
    }

run9_notes = {
    "runname": "test_69old_vs_69new",
    "machine": "prawn",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "new and old treatments up to session 69",
    "speedup": "on"
    }
batch9_param = {
    "path": rootdir + "test_69old_vs_69new",
    "reallyoldway:file": "opt_x_radio_69old_kiyo_noconv_sep",
    "oldway:file": "opt_x_radio_69old_noconv_sep",
    "newway:file": "opt_x_radio_69_noconv_sep",
    "notes": run9_notes
    }

run10_notes = {
    "runname": "run10",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "notes": "should match run7 but using counts as corr bin weights"
    }
batch10_param = {
    "path": rootdir + "xcorr_15modes",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "notes": run10_notes
    }

run11_notes = {
    "runname": "run11",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "convolution": "on",
    "notes": "same as 10, but with convolution of the optical overdensity"
    }
batch11_param = {
    "path": rootdir + "xcorr_withconv",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_conv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_conv_sep",
    "notes": run11_notes
    }

run12_notes = {
    "runname": "run12",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "meansubtract": "off",
    "notes": "same as 10, but without subtracting means"
    }
batch12_param = {
    "path": rootdir + "xcorr_nomean",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_nomean_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_nomean_sep",
    "notes": run12_notes
    }

# /mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/xcorr_15modes_12split/opt_x_radio_15CD_noconv_sep.shelve
run13_notes = {
    "runname": "xcorr_15modes_12split",
    "machine": "prawn",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "12 ABCD maps",
    "speedup": "on"
    }
batch13_param = {
    "path": rootdir + "xcorr_15modes_12split",
    "xcorr15mAB:file": "opt_x_radio_15AB_noconv_sep",
    "xcorr15mAC:file": "opt_x_radio_15AC_noconv_sep",
    "xcorr15mAD:file": "opt_x_radio_15AD_noconv_sep",
    "xcorr15mBA:file": "opt_x_radio_15BA_noconv_sep",
    "xcorr15mBC:file": "opt_x_radio_15BC_noconv_sep",
    "xcorr15mBD:file": "opt_x_radio_15BD_noconv_sep",
    "xcorr15mCA:file": "opt_x_radio_15CA_noconv_sep",
    "xcorr15mCB:file": "opt_x_radio_15CB_noconv_sep",
    "xcorr15mCD:file": "opt_x_radio_15CD_noconv_sep",
    "xcorr15mDA:file": "opt_x_radio_15DA_noconv_sep",
    "xcorr15mDB:file": "opt_x_radio_15DB_noconv_sep",
    "xcorr15mDC:file": "opt_x_radio_15DC_noconv_sep",
    "notes": run13_notes
    }

run14_notes = {
    "runname": "run14",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "speedup": "on",
    "meansubtract": "on",
    "notes": "optical autocorr, DD, DR, RR"
    }
batch14_param = {
    "path": rootdir + "opt_auto_corr",
    "rand:list": {"prefix": "opt_x_opt_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "RR"},
    "DR:list": {"prefix": "opt_x_opt_DR",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "DR"},
    "signal:file": "opt_x_opt_noconv_sep",
    "notes": run14_notes
    }

if __name__ == '__main__':
    #splt.process_batch_correlations("run1_correlations.shelve", batch1_param)
    #splt.process_batch_correlations("run2_correlations.shelve", batch2_param)
    #splt.process_batch_correlations("run3_correlations.shelve", batch3_param)
    #splt.process_batch_correlations("run4_correlations.shelve", batch4_param)
    #splt.process_batch_correlations("run5_correlations.shelve", batch5_param)
    #splt.process_batch_correlations("run6_correlations.shelve", batch6_param)
    #splt.process_batch_correlations("run7_correlations.shelve", batch7_param)
    #splt.process_batch_correlations("run5_correlations_newcorr.shelve", batch5_param)
    #splt.process_batch_correlations("run7_correlations_newcorr.shelve", batch7_param)
    #splt.process_batch_correlations("run8_correlations_modes.shelve", batch8_param)
    #splt.process_batch_correlations("run9_correlations_modes.shelve", batch9_param)
    #splt.process_batch_correlations("run10_correlations_modes.shelve", batch10_param)
    #splt.process_batch_correlations("run11_correlations_modes.shelve", batch11_param)
    #splt.process_batch_correlations("run12_correlations_modes.shelve", batch12_param)
    #splt.process_batch_correlations("run13_correlations.shelve", batch13_param)
    splt.process_batch_correlations("run14_correlations.shelve", batch14_param)

    #print splt.compare_corr(batch2_param, batch3_param)
    #print splt.compare_corr(batch1_param, batch2_param)
    #print splt.compare_corr(batch6_param, batch7_param)

    #splt.plot_batch_correlations("run1_correlations.shelve", batch1_param,
    #                        dir_prefix="plots/run1b/")
    #splt.plot_batch_correlations("run2_correlations.shelve", batch2_param,
    #                        dir_prefix="plots/run2/")
    #splt.plot_batch_correlations("run3_correlations.shelve", batch3_param,
    #                        dir_prefix="plots/run3/")
    #splt.plot_batch_correlations("run4_correlations.shelve", batch4_param,
    #                        dir_prefix="plots/run4/")
    #splt.plot_batch_correlations("run5_correlations.shelve", batch5_param,
    #                        dir_prefix="plots/run5/",
    #                        color_range=[-0.04, 0.04])
    #splt.plot_batch_correlations("run6_correlations.shelve", batch6_param,
    #                        dir_prefix="plots/run6/",
    #                        color_range=[-0.04, 0.04])
    #splt.plot_batch_correlations("run7_correlations.shelve", batch7_param,
    #                        dir_prefix="plots/run7/",
    #                        color_range=[-0.04, 0.04])
    #splt.plot_batch_correlations("run7_correlations_newcorr.shelve", batch7_param,
    #                        dir_prefix="plots/run7n/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run8_correlations_modes.shelve", batch8_param,
    #                        dir_prefix="plots/run8/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run9_correlations_modes.shelve", batch9_param,
    #                        dir_prefix="plots/run9/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run10_correlations_modes.shelve", batch10_param,
    #                        dir_prefix="plots/run10/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run11_correlations_modes.shelve", batch11_param,
    #                        dir_prefix="plots/run11/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run12_correlations_modes.shelve", batch12_param,
    #                        dir_prefix="plots/run12/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run13_correlations.shelve", batch13_param,
    #                        dir_prefix="plots/run13/",
    #                        color_range=[-10, 10])
    #splt.plot_batch_correlations("run14_correlations.shelve", batch14_param,
    #                        dir_prefix="plots/run14/",
    #                        color_range=[-10, 10])

    #splt.batch_correlations_statistics("run1_correlations.shelve", batch1_param)
    #splt.batch_correlations_statistics("run2_correlations.shelve", batch2_param)
    #splt.batch_correlations_statistics("run3_correlations.shelve", batch3_param)
    #splt.batch_correlations_statistics("run4_correlations.shelve", batch4_param)
    #splt.batch_correlations_statistics("run5_correlations.shelve", batch5_param)
    #splt.batch_correlations_statistics("run6_correlations.shelve", batch6_param)
    #splt.batch_correlations_statistics("run7_correlations.shelve", batch7_param)
    #splt.batch_correlations_statistics("run5_correlations_newcorr.shelve", batch5_param)
    #splt.batch_correlations_statistics("run7_correlations_newcorr.shelve", batch7_param)
    #splt.batch_correlations_statistics("run10_correlations_modes.shelve", batch10_param)
    #splt.batch_correlations_statistics("run11_correlations_modes.shelve", batch11_param)
    #splt.batch_correlations_statistics("run12_correlations_modes.shelve", batch12_param)
    splt.batch_correlations_statistics("run14_correlations.shelve",
                                        batch14_param, randtoken="RR")
    splt.batch_correlations_statistics("run14_correlations.shelve",
                                        batch14_param, randtoken="DR")
