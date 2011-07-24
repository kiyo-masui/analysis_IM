from plotting import summary_plot as splt
from simulations import corr21cm

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

# /mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/simmodeloss_x_optsim/simmodeloss_x_optsim_0.shelve
run15_notes = {
    "runname": "run15",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "speedup": "on",
    "meansubtract": "on",
    "notes1": "mode loss test relevant for cross-correlation",
    "notes2": "sim map with N modes removed x sim map with zero modes removed",
    "notes3": "use radio noise as weight on mode-removed side",
    "notes4": "use nbar weight on signal side"
    }
batch15_param = {
    "path": rootdir + "simmodeloss_x_optsim",
    "mode:list": {"prefix": "simmodeloss_x_optsim_",
                  "suffix": "",
                  "indices": range(0,26),
                  "indexfmt": "%d",
                  "id_prefix": "mode"},
    "notes": run15_notes
    }


# use this for repairing files
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


if __name__ == '__main__':
    # find the mean brightness used in simulations
    corrobj = corr21cm.Corr21cm()
    T_b_sim = corrobj.T_b(1420./800.-1)

    #splt.repair_shelve_files(batch15_param, "sim_xloss_correlate_mode",
    #                         params_default, prefix)
    #print splt.compare_corr(batch6_param, batch7_param)

    splt.process_batch_correlations(batch10_param)
    splt.process_batch_correlations(batch15_param, multiplier=1./T_b_sim*1.e-3)

    splt.plot_batch_correlations(batch10_param,
                            dir_prefix="plots/run10/",
                            color_range=[-10, 10], cross_power=True)
    splt.plot_batch_correlations(batch15_param,
                            dir_prefix="plots/run15/",
                            color_range=[-10, 10], cross_power=True)

    #splt.batch_correlations_statistics(batch14_param, randtoken="DR")

    #splt.batch_compensation_function(batch15_param)
