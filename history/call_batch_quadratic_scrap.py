def find_crossbeam_trans():
    # TODO: should not use the noconv weight here... need to separate these
    sim_15hr = gather_batch_genericsim_run("sim_15hr", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr")

    sim_15hr_beam = gather_batch_genericsim_run("sim_15hr_beam", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr_beam")

    sim_15hr_beam_mean = gather_batch_genericsim_run("sim_15hr_beam_meansub",
                         "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo",
                         "sim_wigglezxGBT15hr_beam_meansub")

    sim_15hr_beam_meanconv = gather_batch_genericsim_run("sim_15hr_beam_meansubconv",
                         "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo",
                         "sim_wigglezxGBT15hr_beam_meansubconv")

    trans_beam = tf.calculate_2d_transfer_function(sim_15hr_beam[2],
                                                sim_15hr[2],
                                                "crosstrans_beam")

    trans_beam_mean = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_mean[2],
                                                sim_15hr[2],
                                                "crosstrans_beam_mean")

    trans_beam_meanconv = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_meanconv[2],
                                                sim_15hr[2],
                                                "crosstrans_beam_meanconv")

    # TODO: check that it makes sense by reconstructing the simulated power without
    # these effects:
    #sim_15hr_corrected = gather_batch_sim_run("sim_15hr_beam",
    #                                     "sim_15hr_beam_meanconv_corrected",
    #                                     degrade_resolution=True,
    #                                     subtract_mean=True,
    #                                     transfer=trans_beam_meanconv)

    return (trans_beam, trans_beam_mean, trans_beam_meanconv)


def wigglez_x_GBT():
    gather_batch_genericsim_run("simvel_15hr", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "simvel_wigglezxGBT15hr")

    gather_batch_genericsim_run("simvel_15hr_beam", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "simvel_wigglezxGBT15hr_beam")

    gather_batch_genericsim_run("sim_15hr", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr")

    gather_batch_genericsim_run("sim_15hr_beam", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr_beam")


if __name__ == '__main__':
    #find_crossbeam_trans()
    wigglez_x_GBT()
    sys.exit()


#-----------------------------BARF---------------

def process_wigglez_fluxcal():
    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_nocomp_fluxcal",
                               "GBT_15hr_map_combined_cleaned_fluxcal",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")

    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_nocomp_noconv_fluxcal",
                               "GBT_15hr_map_combined_cleaned_noconv_fluxcal",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")


def process_wigglez(noconv=False, fluxcal=False):

    option_flag = ""
    if noconv:
        option_flag += "_noconv"
    if fluxcal:
        option_flag += "_fluxcal"

    # TODO: REVERT ME?
    (crosstrans_beam, crosstrans_beam_mean, crosstrans_beam_meanconv) = cwt.find_crossbeam_trans()

    # the "theory curve"
    # TODO: run this in the case with the common res convolution (impact in
    # weight function)
    #theory_curve = cwt.gather_batch_genericsim_run("sim_15hr", "sim_15hr_delta",
    #                                "GBT_15hr_map_combined_cleaned%s_0mode_weight" % option_flag,
    #                                       "WiggleZ_15hr_montecarlo",
    #                                       "sim_wigglezxGBT15hr_theory_curve")
    theory_curve = cwt.gather_batch_genericsim_run("sim_15hr", "sim_15hr_delta",
                                    "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                                           "WiggleZ_15hr_montecarlo",
                                           "sim_wigglezxGBT15hr_theory_curve")
    # pick the version from 3D->2D->1D
    theory_curve = theory_curve[1]
    (theory_curve, std_theory, corrmat_theory, covmat_theory) = pe.agg_stat_1d_pwrspec(theory_curve)

    wigglez_transfer_functions_withbeam = cwt.gather_batch_gbtxwigglez_trans_run(
                                    "WiggleZxGBT_withbeam%s" % option_flag,
                                    "sim_15hr_combined_cleaned%s" % option_flag,
                                    "sim_15hr_delta",
                                    "sim_15hr",
                                    "GBT_15hr_map_combined_cleaned%s" % option_flag,
                                    "WiggleZ_15hr_montecarlo")

    wigglez_transfer_functions_nobeam = cwt.gather_batch_gbtxwigglez_trans_run(
                                    "WiggleZxGBT_nobeam%s" % option_flag,
                                    "sim_15hr_combined_cleaned%s" % option_flag,
                                    "sim_15hr_delta",
                                    "sim_15hr_beam",
                                    "GBT_15hr_map_combined_cleaned%s" % option_flag,
                                    "WiggleZ_15hr_montecarlo")

    # now find the cross-power spectra
    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_nocomp%s" % option_flag,
                               "GBT_15hr_map_combined_cleaned%s" % option_flag,
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")

    # adding 2D beam compensation
    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_2dbeamcomp%s" % option_flag,
                               "GBT_15hr_map_combined_cleaned%s" % option_flag,
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo",
                               beam_transfer=crosstrans_beam)

    # adding 1D mode compensation
    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_2dbeam1dmodecomp%s" % option_flag,
                               "GBT_15hr_map_combined_cleaned%s" % option_flag,
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo",
                               theory_curve=theory_curve,
                               beam_transfer=crosstrans_beam,
                               mode_transfer_1d = wigglez_transfer_functions_nobeam)

    # try 2D mode compensation
    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_2dbeam2dmodecomp%s" % option_flag,
                               "GBT_15hr_map_combined_cleaned%s" % option_flag,
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo",
                               theory_curve=theory_curve,
                               beam_transfer=crosstrans_beam,
                               mode_transfer_2d = wigglez_transfer_functions_nobeam)

    # doing the whole compensation in 1D
    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_1dbeam1dmodecomp%s" % option_flag,
                               "GBT_15hr_map_combined_cleaned%s" % option_flag,
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo",
                               theory_curve=theory_curve,
                               mode_transfer_1d = wigglez_transfer_functions_withbeam)


def process_autopower_noconvsig():
    # all of the beam, meansub and conv transfer functions
    (trans_beam, trans_beam_mean, trans_beam_meanconv) = ct.find_beam_trans()

    # the data with only mean subtraction
    transfer_functions_noconv = ct.gather_batch_datasim_run("GBT15hr_sim_noconv",
                                                     alt="noconv_")

    gather_batch_data_run("GBT15hr_noconvsig_nocomp", alt_sig="noconv_")

    gather_batch_data_run("GBT15hr_noconvsig_2dbeamcomp", alt_sig="noconv_",
                          beam_transfer=trans_beam_mean)

    gather_batch_data_run("GBT15hr_noconvsig_1dmodecomp", alt_sig="noconv_",
                          mode_transfer_1d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconvsig_2dbeam1dmodecomp", alt_sig="noconv_",
                          beam_transfer=trans_beam_mean,
                          mode_transfer_1d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconvsig_2dmodecomp", alt_sig="noconv_",
                          mode_transfer_2d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconvsig_2dbeam2dmodecomp", alt_sig="noconv_",
                          beam_transfer=trans_beam_mean,
                          mode_transfer_2d=transfer_functions_noconv)


def process_fluxcal():
    gather_batch_data_run("GBT15hr_fluxcal_nocomp",
                          alt_sig="fluxcal_",
                          alt_noise="fluxcal_")
    gather_batch_data_run("GBT15hr_noconv_fluxcal_nocomp",
                          alt_sig="noconv_fluxcal_",
                          alt_noise="noconv_fluxcal_")


# TODO: ADD alt_noise, alt_sig here
def process_autopower_new():
    # all of the beam, meansub and conv transfer functions
    (trans_beam, trans_beam_mean, trans_beam_meanconv) = ct.find_beam_trans()

    # the data with only mean subtraction
    #transfer_functions_noconv = ct.gather_batch_datasim_run("GBT15hr_sim_noconv",
    #                                                 alt="noconv_")

    #gather_batch_data_run("GBT15hr_noconv_nocomp", alt="noconv_",
    #                      subtract_mean=True)

    #gather_batch_data_run("GBT15hr_noconv_2dbeamcomp", alt="noconv_",
    #                      beam_transfer=trans_beam_mean,
    #                      subtract_mean=True)

    #gather_batch_data_run("GBT15hr_noconv_1dmodecomp", alt="noconv_",
    #                      subtract_mean=True,
    #                      mode_transfer_1d=transfer_functions_noconv)

    #gather_batch_data_run("GBT15hr_noconv_2dbeam1dmodecomp", alt="noconv_",
    #                      subtract_mean=True,
    #                      beam_transfer=trans_beam_mean,
    #                      mode_transfer_1d=transfer_functions_noconv)

    #gather_batch_data_run("GBT15hr_noconv_2dmodecomp", alt="noconv_",
    #                      subtract_mean=True,
    #                      mode_transfer_2d=transfer_functions_noconv)

    #gather_batch_data_run("GBT15hr_noconv_2dbeam2dmodecomp", alt="noconv_",
    #                      subtract_mean=True,
    #                      beam_transfer=trans_beam_mean,
    #                      mode_transfer_2d=transfer_functions_noconv)

    transfer_functions_noconv_fluxcal = ct.gather_batch_datasim_run("GBT15hr_sim_noconv_fluxcal",
                                                     alt="noconv_fluxcal_")

    gather_batch_data_run("GBT15hr_noconv_fluxcal_nocomp", alt="noconv_fluxcal_",
                          subtract_mean=False)

    gather_batch_data_run("GBT15hr_noconv_fluxcal_2dbeamcomp", alt="noconv_fluxcal_",
                          beam_transfer=trans_beam_mean,
                          subtract_mean=False)

    gather_batch_data_run("GBT15hr_noconv_fluxcal_1dmodecomp", alt="noconv_fluxcal_",
                          subtract_mean=False,
                          mode_transfer_1d=transfer_functions_noconv_fluxcal)

    gather_batch_data_run("GBT15hr_noconv_fluxcal_2dbeam1dmodecomp", alt="noconv_fluxcal_",
                          subtract_mean=False,
                          beam_transfer=trans_beam_mean,
                          mode_transfer_1d=transfer_functions_noconv_fluxcal)

    gather_batch_data_run("GBT15hr_noconv_fluxcal_2dmodecomp", alt="noconv_fluxcal_",
                          subtract_mean=False,
                          mode_transfer_2d=transfer_functions_noconv_fluxcal)

    gather_batch_data_run("GBT15hr_noconv_fluxcal_2dbeam2dmodecomp", alt="noconv_fluxcal_",
                          subtract_mean=False,
                          beam_transfer=trans_beam_mean,
                          mode_transfer_2d=transfer_functions_noconv_fluxcal)


# TODO: ADD alt_noise, alt_sig here
def process_autopower():
    # all of the beam, meansub and conv transfer functions
    (trans_beam, trans_beam_mean, trans_beam_meanconv) = ct.find_beam_trans()

    # the data with only mean subtraction
    transfer_functions_noconv = ct.gather_batch_datasim_run("GBT15hr_sim_noconv",
                                                     alt="noconv_")

    gather_batch_data_run("GBT15hr_noconv_nocomp", alt="noconv_",
                          subtract_mean=True)

    gather_batch_data_run("GBT15hr_noconv_2dbeamcomp", alt="noconv_",
                          beam_transfer=trans_beam_mean,
                          subtract_mean=True)

    gather_batch_data_run("GBT15hr_noconv_1dmodecomp", alt="noconv_",
                          subtract_mean=True,
                          mode_transfer_1d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconv_2dbeam1dmodecomp", alt="noconv_",
                          subtract_mean=True,
                          beam_transfer=trans_beam_mean,
                          mode_transfer_1d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconv_2dmodecomp", alt="noconv_",
                          subtract_mean=True,
                          mode_transfer_2d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconv_2dbeam2dmodecomp", alt="noconv_",
                          subtract_mean=True,
                          beam_transfer=trans_beam_mean,
                          mode_transfer_2d=transfer_functions_noconv)

    #gather_batch_data_run("GBT15hr_nomeanconv", alt="nomeanconv_",
    #                      subtract_mean=True)

    # data with mean subtraction and degrading
    transfer_functions = ct.gather_batch_datasim_run("GBT15hr_sim")

    gather_batch_data_run("GBT15hr_nocomp")

    gather_batch_data_run("GBT15hr_2dbeamcomp", beam_transfer=trans_beam_meanconv)

    gather_batch_data_run("GBT15hr_1dmodecomp", mode_transfer_1d=transfer_functions)

    gather_batch_data_run("GBT15hr_2dbeam1dmodecomp", beam_transfer=trans_beam_meanconv,
                          mode_transfer_1d=transfer_functions)

    gather_batch_data_run("GBT15hr_2dmodecomp", mode_transfer_2d=transfer_functions)

    gather_batch_data_run("GBT15hr_2dbeam2dmodecomp", beam_transfer=trans_beam_meanconv,
                          mode_transfer_2d=transfer_functions)

    #gather_batch_data_run("GBT15hr", transfer=trans_beam)


if __name__ == '__main__':
    process_autopower_new()
    #process_wigglez(noconv=True, fluxcal=True)
    #process_wigglez(noconv=False, fluxcal=True)
    #process_wigglez(noconv=True)
    #process_wigglez(noconv=False)
    #process_wigglez_fluxcal()
    #process_fluxcal()
    #process_autopower_noconvsig()
    #process_autopower()

#--------------BARF--------------------------------------------------------------------------------------
def find_one_sided_trans():
    gather_batch_one_sided_trans_run("GBTxGBT15hr",
                                     "sim_15hr_combined_cleaned", "sim_15hr",
                                     "GBT_15hr_map_combined_cleaned")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_noconv",
                                     "sim_15hr_combined_cleaned_noconv", "sim_15hr",
                                     "GBT_15hr_map_combined_cleaned_noconv")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_beam",
                                     "sim_15hr_combined_cleaned", "sim_15hr_beam",
                                     "GBT_15hr_map_combined_cleaned")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_noconv_beam",
                                     "sim_15hr_combined_cleaned_noconv", "sim_15hr_beam",
                                     "GBT_15hr_map_combined_cleaned_noconv")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_beam_meansubconv",
                                     "sim_15hr_combined_cleaned", "sim_15hr_beam_meansubconv",
                                     "GBT_15hr_map_combined_cleaned")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_noconv_beam_meansub",
                                     "sim_15hr_combined_cleaned_noconv", "sim_15hr_beam_meansub",
                                     "GBT_15hr_map_combined_cleaned_noconv")

def find_beam_trans():
    from correlate import compile_simpwrspec as csp
    # derive the beam, preprocessing transfer function:
    sim_15hr = csp.gather_batch_sim_run("sim_15hr", "sim_15hr")

    sim_15hr_beam = csp.gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam")

    sim_15hr_beam_mean = csp.gather_batch_sim_run("sim_15hr_beam",
                                              "sim_15hr_beam_mean",
                                              degrade_resolution=False,
                                              subtract_mean=True)

    sim_15hr_beam_meanconv = csp.gather_batch_sim_run("sim_15hr_beam",
                                                  "sim_15hr_beam_meanconv",
                                                  degrade_resolution=True,
                                                  subtract_mean=True)

    trans_beam = tf.calculate_2d_transfer_function(sim_15hr_beam[2],
                                                sim_15hr[2],
                                                "trans_beam")

    trans_beam_mean = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_mean[2],
                                                sim_15hr[2],
                                                "trans_beam_mean")

    trans_beam_meanconv = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_meanconv[2],
                                                sim_15hr[2],
                                                "trans_beam_meanconv")

    # check that it makes sense by reconstructing the simulated power without
    # these effects:
    sim_15hr_corr = csp.gather_batch_sim_run("sim_15hr_beam",
                                         "sim_15hr_beam_meanconv_corrected",
                                         degrade_resolution=True,
                                         subtract_mean=True,
                                         transfer=trans_beam_meanconv)

    return (trans_beam, trans_beam_mean, trans_beam_meanconv)


#def find_modeloss_transfer(alt=""):
#    # now gather the mode loss transfer functions
#    return gather_batch_datasim_run("GBT15hr_sim", alt=alt)


if __name__ == '__main__':
    find_one_sided_trans()
    #find_crossbeam_trans()
    #check_pipeline()
    #find_beam_trans()
    #find_modeloss_transfer()

#------BARF---------------------------------------------
    sys.exit()

    # TODO: make a loop which runs the various relevant cases of this


    sys.exit()
    # WiggleZ xpower with old calibration
    print "PHASEI"
    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned_noconv_fluxcal",
                                "sim_15hr_delta",
                                "sim_15hr",
                                "GBT_15hr_map_combined_cleaned_noconv_fluxcal",
                                "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned_fluxcal",
                                "sim_15hr_delta",
                                "sim_15hr",
                                "GBT_15hr_map_combined_cleaned_fluxcal",
                                "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned_noconv_fluxcal",
                                "sim_15hr_delta",
                                "sim_15hr_beam",
                                "GBT_15hr_map_combined_cleaned_noconv_fluxcal",
                                "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned_fluxcal",
                                "sim_15hr_delta",
                                "sim_15hr_beam",
                                "GBT_15hr_map_combined_cleaned_fluxcal",
                                "WiggleZ_15hr_montecarlo")

    print "PHASEII"
    bq.run_GBTxGBT()
    # crosspower with old calibration
    print "PHASEIII"
    bq.batch_GBTxwigglez_data_run("GBT_15hr_map_combined_cleaned_noconv_fluxcal",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_data_run("GBT_15hr_map_combined_cleaned_fluxcal",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")


    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned_noconv",
                                "sim_15hr_delta",
                                "sim_15hr",
                                "GBT_15hr_map_combined_cleaned_noconv",
                                "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned",
                                "sim_15hr_delta",
                                "sim_15hr",
                                "GBT_15hr_map_combined_cleaned",
                                "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned_noconv",
                                "sim_15hr_delta",
                                "sim_15hr_beam",
                                "GBT_15hr_map_combined_cleaned_noconv",
                                "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_trans_run("sim_15hr_combined_cleaned",
                                "sim_15hr_delta",
                                "sim_15hr_beam",
                                "GBT_15hr_map_combined_cleaned",
                                "WiggleZ_15hr_montecarlo")

    sys.exit()


    bq.batch_GBTxwigglez_data_run("GBT_15hr_map_combined_cleaned_noconv",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")

    bq.batch_GBTxwigglez_data_run("GBT_15hr_map_combined_cleaned",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")


    bq.batch_wigglez_automock_run("WiggleZ_15hr_mock", "WiggleZ_15hr_montecarlo")

    # stuff that needs more work
    #batch_wigglez_automock_run("WiggleZ_15hr_complete_mock", "WiggleZ_15hr_complete_separable_selection")
    #batch_GBTxwigglez_data_run()

def run_GBTxGBT():
    # autopower
    bq.batch_data_run(alt_sig="noconv_fluxcal_", alt_noise="noconv_fluxcal_")
    bq.batch_data_run(alt_sig="fluxcal_", alt_noise="fluxcal_")
    # autopower signal loss:
    bq.batch_data_run(alt_sig="noconv_fluxcal_", alt_noise="noconv_fluxcal_", sim=True)
    bq.batch_data_run(alt_sig="fluxcal_", alt_noise="fluxcal_", sim=True)

    #bq.batch_data_run(alt_sig="noconv_") # use the convolved noise

    # TODO: update these calls
    #bq.batch_data_run(alt="nomeanconv_", subtract_mean=True)
    #bq.batch_data_run(alt="nomeanconv_")
    #bq.batch_data_run(alt="noconv_", subtract_mean=True)
    #bq.batch_data_run(alt="noconv_", sim=True)
    #bq.batch_data_run(alt="noconv_", sim=True, subtract_mean=True)
    #bq.batch_data_run(sim=True)
    #bq.batch_data_run()

