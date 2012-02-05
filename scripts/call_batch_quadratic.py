from correlate import batch_quadratic as bq

def sim_autopower():
    r"""run all of the auto-power theory cases
    sim = dd, vv
    simvel = dd, vv, streaming
    ideal = dd only
    1800 power spectra
    """

    # run autopower on the physical volumes
    bq.batch_physical_sim_run("sim_15hr_physical")
    bq.batch_physical_sim_run("simvel_15hr_physical")
    bq.batch_physical_sim_run("simideal_15hr_physical")

    # run autopower on the observed volumes
    bq.batch_sim_run("sim_15hr")
    bq.batch_sim_run("simvel_15hr")
    bq.batch_sim_run("simideal_15hr")

    # run autopower on the observed volumes with beam applied
    bq.batch_sim_run("sim_15hr_beam")
    bq.batch_sim_run("simvel_15hr_beam")
    bq.batch_sim_run("simideal_15hr_beam")

    # run autopower on the observed volumes with beam and meansub applied
    bq.batch_sim_run("sim_15hr_beam", degrade_resolution=False, subtract_mean=True)
    bq.batch_sim_run("simvel_15hr_beam", degrade_resolution=False, subtract_mean=True)
    bq.batch_sim_run("simideal_15hr_beam", degrade_resolution=False, subtract_mean=True)

    # run autopower on the observed volumes with beam and conv applied
    bq.batch_sim_run("sim_15hr_beam", degrade_resolution=True, subtract_mean=False)
    bq.batch_sim_run("simvel_15hr_beam", degrade_resolution=True, subtract_mean=False)
    bq.batch_sim_run("simideal_15hr_beam", degrade_resolution=True, subtract_mean=False)

    # run autopower on the observed volumes with beam, meansub and conv applied
    bq.batch_sim_run("sim_15hr_beam", degrade_resolution=True, subtract_mean=True)
    bq.batch_sim_run("simvel_15hr_beam", degrade_resolution=True, subtract_mean=True)
    bq.batch_sim_run("simideal_15hr_beam", degrade_resolution=True, subtract_mean=True)


def sim_crosspower():
    r"""Do not calculate the simideal case"""
    bq.batch_genericsim_run("sim_15hr", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("sim_15hr_beam", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("sim_15hr_beam_meansub", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("sim_15hr_beam_meansubconv", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("sim_15hr_beam_conv", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("simvel_15hr", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("simvel_15hr_beam", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("simvel_15hr_beam_meansub", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("simvel_15hr_beam_meansubconv", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")

    bq.batch_genericsim_run("simvel_15hr_beam_conv", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo")


def sim_one_sided_trans():
    bq.batch_one_sided_trans_run("sim_15hr_combined_cleaned", "sim_15hr",
                              "GBT_15hr_map_combined_cleaned")
    bq.batch_one_sided_trans_run("sim_15hr_combined_cleaned_noconv", "sim_15hr",
                              "GBT_15hr_map_combined_cleaned_noconv")
    bq.batch_one_sided_trans_run("sim_15hr_combined_cleaned", "sim_15hr_beam",
                              "GBT_15hr_map_combined_cleaned")
    bq.batch_one_sided_trans_run("sim_15hr_combined_cleaned_noconv", "sim_15hr_beam",
                              "GBT_15hr_map_combined_cleaned_noconv")
    bq.batch_one_sided_trans_run("sim_15hr_combined_cleaned", "sim_15hr_beam_meansubconv",
                              "GBT_15hr_map_combined_cleaned")
    bq.batch_one_sided_trans_run("sim_15hr_combined_cleaned_noconv", "sim_15hr_beam_meansub",
                              "GBT_15hr_map_combined_cleaned_noconv")


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


if __name__ == '__main__':
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

    sys.exit()
    sim_crosspower()

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

