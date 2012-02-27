from correlate import batch_quadratic as bq
import sys


def sim_autopower(inifile=None, inifile_phys=None):
    r"""run all of the auto-power theory cases
    """
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str', 'sim_22hr_oldmap_str',
                'sim_1hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']

    for base in basesims:
        bq.batch_physical_sim_run("%s_physical" % base,
                                  inifile=inifile_phys)

        for treatment in treatments:
            mapname = base + treatment
            weight = "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes"
            bq.batch_sim_run(mapname, mapname,
                             weight, weight, inifile=inifile)


def sim_crosspower(inifile=None):
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str', 'sim_22hr_oldmap_str',
                'sim_1hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']

    for base in basesims:
        for treatment in treatments:
            map1 = base + treatment
            map2 = base + "_delta"
            print map1, map2
            bq.batch_sim_run(map1, map2,
                         "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes",
                         "WiggleZ_15hr_montecarlo", inifile=inifile)


def sim_one_sided_trans(basemap, basesim, inifile=None):
    left_treatments = ["", "_noconv"]
    right_treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                        '_beam_meansubconv']

    for lt in left_treatments:
        for rt in right_treatments:
            lmap = "%s_cleaned_sims%s_combined" % (basemap, lt)
            rmap = "%s%s" % (basesim, rt)
            wmap = "%s_cleaned%s_combined" % (basemap, lt)
            bq.batch_one_sided_trans_run(lmap, rmap, wmap, inifile=inifile)


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
    inifile = "input/ers/batch_quadratic/default.ini"
    inifile_phys = "input/ers/batch_quadratic/default_physical.ini"

    sim_autopower(inifile=inifile, inifile_phys=inifile_phys)

    # TODO: make a loop which runs the various relevant cases of this
    sim_one_sided_trans("GBT_15hr_map_fluxpolcal",
                        "sim_15hr_oldmap_str", inifile=inifile)

    sim_crosspower(inifile=inifile)

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

