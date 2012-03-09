from correlate import compile_simpwrspec as csp
from correlate import compile_pwrspec as cp
from correlate import compile_wigglez_xspec as cwx
from utils import data_paths
import sys


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


def sim_autopwr(inifile_phys=None, inifile=None, generate=False):
    # add 'sim_22hr_oldmap_str', 'sim_1hr_oldmap_str'
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']
    weight = "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes"

    csp.call_sim_autopower(basesims, treatments, weight, inifile=inifile,
                           inifile_phys=inifile_phys, generate=generate)


def sim_crosspwr(inifile=None, generate=False):
    # add 'sim_22hr_oldmap_str', 'sim_1hr_oldmap_str'
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']
    weight = "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes"
    selection_function = "WiggleZ_15hr_montecarlo"

    csp.call_sim_crosspower(basesims, treatments, weight, selection_function,
                            inifile=inifile, generate=generate)


def wigglez_xspec(inifile=None, generate=False):
    basemaps = ["GBT_15hr_map_fluxpolcal_cleaned",
                "GBT_15hr_map_oldcal_cleaned",
                "GBT_15hr_map_fluxcal_cleaned",
                "GBT_15hr_map_fdgcal_cleaned",
                "GBT_15hr_optimalmap_mapv2fdgcal_cleaned"]
    #basemaps = ["GBT_15hr_map_fdgcal_cleaned"]
    treatments = ["", "_noconv"]

    cwx.call_batch_gbtxwigglez_data_run(basemaps, treatments,
                                        "WiggleZ_15hr_delta_binned_data",
                                        "WiggleZ_15hr_delta_mock",
                                        "WiggleZ_15hr_montecarlo",
                                        inifile=inifile, generate=generate,
                                        outdir="./plots/", mode_transfer_1d=None,
                                        mode_transfer_2d=None, beam_transfer=None)


def gbt_autopwr(inifile=None, generate=False):
    # left out (gridding): "GBT_15hr_optimalmap_fluxpolcal_cleaned"
    # left out (D->C map): "GBT_15hr_optimalmap_mapv2oldcal_cleaned"
    # not done yet: "GBT_15hr_optimalmap_mapv2fdgcalmoderm_cleaned"
    #basemaps = ["GBT_22hr_map_fluxpolcal_cleaned"]
    #basemaps = ["GBT_1hr_map_fluxpolcal_cleaned"]
    basemaps = ["GBT_15hr_map_fluxpolcal_cleaned",
                "GBT_15hr_map_oldcal_cleaned",
                "GBT_15hr_map_fluxcal_cleaned",
                "GBT_15hr_map_fdgcal_cleaned",
                "GBT_15hr_optimalmap_mapv2fdgcal_cleaned"]
    #basemaps = ["GBT_15hr_map_fdgcal_cleaned"]
    treatments = ["", "_sims", "_noconv", "_sims_noconv"]

    cp.call_data_autopower(basemaps, treatments, inifile=inifile,
                           generate=generate,
                           outdir="./plots/", mode_transfer_1d=None,
                           mode_transfer_2d=None, beam_transfer=None)


if __name__ == '__main__':
    inifile = "input/ers/batch_quadratic/default.ini"
    inifile_phys = "input/ers/batch_quadratic/default_physical.ini"

    #sim_autopwr(inifile=inifile, inifile_phys=inifile_phys, generate=False)
    #sim_crosspwr(inifile=inifile, generate=False)
    wigglez_xspec(inifile=inifile, generate=False)
    gbt_autopwr(inifile=inifile, generate=False)

    #sim_crosspower(inifile=inifile)
    #sim_one_sided_trans("GBT_15hr_map_fluxpolcal",
    #                    "sim_15hr_oldmap_str", inifile=inifile)


