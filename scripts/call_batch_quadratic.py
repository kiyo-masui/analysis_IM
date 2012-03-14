r"""
"GBT_15hr_map_fluxpolcal_cleaned",
"GBT_15hr_map_oldcal_cleaned",
"GBT_15hr_map_fluxcal_cleaned",
"GBT_15hr_map_fdgcal_cleaned"
"GBT_15hr_optimalmap_mapv2fdgcal_cleaned"
"""
from correlate import compile_simpwrspec as csp
from correlate import compile_pwrspec as cp
from correlate import compile_wigglez_xspec as cwx
from correlate import compile_crosspwr_transfer as cxt
from utils import data_paths
import sys

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


def wigglez_xspec(inifile):
    r"""required parameters:
    general:
    `inifile`: .ini file for the power spectrum estimator (bins, etc.) fixed
    `outdir`: place to dump plot data

    for the cross power:
    `selection_function`
    `crosspwr_tag`: unique tag (on top of the input data name) for this run

    for the transfer function:
    `cleaned_simkey`
    `truesignal_simkey`
    `truesignal_weightkey`
    `reference_simkey`
    `reference_weightkey`
    `crosspwr_trans_tag`: unique tag (on top of the input data name) for this run
    `generate`: True is regenerate the cache for this
    applied as:

    cleaned_simkey(map) * cleaned_simkey(weight) x
    truesignal_weightkey * truesignal_simkey
    divided by:
    reference_simkey * reference_weightkey x truesignal_weightkey * truesignal_simkey

    """
    datapath_db = data_paths.DataPath()

    #inifile=None, generate=False, alttag=None,
    #              mode_transfer_1d=None, mode_transfer_2d=None,
    #              beam_transfer=None):
    #datapath_db = data_paths.DataPath()
    # cleaned_mapkey "GBT_15hr_map_fdgcal_cleaned_noconv_combined"
    # tack "_xWigglez" onto the output_tag name and an optional alttag after that

    #cleaned_mapkey, wigglez_signalkey, wigglez_mockkey, wigglez_selection

    mapname = "GBT_15hr_map_fdgcal_cleaned_noconv_combined"
    wigglez_map_key = "WiggleZ_15hr_binned_delta"
    wigglez_mock_key = "WiggleZ_15hr_delta_mock"
    wigglez_selection_key = "WiggleZ_15hr_montecarlo"

    cwx.batch_gbtxwigglez_data_run(mapname, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=inifile, datapath_db=datapath_db,
                               outdir="./plots/",
                               output_tag=output_tag,
                               beam_transfer=beam_transfer,
                               mode_transfer_1d=mode_transfer_1d,
                               mode_transfer_2d=mode_transfer_2d,
                               theory_curve=None)


def wigglez_crosspwr_permutations(inifile=None, generate=False):
    basemaps = ["GBT_15hr_map_fdgcal_plussim_cleaned",
                "GBT_15hr_map_fdgcal_cleanedplussim_cleaned",
                "GBT_15hr_map_fdgcal_cleaned"]
    treatments = ["", "_noconv"]

    cwx.call_batch_gbtxwigglez_data_run(basemaps, treatments,
                                        "sim_15hr_oldmap_str_delta:1",
                                        "WiggleZ_15hr_delta_mock",
                                        "WiggleZ_15hr_montecarlo",
                                        inifile=inifile, generate=generate,
                                        outdir="./plots/",
                                        alttag="clfgcorrtest")


def wigglez_transfer_permutations(inifile=None, generate=False):
    rootsim = "sim_15hr_oldmap_str"
    selection_function = "WiggleZ_15hr_montecarlo"

    cleaned_simkey = "GBT_15hr_map_fdgcal_cleaned_sims_noconv_combined"
    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=True, inifile=inifile,
                                            generate=generate, simindex="0",
                                            alttag="clsim-xtrans-withbeam")

    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=False, inifile=inifile,
                                            generate=generate, simindex="0",
                                            alttag="clsim-xtrans-nobeam")

    cleaned_simkey = "GBT_15hr_map_fdgcal_cleaned_sims_combined"
    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=True, inifile=inifile,
                                            generate=generate, simindex="0",
                                            alttag="clsim-xtrans-withbeam")

    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=False, inifile=inifile,
                                            generate=generate, simindex="0",
                                            alttag="clsim-xtrans-nobeam")

    cleaned_simkey = "GBT_15hr_map_fdgcal_plussim_cleaned_noconv_combined"
    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=True, inifile=inifile,
                                            generate=generate,
                                            alttag="corrfg-xtrans-withbeam")

    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=False, inifile=inifile,
                                            generate=generate,
                                            alttag="corrfg-xtrans-nobeam")

    cleaned_simkey = "GBT_15hr_map_fdgcal_plussim_cleaned_combined"
    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=True, inifile=inifile,
                                            generate=generate,
                                            alttag="corrfg-xtrans-withbeam")

    transfer = cxt.wigglez_crosspwr_transfer_run(cleaned_simkey, rootsim, selection_function,
                                            include_beam=False, inifile=inifile,
                                            generate=generate,
                                            alttag="corrfg-xtrans-nobeam")


if __name__ == '__main__':
    inifile = "input/ers/batch_quadratic/default.ini"
    inifile_phys = "input/ers/batch_quadratic/default_physical.ini"

    # run simulations
    #sim_autopwr(inifile=inifile, inifile_phys=inifile_phys, generate=False)
    #sim_crosspwr(inifile=inifile, generate=False)

    # calculate the transfer function using cleaned(realmap + sim) x sim
    # for the cross-power
    wigglez_transfer_permutations(inifile=inifile, generate=False)
    #wigglez_crosspwr_permutations(inifile=inifile, generate=False)

