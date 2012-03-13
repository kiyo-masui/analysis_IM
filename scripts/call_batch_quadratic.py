from correlate import compile_simpwrspec as csp
from correlate import compile_pwrspec as cp
from correlate import compile_wigglez_xspec as cwx
from correlate import compile_corrfg_transfer as cct
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


def corrfg_crosspwr_transfer_run(mapname, rootsim, selection_function,
                                 simindex="1", weightmap="15modes",
                                 include_beam=True,
                                 inifile=None, generate=False, alttag=None):

    datapath_db = data_paths.DataPath()

    delta_simkey = "%s_delta:%s" % (rootsim, simindex)
    deltasim_weightkey = selection_function

    if include_beam:
        basesim_key = "%s_temperature:%s" % (rootsim, simindex)
        output_tag = mapname + "_corrfgtranswithbeam"
    else:
        basesim_key = "%s_beam:%s" % (rootsim, simindex)
        output_tag = mapname + "_corrfgtransnobeam"

    if alttag:
        output_tag += "_" + alttag

    basesim_weightkey = "%s:weight;%s" % (mapname, weightmap)

    if generate:
        output_tag = None

    return cct.batch_corrfg_transfer_run(mapname, delta_simkey,
                                      deltasim_weightkey,
                                      basesim_key, basesim_weightkey,
                                      inifile=inifile, datapath_db=datapath_db,
                                      outdir="./plots/",
                                      output_tag=output_tag)


def wigglez_crosspwr(inifile=None, generate=False):
    rootsim = "sim_15hr_oldmap_str"
    selection_function = "WiggleZ_15hr_montecarlo"

    mapname = "GBT_15hr_map_fdgcal_plussim_cleaned_noconv_combined"
    transfer = corrfg_crosspwr_transfer_run(mapname, rootsim, selection_function,
                                            include_beam=True, inifile=inifile,
                                            generate=generate)

    transfer = corrfg_crosspwr_transfer_run(mapname, rootsim, selection_function,
                                            include_beam=False, inifile=inifile,
                                            generate=generate)

    mapname = "GBT_15hr_map_fdgcal_plussim_cleaned_combined"
    transfer = corrfg_crosspwr_transfer_run(mapname, rootsim, selection_function,
                                            include_beam=True, inifile=inifile,
                                            generate=generate)

    transfer = corrfg_crosspwr_transfer_run(mapname, rootsim, selection_function,
                                            include_beam=False, inifile=inifile,
                                            generate=generate)

if __name__ == '__main__':
    inifile = "input/ers/batch_quadratic/default.ini"
    inifile_phys = "input/ers/batch_quadratic/default_physical.ini"

    # run simulations
    #sim_autopwr(inifile=inifile, inifile_phys=inifile_phys, generate=False)
    #sim_crosspwr(inifile=inifile, generate=False)

    # calculate the transfer function using cleaned(realmap + sim) x sim
    # for the cross-power
    wigglez_crosspwr(inifile=inifile)

