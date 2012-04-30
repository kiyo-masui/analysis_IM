from correlate import compile_simpwrspec as csp
from utils import data_paths
import sys

def sim_crosspower(inifile=None):
    #basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
    #            'sim_15hr_oldmap_str', 'sim_22hr_oldmap_str',
    #            'sim_1hr_oldmap_str']
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str']
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


if __name__ == '__main__':
    inifile = "input/ers/batch_quadratic/default.ini"
    inifile_phys = "input/ers/batch_quadratic/default_physical.ini"

    # add 'sim_22hr_oldmap_str', 'sim_1hr_oldmap_str'
    basesims = ['sim_15hr_oldmap_ideal']
    #basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
    #            'sim_15hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']
    weight = "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes"

    csp.sim_autopower(basesims, treatments, weight, inifile=inifile,
                      inifile_phys=inifile_phys, generate=True)

    #sim_autopower(inifile=inifile, inifile_phys=inifile_phys)
    #sim_crosspower(inifile=inifile)
    #sim_one_sided_trans("GBT_15hr_map_fluxpolcal",
    #                    "sim_15hr_oldmap_str", inifile=inifile)

