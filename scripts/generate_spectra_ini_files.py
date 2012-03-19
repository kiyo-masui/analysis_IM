def call_sim_autopower(basesims, treatments, weight, inifile=None,
                       inifile_phys=None, generate=False,
                       outdir="./plots/"):
    r"""run all of the auto-power theory cases
    """
    datapath_db = data_paths.DataPath()

    for base in basesims:
        mapname = "%s_physical" % base

        usecache_output_tag = None
        if not generate:
            usecache_output_tag = mapname

        batch_physical_sim_run(mapname,
                               inifile=inifile_phys,
                               datapath_db=datapath_db,
                               outdir=outdir,
                               usecache_output_tag=usecache_output_tag)

        for treatment in treatments:
            mapname = base + treatment

            usecache_output_tag = None
            if not generate:
                usecache_output_tag = mapname

            batch_sim_run(mapname, mapname,
                          weight, weight, inifile=inifile,
                          datapath_db=datapath_db,
                          outdir=outdir,
                          usecache_output_tag=usecache_output_tag)


def call_sim_crosspower(basesims, treatments, weight,
                        selection_function,
                        inifile=None,
                        inifile_phys=None, generate=False,
                        outdir="./plots/"):
    r"""run all of the cross-power theory cases
    """
    datapath_db = data_paths.DataPath()

    for base in basesims:
        for treatment in treatments:
            mapname = base + treatment
            deltaname = base + "_delta"

            usecache_output_tag = None
            if not generate:
                usecache_output_tag = mapname + "_xWigglez"

            batch_sim_run(mapname, deltaname,
                          weight, selection_function, inifile=inifile,
                          datapath_db=datapath_db,
                          outdir=outdir,
                          usecache_output_tag=usecache_output_tag)

def sim_autopwr(inifile_phys=None, inifile=None, generate=False):
    # add 'sim_22hr_oldmap_str', 'sim_1hr_oldmap_str'
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']
    weight = "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes"

    call_sim_autopower(basesims, treatments, weight, inifile=inifile,
                           inifile_phys=inifile_phys, generate=generate)


def sim_crosspwr(inifile=None, generate=False):
    # add 'sim_22hr_oldmap_str', 'sim_1hr_oldmap_str'
    basesims = ['sim_15hr_oldmap_ideal', 'sim_15hr_oldmap_nostr',
                'sim_15hr_oldmap_str']
    treatments = ['_temperature', '_beam', '_beam_conv', '_beam_meansub',
                  '_beam_meansubconv']
    weight = "GBT_15hr_map_fluxpolcal_cleaned_combined:weight;0modes"
    selection_function = "WiggleZ_15hr_montecarlo"

    call_sim_crosspower(basesims, treatments, weight, selection_function,
                            inifile=inifile, generate=generate)


#if __name__ == '__main__':
#    inifile = "input/ers/batch_quadratic/default.ini"
#    inifile_phys = "input/ers/batch_quadratic/default_physical.ini"
#    sim_autopwr(inifile=inifile, inifile_phys=inifile_phys, generate=False)
#    sim_crosspwr(inifile=inifile, generate=False)



