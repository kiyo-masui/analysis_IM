# TODO: WARNING, this needs to be updated to the new outdir structure
# but ... we may not need this function anymore?
def call_batch_gbtxwigglez_data_run(basemaps, treatments, wigglez_map_key,
                                    wigglez_mock_key, wigglez_selection_key,
                                    inifile=None, generate=False,
                                    outdir="./plots/", alttag=None):
    r"""Call a chunk of batch data runs for e.g. different map products

    No compensation can be applied using this function because
    batch_gbtxwigglez_data_run should be called on a case by case basis for
    transfer functions associated with different map runs

    This is useful for building up the cache of various cross-powers
    """
    datapath_db = data_paths.DataPath()

    for base in basemaps:
        for treatment in treatments:
            mapname = base + treatment + "_combined"

            output_tag = mapname
            if alttag:
                output_tag += "_" + output_tag

            if generate:
                output_tag = None

            print mapname
            batch_gbtxwigglez_data_run(mapname, wigglez_map_key,
                            wigglez_mock_key, wigglez_selection_key,
                            inifile=inifile,
                            datapath_db=datapath_db,
                            output_tag=output_tag,
                            outdir=outdir,
                            beam_transfer=None,
                            mode_transfer_1d=None,
                            mode_transfer_2d=None,
                            theory_curve=None)

def call_data_autopower(basemaps, treatments, inifile=None, generate=False,
                        outdir="./plots/", mode_transfer_1d=None,
                        mode_transfer_2d=None, beam_transfer=None, alttag=None):
    r"""Call a chunk of batch data runs for e.g. different map products
    """
    datapath_db = data_paths.DataPath()

    for base in basemaps:
        for treatment in treatments:
            output_tag = base + treatment
            if alttag:
                output_tag += "_" + output_tag

            output_root = "%s/%s/" % (outdir, output_tag)

            if generate:
                output_tag = None

            batch_gbtpwrspec_data_run(mapname,
                                 inifile=inifile,
                                 datapath_db=datapath_db,
                                 output_tag=output_tag,
                                 outdir=output_root,
                                 beam_transfer=beam_transfer,
                                 mode_transfer_1d=mode_transfer_1d,
                                 mode_transfer_2d=mode_transfer_2d)



