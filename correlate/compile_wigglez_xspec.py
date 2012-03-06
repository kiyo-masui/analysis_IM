from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
import copy
from correlate import batch_quadratic as bq
# TODO: figure out theory curve stuff

def batch_GBTxwigglez_data_run(gbt_map_key, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=None, datapath_db=None,
                               outdir="./plot_data_v2",
                               usecache_output_tag=None,
                               beam_transfer=None,
                               mode_transfer_1d=None,
                               mode_transfer_2d=None,
                               theory_curve=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(gbt_map_key, "type;treatment")
    mock_cases = datapath_db.fileset_cases(wigglez_mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if usecache_output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    for treatment in map_cases['treatment']:
        # TODO: make this more elegant
        transfer_2d = None
        if (mode_transfer_2d is not None) and (beam_transfer is None):
            transfer_2d = mode_transfer_2d[mode_num][2]

        if (mode_transfer_2d is None) and (beam_transfer is not None):
            transfer_2d = beam_transfer

        if (mode_transfer_2d is not None) and (beam_transfer is not None):
            transfer_2d = mode_transfer_2d[mode_num][2] * beam_transfer

        pwr_1d = []
        pwr_2d = []
        pwr_1d_from_2d = []
        for index in mock_cases['realization']:
            input = {}
            input['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
            input['map2_key'] = "%s:%s" % (wigglez_mock_key, index)
            input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
            input['noiseinv2_key'] = wigglez_selection_key
            files = bq.convert_keydict_to_filedict(input, db=datapath_db)

            pwrspec_out = caller.execute(files['map1_key'], files['map2_key'],
                                         files['noiseinv1_key'],
                                         files['noiseinv2_key'],
                                         inifile=inifile)

            if usecache_output_tag:
                pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                      transfer=transfer_2d))

                pwr_2d.append(pwrspec_out[0])
                pwr_1d.append(pwrspec_out[1])

        if usecache_output_tag:
            if mode_transfer_1d is not None:
                transfunc = mode_transfer_1d[mode_num][0]
            else:
                transfunc = None

            mtag = tag + "_%dmodes_mock" % mode_num
            mean1dmock, std1dmock, covmock = pe.summarize_pwrspec(pwr_1d,
                              pwr_1d_from_2d, pwr_2d, mtag,
                              outdir="./plot_data",
                              apply_1d_transfer=transfunc)

        # now recover the xspec with the real data
        input = {}
        input['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
        input['map2_key'] = wigglez_map_key
        input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
        input['noiseinv2_key'] = wigglez_selection_key
        files = bq.convert_keydict_to_filedict(input, db=datapath_db)

        pwrspec_out_signal = caller.execute(files['map1_key'],
                                            files['map2_key'],
                                            files['noiseinv1_key'],
                                            files['noiseinv2_key'],
                                            inifile=inifile)

        if usecache_output_tag:
            pwr1d_from_2d = pe.convert_2d_to_1d(pwrspec_out_signal[0],
                                                transfer=transfer_2d)

            pwr_1d = pwrspec_out_signal[1]['binavg']
            pwr_1d_from_2d = pwr1d_from_2d['binavg']
            if mode_transfer_1d is not None:
                pwr_1d /= mode_transfer_1d[mode_num][0]
                pwr_1d_from_2d /= mode_transfer_1d[mode_num][0]

            # assume that they all have the same binning
            bin_left = pwrspec_out_signal[1]['bin_left']
            bin_center = pwrspec_out_signal[1]['bin_center']
            bin_right = pwrspec_out_signal[1]['bin_right']
            counts_histo = pwrspec_out_signal[1]['counts_histo']

            filename = "./plot_data/%s_%dmodes.dat" % (tag, mode_num)
            outfile = open(filename, "w")
            for specdata in zip(bin_left, bin_center,
                                bin_right, counts_histo, pwr_1d, pwr_1d_from_2d,
                                mean1dmock, std1dmock):
                outfile.write(("%10.15g " * 8 + "\n") % specdata)
            outfile.close()

            if theory_curve is not None:
                restrict = np.where(np.logical_and(bin_center > 0.09,
                                                   bin_center < 1.1))
                res_slice = slice(min(restrict[0]), max(restrict[0]))

                #restricted_cov = covmock[np.where(restrict)[0][np.newaxis, :]][0]
                print "AMP:", tag, mode_num, utils.ampfit(pwr_1d_from_2d[res_slice],
                                                    covmock[res_slice, res_slice],
                                                    theory_curve[res_slice])
    if not usecache_output_tag:
        caller.multiprocess_stack()

    return None


def call_batch_GBTxwigglez_data_run(basemaps, treatments, wigglez_map_key,
                                    wigglez_mock_key, wigglez_selection_key,
                                    inifile=None, generate=False,
                                    outdir="./plot_data_v2/",
                                    mode_transfer_1d=None,
                                    mode_transfer_2d=None,
                                    beam_transfer=None):

    datapath_db = data_paths.DataPath()
    # TODO: put transfer functions in the naming tags

    for base in basemaps:
        for treatment in treatments:
            mapname = base + treatment

            usecache_output_tag = None
            if not generate:
                usecache_output_tag = mapname

            batch_GBTxwigglez_data_run(mapname, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=inifile, datapath_db=datapath_db,
                               outdir="./plot_data_v2",
                               usecache_output_tag=None,
                               beam_transfer=beam_transfer,
                               mode_transfer_1d=mode_transfer_1d,
                               mode_transfer_2d=mode_transfer_2d,
                               theory_curve=None)
