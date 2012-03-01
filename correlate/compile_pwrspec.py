r"""Code to run large batches of quadratic estimators on combinations of
data/sims
"""
from utils import data_paths
from correlate import pwrspec_estimation as pe
from correlate import transfer_function as tf
from utils import batch_handler
import copy
from correlate import batch_quadratic as bq


def batch_data_run(map_key, inifile=None, datapath_db=None,
                   compile_tag=None, beam_transfer=None,
                   outdir="./plot_data_v2",
                   mode_transfer_1d=None,
                   mode_transfer_2d=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(map_key, "pair;type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    for treatment in map_cases['treatment']:
        uniq_pairs = data_paths.GBTauto_cross_pairs(map_cases['pair'],
                                                    map_cases['pair'],
                                                    cross_sym="_with_")

        # TODO: make this more elegant
        transfer_2d = None
        if (mode_transfer_2d is not None) and (beam_transfer is None):
            transfer_2d = mode_transfer_2d[mode_num][1]

        if (mode_transfer_2d is None) and (beam_transfer is not None):
            transfer_2d = beam_transfer

        if (mode_transfer_2d is not None) and (beam_transfer is not None):
            transfer_2d = mode_transfer_2d[mode_num][1] * beam_transfer

        pwr_1d = []
        pwr_2d = []
        pwr_1d_from_2d = []
        for item in uniq_pairs:
            input = {}
            input['map1_key'] = "%s:%s;map;%s" % (map_key, item[0], treatment)
            input['map2_key'] = "%s:%s;map;%s" % (map_key, item[1], treatment)
            input['noiseinv1_key'] = "%s:%s;weight;%s" % (map_key, item[0], treatment)
            input['noiseinv2_key'] = "%s:%s;weight;%s" % (map_key, item[1], treatment)
            files = convert_keydict_to_filedict(input, db=datapath_db)

            pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                                  files['map2_key'],
                                                  files['noiseinv1_key'],
                                                  files['noiseinv2_key'],
                                                  inifile=inifile)

            if compile_tag:
                pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                      transfer=transfer_2d))

                pwr_2d.append(pwr2d_run)
                pwr_1d.append(pwr1d_run)

                mtag = compile_tag + "_%s" % treatment
                if mode_transfer_1d is not None:
                    transfunc = mode_transfer_1d[mode_num][0]
                else:
                    transfunc = None

                pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                                     outdir="./plot_data",
                                     apply_1d_transfer=transfunc)

    if compile_tag:
        return None
    else:
        caller.multiprocess_stack()
        return None


def batch_GBTxwigglez_data_run(gbt_map_key, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=None, datapath_db=None,
                               compile_tag=None):
                               beam_transfer=None,
                               outdir="./plot_data", mode_transfer_1d=None,
                               mode_transfer_2d=None,
                               theory_curve=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(gbt_map_key, "type;treatment")
    mock_cases = datapath_db.fileset_cases(wigglez_mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    for treatment in map_cases['treatment']:
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
            files = convert_keydict_to_filedict(input, db=datapath_db)

            caller.execute(files['map1_key'], files['map2_key'],
                           files['noiseinv1_key'], files['noiseinv2_key'],
                           inifile=inifile)

            if compile_tag:
                pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                      transfer=transfer_2d))

                pwr_2d.append(pwr2d_run)
                pwr_1d.append(pwr1d_run)

        if compile_tag:
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
            files = convert_keydict_to_filedict(input, db=datapath_db)

            pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                              files['map2_key'],
                                              files['noiseinv1_key'],
                                              files['noiseinv2_key'],
                                              inifile=inifile)

            pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run, transfer=transfer_2d)

            pwr_1d = pwr1d_run['binavg']
            pwr_1d_from_2d = pwr1d_from_2d['binavg']
            if mode_transfer_1d is not None:
                pwr_1d /= mode_transfer_1d[mode_num][0]
                pwr_1d_from_2d /= mode_transfer_1d[mode_num][0]

            # assume that they all have the same binning
            bin_left = pwr1d_run['bin_left']
            bin_center = pwr1d_run['bin_center']
            bin_right = pwr1d_run['bin_right']
            counts_histo = pwr1d_run['counts_histo']

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
    if not compile_tag:
        caller.multiprocess_stack()
        return None


def batch_wigglez_automock_run(mock_key, sel_key,
                               inifile=None, datapath_db=None,
                               compile_tag=None):
    r"""TODO: make this work; wrote this but never really needed it yet
    """
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    for index in mock_cases['realization']:
        input = {}
        input['map1_key'] = "%s:%s" % (mock_key, index)
        input['map2_key'] = "%s:%s" % (mock_key, index)
        input['noiseinv1_key'] = sel_key
        input['noiseinv2_key'] = sel_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    caller.multiprocess_stack()
