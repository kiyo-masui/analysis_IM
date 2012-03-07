from utils import data_paths
from correlate import pwrspec_estimation as pe
from correlate import transfer_function as tf
from utils import batch_handler
import copy
from correlate import batch_quadratic as bq


def batch_one_sided_trans_run(modeloss_simkey, sim_key,
                              modeloss_weight_root,
                              inifile=None, datapath_db=None,
                              compile_tag=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(sim_key, "type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    dbkeydict = {}
    dbkeydict['map1_key'] = "%s:0" % sim_key
    dbkeydict['map2_key'] = "%s:0" % sim_key
    dbkeydict['noiseinv1_key'] = "%s:weight;0modes" % (modeloss_weight_root)
    # TODO: should these be ones or weight
    dbkeydict['noiseinv2_key'] = "%s:ones;0modes" % (modeloss_weight_root)
    files = data_paths.convert_dbkeydict_to_filedict(dbkeydict, datapath_db=datapath_db)

    pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                          files['map2_key'],
                                          files['noiseinv1_key'],
                                          files['noiseinv2_key'],
                                          inifile=inifile)

    pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run)

    pwr_1d_zero = pwr1d_run['binavg']
    pwr_2d_zero = pwr2d_run['binavg']
    pwr_1d_from_2d_zero = pwr1d_from_2d['binavg']

    transfer_functions = {}
    for treatment in map_cases['treatment']:
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:map;%s" % (modeloss_simkey, treatment)
        dbkeydict['map2_key'] = "%s:0" % sim_key

        dbkeydict['noiseinv1_key'] = "%s:weight;%s" % (modeloss_weight_root, \
                                                   treatment)

        dbkeydict['noiseinv2_key'] = "%s:ones;%s" % (modeloss_weight_root, \
                                                 treatment)
        files = data_paths.convert_dbkeydict_to_filedict(dbkeydict, datapath_db=datapath_db)

        pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                              files['map2_key'],
                                              files['noiseinv1_key'],
                                              files['noiseinv2_key'],
                                              inifile=inifile)
        if compile_tag:
            pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run)

            pwr_1d = pwr1d_run['binavg']
            pwr_1d_from_2d = pwr1d_from_2d['binavg']
            trans1d_mode = pwr_1d/pwr_1d_zero
            trans1d_from2d_mode = pwr_1d_from_2d/pwr_1d_from_2d_zero
            trans2d_mode = pwr2d_run['binavg']/pwr_2d_zero

            transfer_functions[mode_num] = (trans1d_mode,
                                            trans1d_from2d_mode,
                                            trans2d_mode)

            # assume that they all have the same binning
            bin_left = pwr1d_run['bin_left']
            bin_center = pwr1d_run['bin_center']
            bin_right = pwr1d_run['bin_right']
            counts_histo = pwr1d_run['counts_histo']

            filename = "./plot_data/%s_%dmodes_trans.dat" % (tag, mode_num)
            outfile = open(filename, "w")
            for specdata in zip(bin_left, bin_center,
                                bin_right, counts_histo, pwr_1d, trans1d_mode,
                                trans1d_from2d_mode):
                outfile.write(("%10.15g " * 7 + "\n") % specdata)
            outfile.close()

    if not compile_tag:
        caller.multiprocess_stack()
        retval = None
    else:
        retval = transfer_functions

    return retval


def gather_batch_datasim_run(tag, subtract_mean=False,
                             degrade_resolution=False, unitless=True,
                             return_3d=False, truncate=False, window=None,
                             n_modes=None, refinement=2, pad=5, order=2,
                             outdir="./plot_data", alt=""):
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_simulations")
    print "reading from: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    transfer_functions = {}
    for mode_num in range(0,55,5):
        mapsim = "sim_15hr"
        map1_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)
        map2_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)
        noise1_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)
        noise2_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)

        (pairlist, pairdict) = \
                data_paths.cross_maps(map1_key, map2_key,
                                      noise1_key, noise2_key,
                                      map_suffix=";map",
                                      noise_inv_suffix=";noise_inv",
                                      cross_sym="_x_",
                                      pair_former="GBTauto_cross_pairs",
                                      ignore=['param'],
                                      tag1prefix=map1_key + "_",
                                      tag2prefix=map2_key + "_",
                                      verbose=False)

        pwr_1d = []
        pwr_2d = []
        pwr_1d_from_2d = []
        for item in pairdict.keys():
            pairrun = pairdict[item]
            print pairrun['tag1'], pairrun['map1'], pairrun['noise_inv1']

            pwr2d_run, pwr1d_run = caller.execute(pairrun['map1'],
                                        pairrun['map2'],
                                        pairrun['noise_inv1'],
                                        pairrun['noise_inv2'],
                                        subtract_mean=subtract_mean,
                                        degrade_resolution=degrade_resolution,
                                        unitless=unitless,
                                        return_3d=return_3d,
                                        truncate=truncate,
                                        window=window, n_modes=n_modes,
                                        refinement=refinement,
                                        pad=pad, order=order)

            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run))
            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

        if (mode_num == 0):
            pwr_1d_from2d_zero_mode = copy.deepcopy(pwr_1d_from_2d)
            pwr_1d_zero_mode = copy.deepcopy(pwr_1d)
            pwr_2d_zero_mode = copy.deepcopy(pwr_2d)

        ttag = tag + "_%dmodes_2dtrans" % mode_num
        trans2d_mode = tf.calculate_2d_transfer_function(
                                            pwr_2d, pwr_2d_zero_mode,
                                            ttag)

        ttag = tag + "_%dmodes_1dtrans" % mode_num
        trans1d_mode = tf.calculate_1d_transfer_function(
                                            pwr_1d, pwr_1d_zero_mode,
                                            ttag)

        ttag = tag + "_%dmodes_1dtrans_from2d" % mode_num
        trans1d_mode = tf.calculate_1d_transfer_function(
                                            pwr_1d_from_2d, pwr_1d_from2d_zero_mode,
                                            ttag)

        transfer_functions[mode_num] = (trans1d_mode, trans2d_mode)

        mtag = tag + "_%dmodes" % mode_num

        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                          outdir="./plot_data")

    return transfer_functions
