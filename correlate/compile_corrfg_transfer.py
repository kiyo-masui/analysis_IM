r"""code to assemble, calculate and write a transfer function

For WiggleZ crosspower:
(cleaned(realdata + Tbeamsim) x deltasim) / (Tbeamsim x deltasim)
deltasim is the original simulation in overdensity and Tbeamsim is the sim
convolved by the beam and multiplied by T_b. This simulates the wiggleZ
crosspower. The transfer function still needs a beam. Calculate the 2D beam
part of the WiggleZ tranfer function using
(Tbeamsim x deltasim) / (Tsim x deltasim)

You can also do the whole transfer function in one swoop as:
(cleaned(realdata + Tbeamsim) x deltasim) / (Tsim x deltasim)

do we want to base this on the auto or cross power (combined?)
"""
import numpy as np
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools

# TODO:
# rename dataplus sim key
# is deltasim_weightkey the wigglez selection function
# id basesim_weightkey the right weight thing

def batch_corrfg_transfer_run(dataplussim_key, delta_simkey,
                              deltasim_weightkey,
                              basesim_key, basesim_weightkey,
                              inifile=None, datapath_db=None,
                              outdir="./plots",
                              output_tag=None):
    r"""take relevant cross-powers
    dataplussim_key(map) * dataplussim_key(weight) x
    deltasim_weightkey * delta_simkey
    divided by:
    basesim_key * basesim_weightkey x deltasim_weightkey * delta_simkey
    """

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(dataplussim_key, "type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    if output_tag:
        output_root = "%s/%s/" % (outdir, output_tag)
        file_tools.mkparents(output_root)

    dbkeydict = {}
    dbkeydict['map1_key'] = basesim_key
    dbkeydict['map2_key'] = delta_simkey
    dbkeydict['noiseinv1_key'] = basesim_weightkey
    dbkeydict['noiseinv2_key'] = deltasim_weightkey
    files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                     datapath_db=datapath_db)

    base_pwrspec_out = caller.execute(files['map1_key'],
                                        files['map2_key'],
                                        files['noiseinv1_key'],
                                        files['noiseinv2_key'],
                                        inifile=inifile)

    if output_tag:
        basepwr_1d_from_2d = pe.convert_2d_to_1d(base_pwrspec_out[0])

        basepwr_1d = base_pwrspec_out[1]['binavg']
        basepwr_1d_from_2d = basepwr_1d_from_2d['binavg']

        bin_left = base_pwrspec_out[1]['bin_left']
        bin_center = base_pwrspec_out[1]['bin_center']
        bin_right = base_pwrspec_out[1]['bin_right']
        counts_histo = base_pwrspec_out[1]['counts_histo']

    transfer_functions = {}
    for treatment in map_cases['treatment']:
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:map;%s" % (dataplussim_key, treatment)
        dbkeydict['map2_key'] = delta_simkey
        dbkeydict['noiseinv1_key'] = "%s:weight;%s" % (dataplussim_key, treatment)
        dbkeydict['noiseinv2_key'] = deltasim_weightkey
        files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                      datapath_db=datapath_db)

        pwrspec_out_signal = caller.execute(files['map1_key'],
                                            files['map2_key'],
                                            files['noiseinv1_key'],
                                            files['noiseinv2_key'],
                                            inifile=inifile)

        if output_tag:
            pwr_1d_from_2d = pe.convert_2d_to_1d(pwrspec_out_signal[0])

            pwr_1d = pwrspec_out_signal[1]['binavg']
            pwr_1d_from_2d = pwr_1d_from_2d['binavg']

            trans1d_mode = pwr_1d / basepwr_1d
            trans1d_from2d_mode = pwr_1d_from_2d / basepwr_1d_from_2d

            trans2d_mode = pwrspec_out_signal[0]['binavg'] / \
                           base_pwrspec_out[0]['binavg']

            transfer_functions[treatment] = (trans1d_mode,
                                             trans1d_from2d_mode,
                                             trans2d_mode)

            filename = "%s/%s_%s.dat" % (output_root, output_tag, treatment)

            outfile = open(filename, "w")
            for specdata in zip(bin_left, bin_center,
                                bin_right, counts_histo,
                                pwr_1d, pwr_1d_from_2d,
                                trans1d_mode, trans1d_from2d_mode):
                outfile.write(("%10.15g " * 8 + "\n") % specdata)
            outfile.close()

    if not output_tag:
        caller.multiprocess_stack()
        return None
    else:
        return transfer_functions
