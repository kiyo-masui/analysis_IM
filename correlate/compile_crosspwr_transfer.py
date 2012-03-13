r"""Code to find the GBT x WiggleZ and one-sided transfer functions.
-------------------------------------------------------------------------------
For WiggleZ crosspower based which includes the anticorrelated cleaned fg:
(cleaned(realdata + Tbeamsim) x deltasim) / (Tbeamsim x deltasim)
cleaned_simkey points to cleaned(realdata + Tbeamsim)
truesignal_simkey points to the overdensity simulation (deltasim)
truesignal_weightkey points to the WiggleZ selection function
reference_simkey points to either the 21cm temperature sim w,w/o beam
reference_weightkey points to the weighting used in cleaned_simkey

deltasim is the original simulation in overdensity and Tbeamsim is the sim
convolved by the beam and multiplied by T_b. The beam transfer function needs
to be calculated separately as:
(Tbeamsim x deltasim) / (Tsim x deltasim)

The beam can be absorbed in this transfer function by calculating:
(cleaned(realdata + Tbeamsim) x deltasim) / (Tsim x deltasim)
-------------------------------------------------------------------------------
For WiggleZ crosspower based on just the cleaned sims:
(cleaned(Tbeamsim) x deltasim) / (Tbeamsim x deltasim)
cleaned_simkey points to cleaned(Tbeamsim)
and everything else is as-above
-------------------------------------------------------------------------------
For a one-sided transfer function to be used in the power spectrum
(cleaned(Tbeamsim) x Tbeamsim) / (Tbeamsim x Tbeansim)
cleaned_simkey points to cleaned(Tbeamsim)
truesignal_simkey points to either the 21cm temperature sim w,w/o beam
truesignal_weightkey is the ones matrix (uniform weight)
reference_simkey points to either the 21cm temperature sim w,w/o beam
reference_weightkey points to the weighting used in cleaned_simkey
"""
import numpy as np
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools


def batch_crosspwr_transfer(cleaned_simkey,
                            truesignal_simkey, truesignal_weightkey,
                            reference_simkey, reference_weightkey,
                            inifile=None, datapath_db=None,
                            outdir="./plots",
                            output_tag=None):
    r"""take relevant cross-powers
    cleaned_simkey(map) * cleaned_simkey(weight) x
    truesignal_weightkey * truesignal_simkey
    divided by:
    reference_simkey * reference_weightkey x truesignal_weightkey * truesignal_simkey
    """

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(cleaned_simkey, "type;treatment")

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
    dbkeydict['map1_key'] = reference_simkey
    dbkeydict['map2_key'] = truesignal_simkey
    dbkeydict['noiseinv1_key'] = reference_weightkey
    dbkeydict['noiseinv2_key'] = truesignal_weightkey
    files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                     datapath_db=datapath_db)

    reference_pwrspec_out = caller.execute(files['map1_key'],
                                        files['map2_key'],
                                        files['noiseinv1_key'],
                                        files['noiseinv2_key'],
                                        inifile=inifile)

    if output_tag:
        ref_pwr_1d_from_2d = pe.convert_2d_to_1d(reference_pwrspec_out[0])

        ref_pwr_1d = reference_pwrspec_out[1]['binavg']
        ref_pwr_2d = reference_pwrspec_out[0]['binavg']
        ref_pwr_1d_from_2d = ref_pwr_1d_from_2d['binavg']

        bin_left = reference_pwrspec_out[1]['bin_left']
        bin_center = reference_pwrspec_out[1]['bin_center']
        bin_right = reference_pwrspec_out[1]['bin_right']
        counts_histo = reference_pwrspec_out[1]['counts_histo']

    transfer_functions = {}
    for treatment in map_cases['treatment']:
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:map;%s" % (cleaned_simkey, treatment)
        dbkeydict['map2_key'] = truesignal_simkey
        dbkeydict['noiseinv1_key'] = "%s:weight;%s" % (cleaned_simkey, treatment)
        dbkeydict['noiseinv2_key'] = truesignal_weightkey
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
            pwr_2d = pwrspec_out_signal[0]['binavg']
            pwr_1d_from_2d = pwr_1d_from_2d['binavg']

            trans1d_mode = pwr_1d / ref_pwr_1d
            trans1d_from2d_mode = pwr_1d_from_2d / ref_pwr_1d_from_2d
            trans2d_mode = pwr_2d / ref_pwr_2d

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
