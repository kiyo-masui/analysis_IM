r"""code to assemble, calculate and write out GBT x WiggleZ"""
import numpy as np
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools
# TODO: figure out theory curve stuff


def batch_gbtxwigglez_data_run(gbt_map_key, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=None, datapath_db=None,
                               outdir="./plots",
                               output_tag=None,
                               beam_transfer=None,
                               mode_transfer_1d=None,
                               mode_transfer_2d=None,
                               theory_curve=None):
    r"""assemble the pairs of GBT and WiggleZ and calculate the cross-power"""

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(gbt_map_key, "type;treatment")
    mock_cases = datapath_db.fileset_cases(wigglez_mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    if output_tag:
        output_root = "%s/%s/" % (outdir, output_tag)
        file_tools.mkparents(output_root)

    for treatment in map_cases['treatment']:
        # TODO: make this more elegant
        # TODO: convert treatment into mode num
        transfer_2d = None
        if (mode_transfer_2d is not None) and (beam_transfer is None):
            transfer_2d = mode_transfer_2d[treatment][2]

        if (mode_transfer_2d is None) and (beam_transfer is not None):
            transfer_2d = beam_transfer

        if (mode_transfer_2d is not None) and (beam_transfer is not None):
            transfer_2d = mode_transfer_2d[treatment][2] * beam_transfer

        pwr_1d = []
        pwr_2d = []
        pwr_1d_from_2d = []
        for index in mock_cases['realization']:
            dbkeydict = {}
            dbkeydict['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
            dbkeydict['map2_key'] = "%s:%s" % (wigglez_mock_key, index)
            dbkeydict['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
            dbkeydict['noiseinv2_key'] = wigglez_selection_key
            files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                      datapath_db=datapath_db)

            pwrspec_out = caller.execute(files['map1_key'], files['map2_key'],
                                         files['noiseinv1_key'],
                                         files['noiseinv2_key'],
                                         inifile=inifile)

            if output_tag:
                pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                      transfer=transfer_2d))

                pwr_2d.append(pwrspec_out[0])
                pwr_1d.append(pwrspec_out[1])

        if output_tag:
            if mode_transfer_1d is not None:
                transfunc = mode_transfer_1d[treatment][1]
            else:
                transfunc = None

            mtag = output_tag + "_%s_mock" % treatment
            mean1dmock, std1dmock, covmock = pe.summarize_pwrspec(pwr_1d,
                              pwr_1d_from_2d, pwr_2d, mtag,
                              outdir=output_root,
                              apply_1d_transfer=transfunc)

        # now recover the xspec with the real data
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
        dbkeydict['map2_key'] = wigglez_map_key
        dbkeydict['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
        dbkeydict['noiseinv2_key'] = wigglez_selection_key
        files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                      datapath_db=datapath_db)

        pwrspec_out_signal = caller.execute(files['map1_key'],
                                            files['map2_key'],
                                            files['noiseinv1_key'],
                                            files['noiseinv2_key'],
                                            inifile=inifile)

        if output_tag:
            pwr_1d_from_2d = pe.convert_2d_to_1d(pwrspec_out_signal[0],
                                                 transfer=transfer_2d)

            pwr_1d = pwrspec_out_signal[1]['binavg']
            pwr_1d_from_2d = pwr_1d_from_2d['binavg']
            if mode_transfer_1d is not None:
                pwr_1d /= mode_transfer_1d[treatment][1]
                pwr_1d_from_2d /= mode_transfer_1d[treatment][1]

            # assume that they all have the same binning
            bin_left = pwrspec_out_signal[1]['bin_left']
            bin_center = pwrspec_out_signal[1]['bin_center']
            bin_right = pwrspec_out_signal[1]['bin_right']
            counts_histo = pwrspec_out_signal[1]['counts_histo']

            filename = "%s/%s_%s.dat" % (output_root,
                                         output_tag,
                                         treatment)

            outfile = open(filename, "w")
            for specdata in zip(bin_left, bin_center,
                                bin_right, counts_histo,
                                pwr_1d, pwr_1d_from_2d,
                                mean1dmock, std1dmock):
                outfile.write(("%10.15g " * 8 + "\n") % specdata)
            outfile.close()

            if theory_curve is not None:
                restrict = np.where(np.logical_and(bin_center > 0.09,
                                                   bin_center < 1.1))
                res_slice = slice(min(restrict[0]), max(restrict[0]))

                #restrict_alt = np.where(restrict)[0][np.newaxis, :]
                #restricted_cov = covmock[restrict_alt][0]

                #from core import utils
                amplitude = utils.ampfit(pwr_1d_from_2d[res_slice],
                                         covmock[res_slice, res_slice],
                                         theory_curve[res_slice])
                print "AMP:", mtag, treatment, amplitude

    if not output_tag:
        caller.multiprocess_stack()

    return None


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
