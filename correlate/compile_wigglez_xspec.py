r"""code to assemble, calculate and write out GBT x WiggleZ"""
import numpy as np
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools
from optparse import OptionParser
from kiyopy import parse_ini
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
        file_tools.mkparents(outdir)

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
                              outdir=outdir,
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

            filename = "%s/%s_%s.dat" % (outdir,
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

def wrap_batch_gbtxwigglez_data_run(inifile, generate=False,
                                    outdir="./plots/"):
    r"""Wrapper to the GBT x WiggleZ calculation"""
    params_init = {"gbt_mapkey": "cleaned GBT map",
                   "wigglez_deltakey": "WiggleZ overdensity map",
                   "wigglez_mockkey": "WiggleZ overdensities from mocks",
                   "wigglez_selectionkey": "WiggleZ selection function",
                   "mode_transfer_1d_ini": "ini file -> 1d trans. function",
                   "mode_transfer_2d_ini": "ini file -> 2d trans. function",
                   "beam_transfer_ini": "ini file -> 2d beam trans. function",
                   "spec_ini": "ini file for the spectral estimation",
                   "output_tag": "tag identifying the output somehow"}
    prefix="cwx_"

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    print params

    output_tag = "%s_%s" % (params['gbt_mapkey'], params['output_tag'])
    output_root = "%s/%s/" % (outdir, output_tag)

    if generate:
        output_tag = None

    print output_root
    print output_tag
    file_tools.mkparents(output_root)
    parse_ini.write_params(params, output_root + 'params.ini',
                           prefix=prefix)

    datapath_db = data_paths.DataPath()

    batch_gbtxwigglez_data_run(params["gbt_mapkey"],
                               params["wigglez_deltakey"],
                               params["wigglez_mockkey"],
                               params["wigglez_selectionkey"],
                               inifile=params["spec_ini"],
                               datapath_db=datapath_db,
                               outdir=output_root,
                               output_tag=output_tag,
                               beam_transfer=None,
                               mode_transfer_1d=None,
                               mode_transfer_2d=None,
                               theory_curve=None)

if __name__ == '__main__':
    parser = OptionParser(usage="usage: %prog [options] filename (-h for help)",
                          version="%prog 1.0")
    parser.add_option("-g", "--generate", action="store_true",
                      dest="generate", default=False,
                      help="regenerate the cache of quadratic products")
    parser.add_option("-o", "--outdir", action="store",
                      dest="outdir", default="./plots/",
                      help="directory to write output data to")
    (optparam, inifile) = parser.parse_args()
    optparam = vars(optparam)

    if len(inifile) != 1:
        parser.error("inifile not specified")

    inifile = inifile[0]
    print optparam

    wrap_batch_gbtxwigglez_data_run(inifile,
                                    generate=optparam['generate'],
                                    outdir=optparam['outdir'])

